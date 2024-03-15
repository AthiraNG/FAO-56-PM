
setwd("C:\\Users\\Athira\\Downloads\\pm-20240314T121210Z-001\\pm")

library(sp)
library(raster)
library(lubridate)


# Inputs -----------------------------------------------------------

Tmin=raster("tmin2019-1.tif")# Maximum temperature
Tmax=raster("tmax2019-1.tif")# Minimum temperature
stack_tminmax=stack(Tmax,Tmin)# Stack of min and max temperature
wind_speed=raster("uh2019-1.tif") # wind speed  
u=raster("u_component_of_wind_10m.tif") # v component of wind
v=raster("v_component_of_wind_10m.tif") # u component of wind
Tdew=raster("tdew2019-1.tif") # Dewpoint temperature
RHmax=raster("RH_max.tif")# Maximum relative humidity
RHmin=raster("RH_min.tif")# Minimun relative humidity
alb=raster("alb2019-1.tif")# Albedo
z=raster("dem.tif") # Elevation above sea level in meter
latitude = raster("rad2019-1.tif") # raster image containing latitude values
date="2019-01-01" # Date 
gsc= 0.0820 #solar constant MJ/m2
h=10 # Height of the wind measurement above the ground surface in meter


# Functions to calculate all the variables --------------------------------


# function to convert DMS to degree decimal

dms_to_decimal = function(degrees, minutes, seconds) {
decimal = degrees + (minutes / 60) + (seconds / 3600)
return(decimal)
}

# function to calculate saturation vapour pressure curve

calc_satvpr_slope <- function(Tmean) {
  d=(4098*(0.6108*(2.718282*(((17.27*Tmean)/(Tmean+237.3))))))/((Tmean+237.2)^2)
  return(d)
}

# function to calculate saturation vapor pressure at the air temperature

calc_eTmax <- function(Tmax) {
  eTmax=0.6108*(2.718282*((17.27*Tmax)/(Tmax+237.3)))
  return(eTmax)
}
calc_eTmin <- function(Tmin) {
  eTmin=0.6108*(2.718282*((17.27*Tmin)/(Tmin+237.3)))
  return(eTmin)
}

# function to calculate actual vapor pressure

calc_ea <- function(Tdew = NULL, RHmax = NULL, RHmin = NULL, eTmin = NULL, eTmax = NULL) {
  if (!is.null(Tdew)) {
    # Calculate ea using Tdew
    ea = 0.6108 * (2.718282 * ((17.27 * Tdew) / (Tdew + 237.3)))
  } else if (!is.null(RHmax) && !is.null(RHmin) && !is.null(eTmin) && !is.null(eTmax)) {
    # Calculate ea using RHmax, RHmin, eTmin, and eTmax
    ea = ((eTmin * (RHmax / 100)) + (eTmax * (RHmin / 100))) / 2
  }
  
  return(ea)
}

# function to calculate inverse relative distance Earth-Sun (dr)

calc_dr <- function(j) {
  dr=1+(0.033*(cos(((2*pi)/365)*j)))
  return(dr)
}

# function to calculate solar declination 

calc_sd <- function(j) {
  sd=0.409*sin((((2*pi)/365)*j)-1.39)
  return(sd)
}

# function to calculate solar hour angle
calc_ωs <- function(radians,sd) {
  ωs=acos((-tan(radians))*tan(sd))
  return(ωs)
}

# function to calculate atmospheric pressure 

calc_atm_pressure <- function(z) {
  atm_pressure= (((293-(0.0065*z))/293)^5.26)*101.3
  return(atm_pressure)
}

# Function to calculate wind speed at 2 meters height (u2)

calc_u2 =  function(ωs, h) {
  u2 = ωs * (4.87 / log((67.8 * h) - 5.42))
  return(u2)
}

# Function to calculate extraterrestrial radiation (Ra)

calc_Ra = function(gsc, dr, ωs, radians, sd) {
  Ra = (1440/pi)*(gsc)*(dr)*((ωs*(sin(radians))*
       (sin(sd)))+((cos(radians))*(cos(sd))*(sin(ωs))))
  return(Ra)
}

# Function to calculate Mean daily solar radiation (Rs)

calc_Rs <- function(Tmax,Tmin,Ra) {
  Rs=(Tmax-Tmin)^0.5*Ra*0.16
  return(Rs)
}
# Function to calculate Net solar or net shortwave radiation

calc_Rns <- function(alb,Rs) {
  Rns=(1-alb)*Rs
  return(Rns)
}

# Function to calculate extraterrestrial radiation (Rnl)

calc_Rnl = function(Tmax,Tmin,ωs,ea,Rs,Rso) {
  Rnl=(((Tmax+273.16)^4+(Tmin+273.16)^4)/2)*4.903*10^(-9)*
      (0.34-(0.14*sqrt(ea)))*((1.35*(Rs/Rso))-0.35)
  return(Rnl)
}

# Function to calculate FAO-56 Reference evapotranspiration (ETo)

calc_ETo = function(Rn,d,r,u2,Tmean,Rso,es,ea) {
  ETo=((0.409*Rn*d)/(d+(r*(1+(0.34*u2)))))+
      ((900*(u2/(Tmean+273)))*((es-ea)*r)*(d+
      (r*(1+(0.34*u2)))))
  return(ETo)
}



###Calculating all the variables using the functions###------------------


# Date of year (DOY)

j=yday(date)

# Radians 


if (file.exists("rad2019-1.tif")) {
  latitude <- raster("rad2019-1.tif")
  
  latitude_radians <- (pi / 180) * latitude
} else {
  # Assuming you have defined 'degrees', 'minutes', and 'seconds'
  decimal_degrees <- dms_to_decimal(degrees, minutes, seconds)
  
  latitude_radians <- (pi / 180) * decimal_degrees
}


# Mean daily temperature

Tmean=calc(stack_tminmax,mean)


# Slope of saturation vapor pressure curve

d=calc_satvpr_slope(Tmean)

#Atmospheric Pressure

atm_pressure= calc_atm_pressure(z)


# Psychrometric constant

r=0.000665*atm_pressure

# Saturation vapor pressure at the air temperature

eTmax=calc_eTmax(Tmax)
eTmin=calc_eTmin(Tmin)

# The mean saturation vapor pressure
eTmax_Tmin_stcak=stack(eTmax,eTmin)
es=calc(eTmax_Tmin_stcak,mean)

# Actual vapor pressure 

if (file.exists("tdew2019-1.tif")) {
  Tdew = raster("tdew2019-1.tif")
  ea = calc_ea(Tdew)
} else if (file.exists("RH_max.tif") && file.exists("RH_min.tif") && file.exists("tmin2019-1.tif") 
  && file.exists("tmax2019-1.tif")) {
  RHmax = raster("RH_max.tif")
  RHmin = raster("RH_min.tif")
  Tmax = raster("tmax2019-1.tif")
  Tmin = raster("tmin2019-1.tif")
  # Calculate actual vapor pressure from RHmax and RHmin, if available
  ea <- calc_u2(RHmax,Rmin,Tmax,Tmin)
}


# The inverse relative distance Earth-Sun (dr)

dr=calc_dr(j)

# Solar declination

sd=calc_sd(j)

# Sunset hour angle

ωs=calc_ωs(latitude_radians,sd)

# Extraterrestrial radiation

Ra=calc_Ra(gsc, dr, ωs, latitude_radians, sd)

# Mean daily solar radiation

Rs=calc_Rs(Tmax,Tmin,Ra)

# Clear sky solar radiation

Rso=0.759*Ra

# Net solar or net shortwave radiation

Rns=calc_Rns(alb,Rs)

# Incoming net long wave radiation

Rnl=calc_Rnl(Tmax,Tmin,ωs,ea,Rs,Rso)

# Net radiation

Rn=Rns-Rnl

# Average wind speed at 2 meter

if (file.exists("uh2019-1.tif")) {
  WS <- raster("uh2019-1.tif")
  u2 <- calc_u2(WS, h)
} else if (file.exists(u) && file.exists(v)) {
  u <- raster(u)
  v <- raster(v)
  # Calculate wind speed from u and v components, if available
  WS <- sqrt(u^2 + v^2)
  u2 <- calc_u2(WS, h)
}


#Penman-Monteith Reference Evapotranspiration..............................

ETo=calc_ETo(Rn,d,r,u2,Tmean,Rso,es,ea)


