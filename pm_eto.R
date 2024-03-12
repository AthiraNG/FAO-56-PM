
setwd()

library(sp)
library(raster)
library(lubridate)


# Inputs -----------------------------------------------------------

Tmax=raster("temp_max.tif") # Maximum temperature
Tmin=raster("temp_min.tif") # Minimum temperature
wind_speed=raster("wind_speed.tif") # wind speed  
u=raster("u_component_of_wind_10m.tif") # v component of wind
v=raster("v_component_of_wind_10m.tif") # u component of wind
Tdew=raster("dewpoint_temperature_2m.tif") # Dewpoint temperature
RHmax=raster("RH_max.tif")
RHmin=raster("RH_min.tif")
alb=raster("albedo.tif") # Albedo
z=raster("dem.tif") # Elevation above sea level in meter
latitude = raster("latitude.tif") # raster image containing latitude values
date="2020-01-01" # Date 
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
  return(decimal)
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
  return(ra)
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

calc_ETo = function(Rn,satvpr_slope,r,u2,Tmean,Rso,es,ea) {
  ETo=((0.409*Rn*satvpr_slope)/(satvpr_slope+(r*(1+(0.34*u2)))))+
      ((900*(u2/(Tmean+273)))*((es-ea)*r)*(satvpr_slope+
      (r*(1+(0.34*u2)))))
  return(ETo)
}



###Calculating all the variables using the functions###------------------


# Date of year (DOY)

j=yday(date)

# Radians 

if (file.exists(latitude)) {
  latitude_raster = raster(latitude)
  
  latitude_radians = (pi / 180) * latitude_raster
} else {
  
  decimal_degrees = dms_to_decimal(degrees, minutes, seconds)
  
  latitude_radians = (pi / 180) * decimal_degrees
}

# Mean daily temperature

Tmean=(Tmax+Tmin)/2

# Slope of saturation vapor pressure curve

satvpr_slope=calc_satvpr_slope(Tmean)

#Atmospheric Pressure

calc_atm_pressure= calc_atm_pressure(z)

# Psychrometric constant

r=0.000665*p

# Saturation vapor pressure at the air temperature

eTmax=calc_eTmax(Tmax)
eTmin=calc_eTmin(Tmin)

# The mean saturation vapor pressure

es=(eTmax+eTmin)/2

# Actual vapor pressure 

if (file.exists(Tdew)) {
  Tdew = raster(Tdew)
  ea = calc_ea(Tdew)
} else if (file.exists(RHmax) && file.exists(RHmin) && file.exists(Tmin)
  && file.exists(Tmax)) {
  RHmax = raster(RHmax)
  RHmin = raster(RHmin)
  Tmax = raster(Tmax)
  Tmin = raster(Tmin)
  # Calculate actual vapor pressure from RHmax and RHmin, if available
  ea <- calc_u2(RHmax,Rmin,Tmax,Tmin)
}


# The inverse relative distance Earth-Sun (dr)

dr=calc_dr(j)

# Solar declination

sd=calc_sd(j)

# Sunset hour angle

ωs=calc_ωs(radians,sd)

# Extraterrestrial radiation

Ra=calc_Ra(gsc, dr, ωs, radians, sd)

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

if (file.exists(wind_speed)) {
  WS <- raster(wind_speed)
  u2 <- calc_u2(WS, h)
} else if (file.exists(u) && file.exists(v)) {
  u <- raster(u)
  v <- raster(v)
  # Calculate wind speed from u and v components, if available
  WS <- sqrt(u^2 + v^2)
  u2 <- calc_u2(WS, h)
}


#Penman-Monteith Reference Evapotranspiration..............................

ETo=calc_ETo(Rn,satvpr_slope,r,u2,Tmean,Rso,es,ea)



