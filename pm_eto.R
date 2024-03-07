
setwd()

library(sp)
library(raster)
library(lubridate)


# Raster layers -----------------------------------------------------------

tmax=raster("temp_max.tif") # Maximum temperature
tmin=raster("temp_min.tif") # Minimum temperature
wind_speed=raster("wind_speed.tif") # wind speed  
u=raster("u_component_of_wind_10m.tif") # v component of wind
v=raster("v_component_of_wind_10m.tif") # u component of wind
dew=raster("dewpoint_temperature_2m.tif") # Dewpoint temperature
alb=raster("albedo.tif") # Albedo
z=raster("dem.tif") # Elevation above sea level in meter
latitude = raster("latitude.tif") # raster image containing latitude values


# Other inputs ------------------------------------------------------------

date="2020-01-01" # Date 
gsc= 0.0820 #solar constant MJ/m2
h=10 # Height of the wind measurement above the ground surface in meter

# function to convert DMS to degree decimal

dms_to_decimal = function(degrees, minutes, seconds) {
decimal = degrees + (minutes / 60) + (seconds / 3600)
return(decimal)
}

# function to calculate saturation vapour pressure curve

calc_satvpr_slope <- function(tmean) {
  d=(4098*(0.6108*(2.718282*(((17.27*tmean)/(tmean+237.3))))))/((tmean+237.2)^2)
  return(decimal)
}

# function to calculate saturation vapor pressure at the air temperature

calc_etmax <- function(tmax) {
  etmax=0.6108*(2.718282*((17.27*tmax)/(tmax+237.3)))
  return(etmax)
}
calc_etmin <- function(tmin) {
  etmin=0.6108*(2.718282*((17.27*tmin)/(tmin+237.3)))
  return(etmin)
}


calc_ea <- function(tmin) {
  etmin=0.6108*(2.718282*((17.27*tmin)/(tmin+237.3)))
  return(ea)
}

# function to calculate atmospheric pressure 

calc_atm_pressure <- function(z) {
  atm_pressure= (((293-(0.0065*z))/293)^5.26)*101.3
  return(atm_pressure)
}

# Function to calculate wind speed at 2 meters height (u2)
calc_u2 =  function(ws, h) {
  u2 = ws * (4.87 / log((67.8 * h) - 5.42))
  return(u2)
}

# Function to calculate extraterrestrial radiation (Ra)
calc_Ra = function(gsc, dr, ws, radians, sd) {
  Ra = (1440/pi)*(gsc)*(dr)*((ws*(sin(radians))*
       (sin(sd)))+((cos(radians))*(cos(sd))*(sin(ws))))
  return(ra)
}

# Function to calculate extraterrestrial radiation (Rnl)
calc_Rnl = function(tmax,tmin,ws,ea,rs,rso) {
  Rnl=(((tmax+273.16)^4+(tmin+273.16)^4)/2)*4.903*10^(-9)*
      (0.34-(0.14*sqrt(ea)))*((1.35*(rs/rso))-0.35)
  return(Rnl)
}

# Function to calculate FAO-56 Reference evapotranspiration (ETo)
calc_ETo = function(rn,satvpr_slope,r,u2,tmean,rso,es,ea) {
  ETo=((0.409*rn*satvpr_slope)/(satvpr_slope+(r*(1+(0.34*u2)))))+
      ((900*(u2/(tmean+273)))*((es-ea)*r)*(satvpr_slope+
      (r*(1+(0.34*u2)))))
  return(ETo)
}

#Date of year (DOY)

j=yday(date)

#Radians 

if (file.exists(latitude)) {
  latitude_raster = raster(latitude)
  
  latitude_radians = (pi / 180) * latitude_raster
} else {
  
  decimal_degrees = dms_to_decimal(degrees, minutes, seconds)
  
  latitude_radians = (pi / 180) * decimal_degrees
}

# Mean daily temperature

tmean=(tmax+tmin)/2

# Slope of saturation vapor pressure curve

satvpr_slope=calc_satvpr_slope(tmean)

#Atmospheric Pressure

calc_atm_pressure= calc_atm_pressure(z)

# Psychrometric constant

r=0.000665*p

#saturation vapor pressure at the air temperature

etmax=calc_etmax(tmax)

etmin=calc_etmin(tmin)

#The mean saturation vapor pressure

es=(etmax+etmin)/2

#Actual vapor pressure

ea=0.6108*(2.718282*((17.27*tmin)/(tmin+237.3)))

#The inverse relative distance Earth-Sun (dr)

dr=1+(0.033*(cos(((2*pi)/365)*j)))

# solar declination

sd=0.409*sin(((2*pi)/365)-1.39)

#Sunset hour angle

ws=acos((-tan(radians))*tan(sd))


#Extraterrestrial radiation

Ra=calc_Ra(gsc, dr, ws, radians, sd)

#Mean daily solar radiation

rs=(tmax-tmin)^0.5*Ra*0.16

#Clear sky solar radiation

rso=0.759*Ra

#Net solar or net shortwave radiation

rns=(1-alb)*rs

#Incoming net long wave radiation

Rnl=calc_Rnl(tmax,tmin,ws,ea,rs,rso)

#Net radiation

rn=rns-Rnl

#Average wind speed at 2 meter
if (file.exists(wind_speed)) {
  ws <- raster(wind_speed)
  u2 <- calc_u2(ws, h)
} else if (file.exists(u) && file.exists(v)) {
  u <- raster(u)
  v <- raster(v)
  # Calculate wind speed from u and v components, if available
  ws <- sqrt(u^2 + v^2)
  u2 <- calc_u2(ws, h)
}


#.........PM ET_o...............................................................

eto=calc_ETo(rn,satvpr_slope,r,u2,tmean,rso,es,ea)



