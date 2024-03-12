
library(sp)
library(raster)

#...................................................................................................................

wd=setwd("E:\\jan_2018")
lst=raster('lst_jan.tif')
ndvi=raster('ndvi_jan.tif')
alb_24=raster('alb_jan.tif')
lst
ndvi
alb_24

kc=exp((1.8)+((-0.008)*((lst-273.15)/(alb_24*ndvi))))

#.....................................................................................................................


values=read.csv("meter_values.csv")
tmax=mean(values$max)-273.15
tmin=mean(values$min)-273.15
u=mean(values$u_component_of_wind_10m)
v=mean(values$v_component_of_wind_10m)
dew=mean(values$dewpoint_temperature_2m)
alb=mean(values$alb)

z=349
j=1
gsc= 0.0820
as=0.25
bs=0.50

##latitude
#degrees=30
#minutes=7
#seconds=27.80
radians=38.85154

##wind component
h=10


##########################################################################################################################################

# Mean daily temperature

tmean=(tmax+tmin)/2

# Slope of saturation vapor pressure curve

d=(4098*(0.6108*(2.718282*(((17.27*tmean)/(tmean+237.3))))))/((tmean+237.2)^2)

#Atmospheric Pressure

p= (((293-(0.0065*z))/293)^5.26)*101.3

# Psychrometric constant

r=0.000665*p

#saturation vapor pressure at the air temperature

etmax=0.6108*(2.718282*((17.27*tmax)/(tmax+237.3)))

etmin=0.6108*(2.718282*((17.27*tmin)/(tmin+237.3)))

#The mean saturation vapor pressure

es=(etmax+etmin)/2

#Actual vapor pressure

ea=0.6108*(2.718282*((17.27*dew)/(dew+237.3)))

#The inverse relative distance Earth-Sun (dr)

dr=1+(0.033*(cos(((2*pi)/365)*j)))

# solar declination

sd=0.409*sin(((2*pi)/365)-1.39)

#radians

#Degree_decimals =degrees+(minutes/60)+(seconds/3600)
#radians=(pi/180)*(Degree_decimals)

#Sunset hour angle

ws=acos((-tan(radians))*tan(sd))

sin(radians)
sin(sd)
cos(radians)
cos(sd)
sin(ws)

#Extraterrestrial radiation

ra=(1440/pi)*(gsc)*(dr)*((ws*(sin(radians))*(sin(sd)))+((cos(radians))*(cos(sd))*(sin(ws))))

#Mean daily solar radiation

rs=(tmax-tmin)^0.5*ra*0.16

# maximum possible daily duration of sunshine in the same day (hour)

# N =- ws * (24 / pi)

# actual daily duration of sunshine (hour)

# n = tmean * 0.232 + 4.352

# rs <- ((n/N) * bs + as) * Ra

#Clear sky solar radiation

rso=0.759*ra

#Net solar or net shortwave radiation

rns=(1-alb)*rs

#Incoming net long wave radiation

rnl=(((tmax+273.16)^4+(tmin+273.16)^4)/2)*4.903*10^(-9)*(0.34-(0.14*sqrt(ea)))*((1.35*(rs/rso))-0.35)

#Net radiation

rn=rns-rnl

#Average wind speed at 2 meter

uh=sqrt((u)^2+(v)^2)
u2=uh*(4.87/log((67.8*h)-5.42))

#.........PM ET_o.......................................................................................................

eto=((0.409*rn*d)/(d+(ra*(1+(0.34*u2)))))+((900*(u2/(tmean+273)))*((es-ea)*r)*(d+(r*(1+(0.34*u2)))))
eto

#.......................................................................................................................

eta=kc*eto

