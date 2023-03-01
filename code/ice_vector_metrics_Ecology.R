#### Get location data and attribute with ice vector rasters and compute vector metrics - ice support, drift, ice and penguin speeds and directions
#Set 5 day bins and aggregate at 5 day bins
#Globals
library(ncdf4)
library(raster)
library(rgdal)
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(skimr)
library(ggsci)




############################
#Set functions
############################

#EPSG code 3409 this is the sea ice vector projection
ease<- "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs"
easecrs<-crs(ease)


#This (EPSG 3976) is the geodetically preferred alternative to NSIDC polar sterographic projection with the hughes ellipsoid (EPSG 3412), replacing it with with WGS84
pss84<-"+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"



#Angle functions
deg2rad = pi/180  # Converts degrees to radians when you multiply by the degrees (radius of a circle is 2pi*radius)
rad2deg = 180/pi

deg_r = function(radians) {
  180 * radians / pi
}


rad_r = function(degrees) {
  degrees * pi / 180
}

## ===========================================================================
## convert delta x and delta y to direction and length (or u & v to wd & ws)
## ===========================================================================
 #This calculates direction TOWARDS (i.e bearing)
calc_dir = function(x, y) {
  w = deg_r(atan2(y, x))
  w = ifelse(w <= 0, w + 360, w)
  dir = 90 - w
  dir = ifelse(w >= 90, dir + 360, dir)
  return(dir)
}

############################





####Caluclating the ice assistance for tracks####

############################
#Load files and set directories
############################



#Read in your location data
rawFile<-read.csv("V:/MyFolder/LocationData.csv")


##adjust LON to true unprojected values (i.e. negative once you cross the 180 line going East)
##Raw has continuous LON that goes beyond 180.  Lat is already in negative values
rawFile$LonFix <- ifelse(rawFile$Lon.50. > 180, rawFile$Lon.50. - 360, rawFile$Lon.50.)
rawFile$Lon97Fix <- ifelse(rawFile$Lon.97.5. > 180, rawFile$Lon.97.5. - 360, rawFile$Lon.97.5.)
rawFile$Lon25Fix <- ifelse(rawFile$Lon.2.5. > 180, rawFile$Lon.2.5. - 360, rawFile$Lon.2.5.)


#add attributes.  Keep "RLong" here for mapping in R
data2 <- data.frame(Long=rawFile$LonFix ,
                    Lat=rawFile$Lat.50.,
                    RLong=rawFile$Lon.50.,
                    RLat=(rawFile$Lat.50.)*-1,
                    Lon25=rawFile$Lon25Fix,
                    Lat25=rawFile$Lat.2.5.,
                    Lon97=rawFile$Lon97Fix,
                    Lat97=rawFile$Lat.97.5.,
                    mydate=rawFile$Time1,
                    day=day(rawFile$Time1),
                    month=month(rawFile$Time1),
                    year=year(rawFile$Time1),
                    doy=julian(as.Date(rawFile$Time1),origin = as.Date(yearorigin)))


#add simple ymd column
data2$plotdate<-ymd(paste(data2$year,data2$month,data2$day,sep="-"))




###Convert table to an SF spatial object
projectedxy <- st_as_sf(data2, coords = c("Long", "Lat"), remove = FALSE, crs = 4326, agr = "constant")


###Project to polar stereo wgs84
df.SP<-st_transform(x = projectedxy, crs = pss84)
###Add polar stereo wgs84 XY coords then export to table
df.SP$Xpsgwgs84<-st_coordinates(df.SP)[,1] # get coordinates
df.SP$Ypsgwgs84<-st_coordinates(df.SP)[,2] # get coordinates



###Project to EASE
df.E<-st_transform(x = df.SP, crs = ease)
###Add EASE XY coords then export to table
df.E$Xease<-st_coordinates(df.E)[,1] # get coordinates
df.E$Yease<-st_coordinates(df.E)[,2] # get coordinates


###convert to table
alltracks<-st_set_geometry(df.E, NULL)


#This creates a column with 5 day bins that can be used to aggregate by bird, monthly bin, and year.
#create a day column based on the date column
alltracks$mnthbin[alltracks$day >= 1 & alltracks$day <= 5]<-1
alltracks$mnthbin[alltracks$day >= 6 & alltracks$day <= 10]<-2
alltracks$mnthbin[alltracks$day >= 11 & alltracks$day <= 15]<-3
alltracks$mnthbin[alltracks$day >= 16 & alltracks$day <= 20]<-4
alltracks$mnthbin[alltracks$day >= 21 & alltracks$day <= 25]<-5
alltracks$mnthbin[alltracks$day >= 26 & alltracks$day <= 31]<-6
alltracks$plotBinBird<-paste0(alltracks$bird_fn,alltracks$month,alltracks$mnthbin,alltracks$year)
alltracks$mnthbinyear<-paste0(alltracks$month,".",alltracks$mnthbin)
alltracks$binorder<-as.numeric(alltracks$mnthbinyear)



##########################################################
###Attribute locs with uv ice data
##########################################################

#load u,v ice vector data for given year indicated by myyear
#This assumes you have downloaded the ice vector data from online and saved them as yearly netcdf files.
ncu<-raster::brick(paste0("C:/YourDirectory/icemotion_daily_sh_25km_",myyear,"0101_",myyear,"1231_v4.1.nc"), var="u")

#Define with EASE grid south nonwgs84
ncu<-projectRaster(ncu,crs=easecrs)


ncv<-raster::brick(paste0("C:/YourDirectory/icemotion_daily_sh_25km_",myyear,"0101_",myyear,"1231_v4.1.nc"), var="v")

ncv<-projectRaster(ncv,crs=easecrs)


ncerr<-raster::brick(paste0("C:/YourDirectory/icemotion_daily_sh_25km_",myyear,"0101_",myyear,"1231_v4.1.nc"), var="icemotion_error_estimate")

projection(ncerr)<-easecrs


#Extract ice vector values by date

for(dd in 1:nlayers(ncu)){

  ncutemp<-ncu[[dd]]

  ncvtemp<-ncv[[dd]]

  ncerrtemp<-ncerr[[dd]]


  nametemp<-ymd(stringr::str_sub(names(ncv[[dd]]),start =  2 ))

  datetrack<-alltracks%>%
    filter(plotdate == nametemp)


  #extract to subset of points then bind to full sp data frame using pss84
  u<-raster::extract(ncutemp,datetrack[c("Xease","Yease")],method="simple",sp=TRUE)
  datetrack$u<-u

  #extract to subset of points then bind to full sp data frame using pss84
  v<-raster::extract(ncvtemp,datetrack[c("Xease","Yease")],method="simple",sp=TRUE)
  datetrack$v<-v


  #extract to subset of points then bind to full sp data frame using pss84
  err<-raster::extract(ncerrtemp,datetrack[c("Xease","Yease")],method="simple",sp=TRUE)
  datetrack$uverr<-err




  # ##########################################
  # ##Bind extraction results together
  # ##########################################
  if(dd == 1){

    alldats<-datetrack

  }else{

    alldats<-rbind(alldats,datetrack)

  }

}


###########################################################################
####################### Calculate vector based ice support metrics
######################### Following equations within manuscript
###########################################################################


#list birds ids where bird_fn is unique per individual
  birdids<-unique(alldats$bird_fn)
  birdids<-na.omit(birdids)
  assistdf<-data.frame()

  #For each birdID
  for(bb in 1:length(birdids)){

    print(birdids[bb])

    track<-alldats %>%
      filter(bird_fn ==birdids[bb])

    print(table(track$month))


    #Get observed (ground) penguin BEARING between set of consecutive coordinates.  North is 0 south is 180
    track$obsbearing <- c(geosphere::bearingRhumb(track[-nrow(track),c("Long","Lat")], track[-1,c("Long","Lat")]), NA) #creates pair of df's with offset coords


    #Get observed (ground) penguin SPEED = distance (in meters) then divide by time (in seconds) between consecutive coordinates
    track$obsspeed   <- c(geosphere::distGeo(track[-nrow(track),c("Long","Lat")], track[-1,c("Long","Lat")]), NA)/c(apply(cbind(track$rdatetime[-nrow(track)], track$rdatetime[-1]), 1, diff),NA)

    #convert to uv components of penguin ground vector
    #http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    track$Xg <- track$obsspeed * sin(track$obsbearing*deg2rad)
    track$Yg <- track$obsspeed * cos(track$obsbearing*deg2rad)



    #Double check that this Dir matches obsbearing
   track$trackDir<- calc_dir(track$Xg, track$Yg)


    #get observed distance between successive coordinates
    track$distance   <- c(geosphere::distGeo(track[-nrow(track),c("Long","Lat")], track[-1,c("Long","Lat")]), NA)




    ###Run Ice numbers###########################################
    #Apply rotation matrix to convert uv to calculate ice speed and direction
    #https://nsidc.org/support/how/how-convert-horizontal-and-vertical-components-east-and-north
    track$Xi<-(track$u * cos(track$Long*deg2rad))  -  (track$v*sin(track$Long*deg2rad)) #aka easterly
    track$Yi<-(track$u * sin(track$Long*deg2rad))  +  (track$v*cos(track$Long*deg2rad)) #aka northerly


    #Calculate ice velocity
    track$icespd <- sqrt(track$Xi^2 + track$Yi^2)*0.01 #0.01 to convert from cm to *meters per second*


    #Calculate ice direction
    track$icedir <- calc_dir(track$Xi, track$Yi)




    ############Calculate support metrics ##################################

    #Theta is the difference in angle between the ice direction and the penguin observed bearing
    track$theta<-180 - abs(180 - abs(track$icedir  - track$obsbearing ) %% 360)


    #Calculate penguin speed in relation to ice movement (aka 'air speed')
    track$airspeed<-(sqrt(track$obsspeed^2 + track$icespd^2 - 2 * track$obsspeed * track$icespd * (cos(track$theta * deg2rad))))


    track$Xa<-track$obsspeed*sin(track$obsbearing*deg2rad) - track$icespd*sin(track$icedir*deg2rad)
    track$Ya<-track$obsspeed*cos(track$obsbearing*deg2rad) - track$icespd*cos(track$icedir*deg2rad)


    track$airdir <- calc_dir(track$Xa, track$Ya)




    #Calculate difference in angle between expected direction and observed direction of penguin
    track$driftangle <- asin((track$icespd*sin(track$theta * deg2rad)/ track$airspeed )) * rad2deg


    ######Now calculate scalar PROJECTION of ice onto penguin air speed aka actual speed or speed relative to       the ice
   #Ice uv components are Xi Yi (Vi) and penguin air speed components are Xa and Ya (Vp)

    track$Projia<-((track$Xa*(track$Xi*0.01)) + (track$Ya*(track$Yi*0.01)))/sqrt(track$Xa^2 + track$Ya^2)

    #bind results
    assistdf<-rbind(assistdf,track)




  }







  ###############
  #aggregate values at 5 day bins
  ################

  #run averages
  mean5day <- aggregate(list(support5day = assistdf$Projia, theta5day = assistdf$theta, drift5day = assistdf$driftangle, obsspeed5day = assistdf$obsspeed, icespd5day = assistdf$icespd, airspeed5day = assistdf$airspeed, sic5day = assistdf$sic, err5day = assistdf$d975errwght), by = list(plotBinBird = assistdf$plotBinBird), mean, na.rm = TRUE)

  #Run 5 day sum for distance
  dist5day <- aggregate(list(dist5day = assistdf$distance) , by = list(plotBinBird = assistdf$plotBinBird), sum, na.rm = TRUE)
  assistdf<-data.table(assistdf)


  #aggregates based on earliest date in the 5 day bin
  fivday <- aggregate(list(rdatetime = assistdf$rdatetime), list(plotBinBird = assistdf$plotBinBird), min)
  fivedaymerge<-merge(fivday,assistdf, by = c('plotBinBird','rdatetime'),all.x = TRUE,all.y=FALSE)
  fivedaymerge <- fivedaymerge[order(fivedaymerge$X),]





