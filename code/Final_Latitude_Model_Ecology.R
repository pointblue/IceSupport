library(geosphere)
library(raster)
library(lubridate)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(sciplot)
library(performance)

DataDirectory<-"Directory where dataframe is saved/"
load(paste0(DataDirectory,"LatitudeDF_latIceDF.RData"))



#Load latitude model data frame "latIceDF" with attributes:
#"year" = Tag year as factor        
#"lonAvgBinom180" = Migration route east or west of 180 deg long 
#"month" = month of location date         
#"Lat"   = mean latitude across all individuals for a given migration type, month, and year         
#"meanicespd" = Average monthly ice speed from previous month    
#"icespdcurrent" = Average monthly ice speed from current month


#####Set plotting parameters
options(scipen=10000)
formatter1000 <- function(){
  function(x)x/1000
}

My_Theme = theme(
  axis.title.x = element_text(size = 25, colour = "black"),
  axis.text.x = element_text(size = 25, colour = "black"),
  axis.title.y = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"), title  = element_text(size = 30))


#Run linear model and evaluate
modlm<-lm(Lat~ meanicespd  + lonAvgBinom180 + year + month + icespdcurrent  , data = latIceDF )
summary(modlm)
r2(modlm)
effectsmixmodtop<-effects::allEffects(modlm)
plot(effectsmixmodtop)



check_model(modlm)



ggemmeans(modlm, terms= c('meanicespd'), back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 1)  + labs( x = "Ice speed",  y = "Latitude", title = "A")



ggemmeans(modlm, terms= c('meanicespd'), back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 1)  + labs( x = "Ice speed",  y = "Latitude", title = "(c)")  + My_Theme + scale_y_continuous(limits = c(-74,-60))


ggemmeans(modlm, terms= c('lonAvgBinom180'), back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 5)  + labs( x = "Migration route",  y = "Latitude", title = "(d)")  + My_Theme + My_Theme + scale_y_continuous(limits = c(-74,-64))



ggemmeans(modlm, terms= c('month'), back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 4)  + labs( x = "Month",  y = "Latitude", title = "(d)")  + My_Theme + My_Theme + scale_y_continuous(limits = c(-74,-60))



ggemmeans(modlm, terms= c('year'), back.transform = FALSE)%>%
  plot(add.data = TRUE, dot.size = 1)  + labs( x = "Year",  y = "Latitude", title = "D")  + My_Theme + My_Theme + scale_y_continuous(limits = c(-74,-64))







