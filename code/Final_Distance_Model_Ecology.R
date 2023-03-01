library(effects)
library(MuMIn)
library(ggplot2)
library(performance)
library(nlme)
library(dplyr)
library(sjPlot)
library(viridis)
library(corrplot)
library(ggeffects)
library(data.table)
library(ggeffects)
library(emmeans)
library(lubridate)

DataDirectory<-"Directory where dataframe is saved/"
load(paste0(DataDirectory,"IceSupportDF_envinddatana.RData"))

#Load distance model data frame "envinddatana" with attributes:
#"birdyear"  = Unique individual bird ID and year of tag
#"dist5day"  = Total distance traveled
#"support5day" = Ice support index
#"drift5day"   = Drift value
#"yearweek"    = Week number within a given year
#"year.f"      = Tag year as factor
#"sic5day"     = Average sea ice concentration
#"bird_fn"     = Unique individual bird ID
#"err5day"     = Location error
#"divesum5"    = Total dive time
#"lonAvgBinom180" = Migration route east or west of 180 deg long
#"age"           = Bird age in years
#"sex.f"         = Bird sex
#"Breeding.f"    = Bird breeding status
#"colony"        = Bird breeding colony

#####Set plotting parameters
options(scipen=10000)
formatter1000 <- function(){
  function(x)x/1000
}

My_Theme = theme(
  axis.title.x = element_text(size = 25, colour = "black"),
  axis.text.x = element_text(size = 25, colour = "black"),
  axis.title.y = element_text(size = 25, colour = "black"), axis.text.y = element_text(size = 25, colour = "black"), title  = element_text(size = 30))




#############################################




##############################
#########################Distance by ice support (aka Projection)
##############################################
## STEP 1
#mixed or non mixed assessment
##without model weights
##Using REML here as per Zuur et al for comparing mixed vs non mixed models


mixmod1<-lme(log(dist5day) ~ yearweek*year.f + sic5day + divesum5 + support5day*lonAvgBinom180 + drift5day*lonAvgBinom180 + age + sex.f + Breeding.f , random=~1|bird_fn, data = envinddatana, na.action=na.omit, method = "REML")
summary(mixmod1)


nonmixmod<-gls(log(dist5day) ~ yearweek*year.f + sic5day + divesum5 + support5day*lonAvgBinom180 + drift5day*lonAvgBinom180 + age + sex.f + Breeding.f ,data = envinddatana,na.action=na.omit , method = "REML")

anova(mixmod1, nonmixmod)



######################################################
###STEP 2
## Find optimal FIXED structure using an ML testing procedure
###In this case make sure to us ML for comparison of nested models
##Full model



mixmod<-lme(log(dist5day) ~ yearweek*year.f  + sic5day + divesum5 + support5day*lonAvgBinom180 + drift5day*lonAvgBinom180 + age + sex.f + Breeding.f, random=~1|bird_fn, data = envinddatana, weights = ~ I(1/err5day), na.action=na.fail, method = "ML")

plot(mixmod)


summary(mixmod)
r2(mixmod)
effectsmixmod<-effects::allEffects(mixmod)
plot(effectsmixmod)
intervals(mixmod,which = 'fixed')

#Selection procedure
selectmixmod <- MuMIn::dredge(mixmod, rank="AIC")
subset(selectmixmod, delta <2) # show the top models within 2 AIC
mysum<-summary(model.avg(selectmixmod, delta<2))


#Run and examine each model within 2 AIC
#top model
mixmod1<-lme(log(dist5day)  ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod1)
summary(mixmod1)
intervals(mixmod1,which = 'fixed')
effectsmixmod<-effects::allEffects(mixmod1)
plot(effectsmixmod)


#2 add age
mixmod2<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f + age, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod2)
summary(mixmod2)
intervals(mixmod2,which = 'fixed')



#3 add support lon interactin
mixmod3<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f + support5day:lonAvgBinom180, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod3)
summary(mixmod3)
intervals(mixmod3,which = 'fixed')
effectsmixmod<-effects::allEffects(mixmod3)
plot(effectsmixmod)


#4 add age and support long interaction
mixmod4<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f +  support5day:lonAvgBinom180 + age, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod4)
summary(mixmod4)
intervals(mixmod4,which = 'fixed')


#5 add drift lon interaction
mixmod5<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f + drift5day:lonAvgBinom180, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod5)
summary(mixmod5)
intervals(mixmod5,which = 'fixed')


#6 add age and drift lon interaction
mixmod6<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f + drift5day:lonAvgBinom180 + age, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmod6)
summary(mixmod6)
intervals(mixmod6,which = 'fixed')


#null model
mixmodNull<-lme(log(dist5day) ~1, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
r2(mixmodNull)

#Model comparison
anova(mixmod1, mixmod2, mixmod3, mixmod4, mixmod5, mixmod6, mixmodNull)
model.sel( mixmod1, mixmod2, mixmod3, mixmod4, mixmod5, mixmodNull)



######################################################
##STEP 3
#Run top model with REML

mixmodtop<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "REML")


summary(mixmodtop)
r2(mixmodtop)
effectsmixmodtop<-effects::allEffects(mixmodtop)
plot(effectsmixmodtop)
intervals(mixmodtop,which = 'fixed')

#Get imapct values for covars
ggpredict(mixmodtop, terms= c('divesum5'), back.transform = TRUE)
ggpredict(mixmodtop, terms= c('lonAvgBinom180'), back.transform = TRUE)
ggpredict(mixmodtop, terms= c('Breeding.f'), back.transform = TRUE)
ggpredict(mixmodtop, terms= c('year.f', 'yearweek'), back.transform = TRUE)
ggpredict(mixmodtop, terms= c('sic5day'), back.transform = TRUE)




###############PLOT results

png(filename="topModel_Distance_supporta.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "support5day", back.transform = TRUE)%>%
  plot(add.data = FALSE, dot.size = 0.1) + scale_y_continuous(limits = c(100000,400000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Ice support",  y = "Distance (km)", title = "(a)") + My_Theme 
dev.off()


png(filename="topModel_Distance_driftb.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "drift5day", back.transform = TRUE)%>%
  plot(add.data = FALSE, dot.size = 0.1) + scale_y_continuous( limits = c(100000,400000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Drift",  y = "Distance (km)", title = "(b)")  + My_Theme
dev.off()


png(filename="topModel_Distance_yearweekyearc.png", width=600, height=600)
ggemmeans(mixmodtop, terms= c('yearweek', 'year.f'), back.transform = TRUE)%>%
  plot() + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Week by Year",  y = "Distance (km)", title = "(c)") + My_Theme
dev.off()


png(filename="topModel_Distance_longAvgBinom180d.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "lonAvgBinom180", back.transform = TRUE)%>%
  plot() + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Migration type",  y = "Distance (km)", title = "(d)") + My_Theme
dev.off()


png(filename="topModel_Distance_divesum5d.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "divesum5", back.transform = TRUE)%>%
  plot(add.data = FALSE, dot.size = 0.1) + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Total Dive Time (hours)",  y = "Distance (km)", title = "(e)") + My_Theme
dev.off()



png(filename="topModel_Distance_sic5day.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "sic5day", back.transform = TRUE)%>%
  plot(add.data = TRUE, dot.size = 0.1) + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Sea ice concentration %",  y = "Distance (km)", title ="(f)") + My_Theme
dev.off()

png(filename="topModel_Distance_sex.f.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "sex.f", back.transform = TRUE)%>%
  plot() + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Sex",  y = "Distance (km)", title ="(g)") + My_Theme
dev.off()


png(filename="topModel_Distance_breeding.png", width=600, height=600)
ggemmeans(mixmodtop, terms= "Breeding.f", back.transform = TRUE)%>%
  plot() + scale_y_continuous(limits = c(100000,375000), breaks = c(100000,150000,200000,250000,300000,350000),labels=formatter1000()) + labs( x = "Breeding status",  y = "Distance (km)", title = "(h)") + My_Theme
dev.off()

######################################################
##STEP 4
#Model validation with REML model
#use lme4 lmer to be able to use check_model function
#Plotting of residuals

mixmodtoplme4<- lme4::lmer(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f + (1|bird_fn), weights = err5day, data = envinddatana, REML = TRUE, na.action = na.fail)
check_model(mixmodtoplme4) #produces plot with several graphs to check assumptions

summary(mixmodtoplme4)

mixmodtop<-lme(log(dist5day) ~Breeding.f + drift5day + lonAvgBinom180  + sic5day + divesum5 + support5day + year.f +  yearweek + yearweek:year.f + sex.f, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "REML")

plot(residuals(mixmodtop, type= "normalized"))
R1<-residuals(mixmodtop, type= "normalized")
F1<-fitted(mixmodtop)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals"
plot(x = F1, y = R1, xlab = "Fitted values", ylab = MyYlab)
boxplot(R1 ~ drift5day, data = envinddatana,
        main = "Drift ", ylab = MyYlab)

boxplot(R1 ~ support5day, data = envinddatana,
        main = "Support ", ylab = MyYlab)

boxplot(R1 ~ year.f, data = envinddatana,
        main = "Year ", ylab = MyYlab)

boxplot(R1 ~ Breeding.f, data = envinddatana,
        main = "Breeding ", ylab = MyYlab)

boxplot(R1 ~ lonAvgBinom180, data = envinddatana,
        main = "Migration ", ylab = MyYlab)

par(op)


op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals"
plot(x = F1, y = R1, xlab = "Fitted values", ylab = MyYlab)
boxplot(R1 ~ sic5day, data = envinddatana,
        main = "SIC ", ylab = MyYlab)

boxplot(R1 ~ sumO5, data = envinddatana,
        main = "Other dive ", ylab = MyYlab)

boxplot(R1 ~ sumF5, data = envinddatana,
        main = "Foraging dive ", ylab = MyYlab)

boxplot(R1 ~ sex.f, data = envinddatana,
        main = "Sex ", ylab = MyYlab)

boxplot(R1 ~ yearweek, data = envinddatana,
        main = "Week ", ylab = MyYlab)

par(op)

