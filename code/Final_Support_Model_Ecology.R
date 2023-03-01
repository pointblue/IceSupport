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
library(nlme)

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


##############################
#########################Distance by ice support (aka Projection)
##############################################
## STEP 1
#BASIC LINEAR MODELS with mixed or non mixed assessment
##Enviro only models without model weights
##Using REML here as per Zuur et al for comparing mixed vs non mixed models


mixmod<-lme(support5day ~  age + sex.f + Breeding.f  + poly(yearweek,2)*lonAvgBinom180  + year.f + divesum5 , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.omit, method = "REML")


nonmixmod<-gls(support5day ~  age + sex.f + Breeding.f  + poly(yearweek,2)*lonAvgBinom180    + year.f + divesum5 , weights = ~ err5day, data = envinddatana, na.action=na.omit , method = "REML")
summary(nonmixmod)

anova(mixmod, nonmixmod)


######################################
###STEP 2
## Find optimal fixed structure using an ML testing procedure
###In this case make sure to us ML for comparison of nested models



mixmodsupport<-lme(support5day ~ age + sex.f + Breeding.f   + poly(yearweek,2)*lonAvgBinom180  + year.f + d , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")


summary(mixmodsupport)
r2(mixmodsupport)
effectsnonmixmod<-effects::allEffects(mixmodsupport)
plot(effectsnonmixmod)
intervals(mixmodsupport)

supportdredge <- MuMIn::dredge(mixmodsupport, rank="AIC")
subset(supportdredge, delta <2) # show the top models within 2 AIC


#1 top
mixmodsupport1<-lme(support5day ~ poly(yearweek,2) + lonAvgBinom180  + year.f , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
intervals(mixmodsupport1,which = 'fixed')
summary(mixmodsupport1)
r2(mixmodsupport1)

#2 add breeding
mixmodsupport2<-lme(support5day ~ poly(yearweek,2) + lonAvgBinom180 + year.f + Breeding.f  , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
intervals(mixmodsupport2,which = 'fixed')
summary(mixmodsupport2)
r2(mixmodsupport2)

#3 add divesum
mixmodsupport3<-lme(support5day ~  poly(yearweek,2) + lonAvgBinom180  + year.f + divesum5 , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
intervals(mixmodsupport3,which = 'fixed')
summary(mixmodsupport3)
r2(mixmodsupport3)


#4 add age
mixmodsupport4<-lme(support5day ~  poly(yearweek,2) + lonAvgBinom180  + year.f + age , random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "ML")
intervals(mixmodsupport4,which = 'fixed')
summary(mixmodsupport4)
r2(mixmodsupport4)

anova(mixmodsupport1, mixmodsupport2, mixmodsupport3, mixmodsupport4)




####################################################
##STEP 3
#Run top model with REML


mixmodsupporttop<-lme(support5day ~ poly(yearweek,2) + lonAvgBinom180 +  year.f, random=~1|bird_fn, weights = ~ I(1/err5day), data = envinddatana, na.action=na.fail, method = "REML")


summary(mixmodsupporttop)
r2(mixmodsupporttop) #Marginal R2 provides the variance explained only by fixed effects and conditional R2 provides the variance explained by the entire model, i.e., both fixed effects and random effects.
effectsmixmodtop<-effects::allEffects(mixmodsupporttop)
plot(effectsmixmodtop)
intervals(mixmodsupporttop,which = 'fixed')


#predict values
ggpredict(mixmodsupporttop, terms= c('year.f'), back.transform = FALSE)
ggpredict(mixmodsupporttop, terms= c('lonAvgBinom180'), back.transform = FALSE)


png(filename="V:/Project/Terrestrial/adpe/nasa_winter_ecology/ice_support/finalAnalysis/images/topModel_Support_yearweeka.png", width=600, height=600)
ggemmeans(mixmodsupporttop, terms= "yearweek", back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 0.1) + labs( x = "Week",  y = "Support", title = "(a)") + scale_y_continuous(limits = c(-0.05,0.0)) + My_Theme 
dev.off()

png(filename="topModel_Support_Yearb.png", width=600, height=600)
ggemmeans(mixmodsupporttop, terms= "year.f", back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 0.1) + labs( x = "Year",  y = "Support", title = "(b)") +  scale_y_continuous(limits = c(-0.05,0.0)) + My_Theme 
dev.off()



png(filename="topModel_Support_Migrationb.png", width=600, height=600)
ggemmeans(mixmodsupporttop, terms= "lonAvgBinom180", back.transform = FALSE)%>%
  plot(add.data = FALSE, dot.size = 5) + labs( x = "Migration Route",  y = "Support", title = "(b)") +  scale_y_continuous(limits = c(-0.05,0.0)) + scale_colour_manual(values = cbbPalette)  + My_Theme 
dev.off()



#####################################################
##STEP 4
#Model validation with REML model
#use lme4 lmer to be able to use check_model function
#Plotting of residuals

mixmodsupporttoplme4<- lme4::lmer(support5day ~ poly(yearweek,2) + lonAvgBinom180 + sumO5 + sumE5 + year.f + (1|bird_fn), weights = err5day, data = envinddatana, REML = TRUE, na.action = na.fail)
check_model(mixmodsupporttoplme4)


summary(mixmodsupporttoplme4)


plot(residuals(mixmodsupporttoplme4, type= "normalized"))
R1<-residuals(mixmodsupporttoplme4, type= "normalized")
F1<-fitted(mixmodsupporttoplme4)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))
MyYlab <- "Residuals"
plot(x = F1, y = R1, xlab = "Fitted values", ylab = MyYlab)
boxplot(R1 ~ sumF5, data = envinddatana,
        main = "sumF5 ", ylab = MyYlab)

boxplot(R1 ~ sumO5, data = envinddatana,
        main = "sumO5 ", ylab = MyYlab)

boxplot(R1 ~ year.f, data = envinddatana,
        main = "Year ", ylab = MyYlab)

boxplot(R1 ~ yearweek, data = envinddatana,
        main = "yearweek ", ylab = MyYlab)


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


