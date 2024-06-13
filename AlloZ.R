#---------------------------------------------------------------------------------------------------------
#Import des données
alloz <- read.csv("CHERYL/dataAlloZ/dataAlloZ.csv", 
                 header = TRUE, 
                 sep = ";", 
                 encoding = "latin1", na.strings = c("", " ", "NA", "NI"))

#Formattage en date
alloz$death.dt <- as.Date(as.character(alloz$deathdt.T),format="%Y-%m-%d")
alloz$rando.dt <- as.Date(as.character(alloz$datrand),format="%Y-%m-%d")
alloz$max.dt <- as.Date(as.character(alloz$daderdt.T),format="%Y-%m-%d")
#qaly$rech.dt <- as.Date(as.character(qaly$relapsdt),format="%Y-%m-%d")



# Délais calculés
alloz$suivi  <- as.numeric(alloz$max.dt-alloz$rando.dt)/ (365.25/12)
alloz$survie <- as.numeric(alloz$death.dt-alloz$rando.dt)/ (365.25/12)
alloz$survie[is.na(alloz$survie)] <- alloz$suivi[is.na(alloz$survie)]
alloz$dc<- ifelse(is.na(alloz$death.dt),0,1)

summary(alloz)
alloz$Arm<-factor(alloz$Arm)


library(survival)
library(survminer)
os<-survfit(Surv(survie,dc)~Arm, data=alloz)
ggsurvplot(os)

# coxph_model<-coxph(Surv(survie,dc)~Arm, data=alloz)
# t<-cox.zph(coxph_model)
# ggcoxzph(t)


library(cuRe)
cm<-fit.cure.model(Surv(survie,dc)~Arm,formula.surv =list(~Arm,~1) , data=alloz, dist = "weibull",link="logit",method="L-BFGS-B")
summary(cm)

plot(cm,newdata=data.frame(Arm = c(unique(alloz$Arm))),ci=F)
