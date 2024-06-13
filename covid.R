#---------------------------------------------------------------------------------------------------------
#Import des données
covid <- read.csv("COVIDICUS/covidicus.csv", 
                header = TRUE, 
                sep = ";", 
                encoding = "latin1", na.strings = c("", " ", "NA", "NI"))


covid$rando.dt <- as.Date(as.character(covid$randodt),format="%Y-%m-%d")
covid$max.dt <- as.Date(as.character(covid$Dader),format="%Y-%m-%d")

covid$survie  <- as.numeric(covid$max.dt-covid$rando.dt)
covid$survie[which(covid$survie>60)] <- 60


# covid60<-data.frame(survie=survie60,dc=dc60,DXM=DXM60)

#covid$DXM<-factor(covid$DXM)
# summary(covid60)
library(survival)
library(survminer)
os<-survfit(Surv(survie,dc)~DXM, data=covid)
os_p<-ggsurvplot(os)
ggsurvplot(os)


library(survRM2)
rmst2(covid$survie,covid$dc,as.numeric(as.factor(covid$DXM))-1,tau = 60)

#----------- Hazard ratio
cox<-coxph(Surv(survie,dc)~DXM,data=covid)
ggsurvplot(cox)

#---- Test hypothhèse de hazard proportionnel
t<-cox.zph(coxph(Surv(survie,dc)~DXM,data=covid))
ggcoxzph(t)

#----------- Modèle de cure
library(cuRe)
cm<-fit.cure.model(Surv(survie,dc)~DXM,data=covid,formula.surv =list(~DXM,~1) , dist = "weibull",link="logit",method="L-BFGS-B")

summary(cm)
plot(cm,newdata=covid)

#----------- Test de ratio de survie
whichTau_A<-max(covid$survie[covid$DXM=="A"])
pi_A<-exp(unlist(cm$coefs)[1])/(1+exp(unlist(cm$coefs)[1]))
scale_A<-exp(unlist(cm$coefs)[3]) 
shape<-exp(unlist(cm$coefs)[5])
ratio_test_A<- pweibull(whichTau_A,shape=shape,scale=1/scale_A,lower.tail=FALSE)/
  (pi_A+(1-pi_A)*pweibull(whichTau_A,shape=shape,scale=1/scale_A,lower.tail=FALSE))

print(pi_A)
print(ratio_test_A)

#----------- Test Maller&Zhou => sufficient follow-up
DMXA<-covid[covid$DXM=="A",]
maxE <- max(DMXA$survie[DMXA$dc==1]) # last event time
maxAll <- max(DMXA$survie) # last observation time
qn<-NA
if (maxAll > maxE){
  numPlat <- sum(DMXA$dc==1 & DMXA$survie > (2*maxE-maxAll) & DMXA$survie <= maxE) # number of events that are between (2*maxE-maxAll) (exclusive) and maxE (inclusive)
    qn<-numPlat/nrow(DMXA) # qn
  }
print(qn)

#----------- Test Maller&Zhou => presence of immunes
kmfit <- survfit(Surv(DMXA$survie,DMXA$dc)~1)

lastObs <- max(DMXA$survie) 
isCens <- (DMXA$dc[which.max(DMXA$survie)]==0) # test if the last observation is a censoring, sanity check that phat < 1
pHat <- 1-summary(kmfit,times=lastObs)$surv
pCens <- 1-mean(DMXA$dc) # proportion of censored

immune_test<-c(pHat,pCens,lastObs,isCens)

(1/immune_test[2])-1


