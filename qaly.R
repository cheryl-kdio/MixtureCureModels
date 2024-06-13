#---------------------------------------------------------------------------------------------------------
#Import des données
qaly <- read.csv("GRAAPH/graaphR_ensai.csv", 
                 header = TRUE, 
                 sep = ";", 
                 encoding = "latin1", na.strings = c("", " ", "NA", "NI"))

#Formattage en date
qaly$death.dt <- as.Date(as.character(qaly$deathdt),format="%Y-%m-%d")
qaly$rando.dt <- as.Date(as.character(qaly$randodt),format="%Y-%m-%d")
qaly$max.dt <- as.Date(as.character(qaly$Datemax),format="%Y-%m-%d")
qaly$rech.dt <- as.Date(as.character(qaly$relapsdt),format="%Y-%m-%d")



# Délais calculés
qaly$suivi  <- as.numeric(qaly$max.dt-qaly$rando.dt)/ (365.25/12)
#qaly$survie <- as.numeric(qaly$death.dt-qaly$rando.dt)/ (365.25/12)
qaly$delrech<- as.numeric(qaly$rech.dt-qaly$rando.dt)/ (365.25/12)
qaly$delpfs<- qaly$suivi
qaly$delpfs[!is.na(qaly$delrech)]<-qaly$delrech[!is.na(qaly$delrech)]
qaly$R1<-factor(qaly$R1,levels=c("Intensive arm (A)","Light arm (B)"))


#Délais médians
summary(qaly$suivi)
summary(qaly$delpfs)

library(survRM2)
rmst2(qaly$delpfs, qaly$pfs, as.numeric(qaly$R1)-1, covariates = NULL, alpha = 0.05)


#---------------------------------------------------------------------------------------------------------
#-------- Estimation des courbes de survie KM
library(survival)
library(survminer)
efs <- survfit(Surv(delpfs,pfs)~R1, data=qaly)

efs_p<-ggsurvplot(efs, data=qaly,
                  risk.table = TRUE,
                  ggtheme = theme_bw(),
                  conf.int = F,
                  xlab = "Durée de suivi (mois)",
                  ylab = "Probabilité de survie")

efs_p$plot+
  coord_cartesian(xlim=c(0,71.49),ylim=c(0,1),expand = F)+
  scale_color_manual(values=c("#009988","#CC6677"))


# plot cum hazard risk
# ggsurvplot(efs, data=qaly,
#            risk.table = TRUE,
#            ggtheme = theme_bw(),
#            conf.int = F,
#            fun = "cumhaz",
#            xlab = "Durée de suivi (mois)",
#            ylab = "Risque cumulé de rechute")

# print odds ratio
coxph_model<-coxph(Surv(delpfs,pfs)~R1, data=qaly)


#Test de log rank
survdiff(Surv(delpfs,pfs)~R1, data=qaly)

#-------- Rapport de risques instantanés
coxph_model<-coxph(Surv(delpfs,pfs)~R1, data=qaly)
summary(coxph_model)

# Test de risques proportionnels
t<-cox.zph(coxph_model)
ggcoxzph(t)


ggsurvplot(survfit(x, newdata = data.frame(R1 = c("Intensive arm (A)","Light arm (B)"))),data=qaly,conf.int = F)

plot(cure_model,newdata=data.frame(R1 = c("Intensive arm (A)","Light arm (B)")),
     ci=F,type="surv",
     lwd=2,lty=2,
     ylab="Probabilité de survie",
     xlab="Durée de suivi (mois)",
     xaxs = "i",
     yaxs = "i",
     col="red")
lines(survfit(x, newdata = data.frame(R1 = c("Intensive arm (A)","Light arm (B)"))),lwd=2,col="blue")
lines(efs,col=c("#009988","#CC6677"),conf.int=F,lwd=3)
legend("topright",legend=c("Intensive arm (A)","Light arm (B)"),col=c("#009988","#CC6677"),lwd=3)


### pvaleur = 0.0316, HR=1.86 
### => Patients du bras B ont un risque de rechute (statistiquement) de 1.86 fois plus élevé que les patients du bras A

#---------------------------------------------------------------------------------------------------------
#-------- Modèle de cure de mélange with cuRe

library(cuRe)
cure_model<-fit.cure.model(Surv(delpfs,pfs)~R1,formula.surv =list(~R1,~1) , data=qaly, dist = "weibull",link="logit",method="L-BFGS-B")

AIC(cure_model)

summary(cure_model)
plot(cure_model,newdata=qaly,ci=F,type="surv")

#-------- Estimations
### Probability of being cured for group A
exp(cure_model$coefs$`1`)/(1+exp(cure_model$coefs$`1`))

### Probability of being cured for group B gamma1 + gamma2
pi_u<-c(1, 1, 0, 0, 0) %*% unlist(cure_model$coefs)
exp(pi_u)

#OR and HR
exp(unlist(cure_model$coefs)[c(2,4)])

### IC
#variance de gamma1 + gamma2
var_u<-c(1, 1, 0, 0, 0) %*% cure_model$covariance %*% c(1, 1, 0, 0, 0)

exp(pi_u- 1.96*var_u)/(1+exp(pi_u- 1.96*var_u))
exp(pi_u+ 1.96*var_u)/(1+exp(pi_u+ 1.96*var_u))

#-------- Courbes de survie estimées superposées aux courbes de survie KM
plot(cure_model,newdata=qaly,ci=F,type="surv",lwd=2,col="red",lty=2,
     ylab="Probabilité de survie",
     xlab="Durée de suivi (mois)",
     xaxs = "i",
     yaxs = "i")
lines(efs,col=c("#009988","#CC6677"),conf.int=F,lwd=3)
legend("topright",legend=c("Intensive arm (A)","Light arm (B)"),col=c("#009988","#CC6677"),lwd=3)

#---------------------------------------------------------------------------------------------------------
#-------- Ratio test in Intensive arm (A)
whichTau_A<-max(qaly$delpfs[qaly$R1=="Intensive arm (A)"])
pi_A<-exp(unlist(cure_model$coefs)[1])/(1+exp(unlist(cure_model$coefs)[1]))
scale_A<-exp(unlist(cure_model$coefs)[3]) 
shape<-exp(unlist(cure_model$coefs)[5])
ratio_test_A<- pweibull(whichTau_A,shape=shape,scale=1/scale_A,lower.tail=FALSE)/
                 (pi_A+(1-pi_A)*pweibull(whichTau_A,shape=shape,scale=1/scale_A,lower.tail=FALSE))

print(pi_A)
print(ratio_test_A)

#-------- Ratio test in Light arm (B)
whichTau_B<-max(qaly$delpfs[qaly$R1=="Light arm (B)"])
pi_B<-exp(unlist(cure_model$coefs)[2]+unlist(cure_model$coefs)[1])/(1+exp(unlist(cure_model$coefs)[2]+unlist(cure_model$coefs)[1]))
scale_B<-exp(unlist(cure_model$coefs)[4]+unlist(cure_model$coefs)[3])
ratio_test_B<- pweibull(whichTau_B,shape=shape,scale=1/scale_B,lower.tail=FALSE)/
                 (pi_B+(1-pi_B)*pweibull(whichTau_B,shape=shape,scale=1/scale_B,lower.tail=FALSE))

print(pi_B)
print(ratio_test_B)


#----------- Test Maller&Zhou => sufficient follow-up
qalyA<-qaly[which(qaly$R1=="Intensive arm (A)"),]
maxE <- max(qalyA$delpfs[qalyA$pfs==1]) # last event time
maxAll <- max(qalyA$delpfs) # last observation time
qn<-NA
if (maxAll > maxE){
  numPlat <- sum(qalyA$pfs==1 & qalyA$delpfs > (2*maxE-maxAll) & qalyA$delpfs <= maxE) # number of events that are between (2*maxE-maxAll) (exclusive) and maxE (inclusive)
  qn<-numPlat/nrow(qalyA) # qn
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



#---------------------------------------------------------------------------------------------------------
#-------- Modèle de cure de mélange with flexsurvcure

library(flexsurvcure)
flexcure_model<-flexsurvcure(Surv(delpfs, pfs)~R1, data=qaly, dist="weibull",link="logistic", mixture=T,anc = list(shape = ~ R1,scale= ~ 1))
print(flexcure_model)

############################################################################################################
# # Courbes de survie des patients non guéris
# predict(cure_model,newdata=qaly,type="survuncured")
# 
# qalyA<-qaly[which(qaly$R1=="Intensive arm (A)"),]
# qalyB<-qaly[which(qaly$R1=="Light arm (B)"),]
# 
# 
# xlim <- c(1e-05, max(cure_model$time))
# time <- seq(xlim[1], xlim[2], length.out = 100)
# 
# predA <- predict(cure_model, newdata = qalyA, type = "survuncured",time=time,var.type="se")
# predB <- predict(cure_model, newdata = qalyB, type = "survuncured",time=time)
# 
# # Compute area under the survival curve
# area_surv <- function(x, y) {
#   step <- x[2] - x[1]
#   len <- length(y)
#   val <- sum(y[2:(len - 1)])
#   val <- val + 0.5 * (y[1] + y[len])
#   val <- step * val
#   
#   return(val)
# }
# 
# 
# # Compare RMST truncated at time 71.49
# ttime <- time[time <= 71.49]
# 
# # RMST given with KM
# library(survRM2)
# rmst2(qaly$delpfs, qaly$pfs, as.numeric(as.factor(qaly$R1))-1, covariates = NULL, alpha = 0.05)
# 
# 
# # RMST given with cure model
# area_surv(ttime,  predA[[1]]$Estimate[1:length(ttime)])
# area_surv(ttime, predB[[1]]$Estimate[1:length(ttime)])
# 
# #bootstrap to compute rmst IC
# compute_rmst<-function(data, indices,cure_model, ttime) {
#     boot_data <- data[indices, ]
#     pred_boot <- predict(cure_model, newdata = boot_data, type = "survuncured")
#     
#     rmst <- area_surv(ttime, pred_boot[[1]]$Estimate[1:length(ttime)])
#     return(rmst)
#   }
#   
# library(boot)
# 
# boot_obj<-boot(qalyA, compute_rmst, R = 50, cure_model=cure_model, ttime=ttime)
# boot.ci(boot_obj, type = "norm")
# 
# boot
############################################################################################################
# Modèle de cure de mélange semi paramétrique



library(smcure)
qaly$R1.f<-(as.numeric(qaly$R1)-1)
cure_model_sp<-smcure(Surv(delpfs,pfs)~R1.f,cureform = ~R1.f, data=qaly[,c("delpfs","pfs","R1.f")], model="ph",Var = T,link="logit")

print(cure_model_sp$logistfit)
summary(cure_model_sp)

data(e1684)
pd <- smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,
             cureform=~TRT+SEX+AGE,data=e1684,model="ph",
             Var = FALSE)
printsmcure(pd,Var = FALSE)