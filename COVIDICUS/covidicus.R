library(dplyr)
library(rpact)
library(survival)
library(nph)
library(foreach)
library(doRNG)
library(flexsurv)

ensai_covid <- read.csv("~/M2 SMDS/Stage M2/COVIDICUS/ensai_covid.csv", sep=";")
str(ensai_covid)
ensai_covid$tt <- ifelse(ensai_covid$suivi<60,ensai_covid$suivi,60)

# Dates de recrutement manquantes ####
new_base <- data.frame()
for(j in 1:117){
i <- as.character(j)
length(i)
if(nchar(i)==1){
  valeur <- paste0("000",i)
}else if(nchar(i)==2){
  valeur <- paste0("00",i)
}else if(nchar(i)==3){
  valeur <- paste0("0",i)
}
new_vector <- ensai_covid[word(ensai_covid$SUBJECT_REF, start = 2, sep="-")==valeur,]
new_base <- rbind(new_base, new_vector)
}
View(new_base)

# quel sample size pour essai fixe ####
end_accrual_after_nb_month <- 10*30
ta <- end_accrual_after_nb_month
end_study_after_nb_month <- 12*30
tps_restant_apres_dernier_recrutement <- end_study_after_nb_month - end_accrual_after_nb_month
beta <- 0.2
alpha <- 0.025

coxph(Surv(tt, dc60)~DXM, data = ensai_covid)$coef
d <- (4*(qnorm(1-alpha)+qnorm(1-beta))^2) / as.numeric(coxph(Surv(tt, dc60)~DXM, data = ensai_covid)$coef)^2 
lambda1_reestime <- sum(filter(ensai_covid,DXM=="A")$dc60)/sum(filter(ensai_covid,DXM=="A")$tt)
lambda0_reestime <- sum(filter(ensai_covid,DXM=="B")$dc60)/sum(filter(ensai_covid,DXM=="B")$tt)
P0 <- 1 - (1/(end_accrual_after_nb_month*lambda0_reestime))*(exp(-lambda0_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda0_reestime*end_accrual_after_nb_month)))
P1 <- 1 - (1/(end_accrual_after_nb_month*lambda1_reestime))*(exp(-lambda1_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda1_reestime*end_accrual_after_nb_month)))
(N_revised <- 2*d / (P0+P1)) #24643.03

filter(ensai_covid,Type_ventil=="MV")
dim(filter(ensai_covid,Type_ventil=="MV"))

d <- (4*(qnorm(1-alpha)+qnorm(1-beta))^2) / as.numeric(coxph(Surv(tt, dc60)~DXM, data = filter(ensai_covid,Type_ventil=="MV"))$coef)^2 
lambda1_reestime <- sum(filter(filter(ensai_covid,Type_ventil=="MV"),DXM=="A")$dc60)/sum(filter(filter(ensai_covid,Type_ventil=="MV"),DXM=="A")$tt)
lambda0_reestime <- sum(filter(filter(ensai_covid,Type_ventil=="MV"),DXM=="B")$dc60)/sum(filter(filter(ensai_covid,Type_ventil=="MV"),DXM=="B")$tt)
P0 <- 1 - (1/(end_accrual_after_nb_month*lambda0_reestime))*(exp(-lambda0_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda0_reestime*end_accrual_after_nb_month)))
P1 <- 1 - (1/(end_accrual_after_nb_month*lambda1_reestime))*(exp(-lambda1_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda1_reestime*end_accrual_after_nb_month)))
(N_revised <- 2*d / (P0+P1))#2646.262



d <- table(ensai_covid$dc60)[2]

# permutation ####
set.seed(27011996)
variable <- ensai_covid$SUBJECT_REF # variable we will resample from 
permut_nbr <- 10000  #number of permutation
obs_nbr <- dim(ensai_covid)[1] # number of observations

#initialize a matrix to store the permutation data 
PermSamples <- matrix(0,nrow=obs_nbr, ncol = permut_nbr)

for(i in 1:permut_nbr){
  PermSamples[,i] <- sample(variable, size = obs_nbr, replace=F)
}
PermSamples[,1]

PermSamples[ensai_covid$dc60==1,1]

ne

order(ensai_covid$SUBJECT_REF==PermSamples[,1])
ensai_covid_v2<-ensai_covid[match(PermSamples[,1], ensai_covid$SUBJECT_REF),]

ensai_covid_v2<-ensai_covid$SUBJECT_REF[PermSamples[,1]]


# Sample size reassessment a 6 mois ####
end_accrual_after_nb_month <- 10*30
ta <- end_accrual_after_nb_month
end_study_after_nb_month <- 12*30
tps_restant_apres_dernier_recrutement <- end_study_after_nb_month - end_accrual_after_nb_month
beta <- 0.2
alpha <- 0.025

data <- new_base

table(new_base$classBMI)
mean(new_base$age)

data <- filter(new_base,classBMI=="Moderate Obesity" | classBMI=="Severe Obesity")
data <- filter(new_base,classBMI=="Normal")
data <- filter(new_base, Type_ventil=="MV")
data <- filter(new_base, Type_ventil=="Optiflow")
data <- filter(new_base, Type_ventil=="O2")
data <- filter(new_base, age>=65)
data <- filter(new_base, age<65)

lambda1 <- sum(filter(data,DXM=="A")$dc60)/sum(filter(data,DXM=="A")$tt)
lambda0 <- sum(filter(data,DXM=="B")$dc60)/sum(filter(data,DXM=="B")$tt)
temps_informatifs <- 1/2
TAU <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                    end_study_after_nb_month = 12*30, end_accrual_after_nb_month = 10*30)$TAU)
IA <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                   end_study_after_nb_month = 12*30, end_accrual_after_nb_month = 10*30)$IA)
vect_early_stop <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                                end_study_after_nb_month = 12*30,
                                                end_accrual_after_nb_month = 10*30)$vect_early_stop)
vect_final_stop <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                                end_study_after_nb_month = 12*30,
                                                end_accrual_after_nb_month = 10*30)$vect_final_stop)
vect_c_alpha2 <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                              end_study_after_nb_month = 12*30,
                                              end_accrual_after_nb_month = 10*30)$vect_c_alpha2)
eventsPerStage <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                               end_study_after_nb_month = 12*30,
                                               end_accrual_after_nb_month = 10*30)$eventsPerStage)
interim_SS <- as.numeric(Get_Analysis_Time(TAU = temps_informatifs, alpha, beta, lambda0, lambda1,
                                           end_study_after_nb_month = 12*30,
                                           end_accrual_after_nb_month = 10*30)$interim_SS)

tau <- 1/2 # étude à 6 mois --> les patients recrutés jusqu'à 4 mois auront été suivis 2 mois
tps_interim <- 30*6
early_stop <- vect_early_stop
final_stop <- vect_final_stop
c_alpha2 <- vect_c_alpha2
ePS <- ceiling(sum(data$dc60)*tau)
n1_planned <- dim(data)[1]*tau
(d <- (4*(qnorm(1-alpha)+qnorm(1-beta))^2) / as.numeric(coxph(Surv(tt, dc60)~DXM, data = data)$coef)^2)
#(d <- (4*(qnorm(1-alpha)+qnorm(1-beta))^2) / log(0.9)^2)
d <- table(data$dc60)[2]
w1 <- sqrt((tau * d) / d) # proportion of events expected to be observed at IA
w2 <- sqrt(1 - w1^2)

# Expected number of events Jorgens
e_11 <- tau * d
e_12 <-  round(n1_planned / N * d) -  e_11#round(dim(filter(data, recrutement <= tps_interim))[1] / N * d) -  e_11
e_22 <- pmax(0, d - e_11 - e_12) # pas en fonction du nouveau d car les pondérations doivent etre préspécifiées pour préserver le risque d'erreur de type 1
# Jorgens preplanned Weights
w11 <- sqrt(e_11 / (e_11 + e_12 + e_22))
if (e_12 <= 0) {
  w12 <- 0
}else{
  w12 <- sqrt((e_12 / e_11) * w11^2)
}
if (e_22 != 0) {
  w22 <- sqrt(1 - (w11^2 + w12^2))
}
nb_simu = 1
prop_rejet_WASSMER <- c(rep(NA, nb_simu))
prop_efficacy <- c(rep(NA, nb_simu))
prop_rejet_DESSEAUX <- c(rep(NA, nb_simu))
prop_rejet_JORGENS <- c(rep(NA, nb_simu))
prop_rejet_NAIVE <- c(rep(NA, nb_simu))
SS2  <- c(rep(NA, nb_simu))
adaptation <- c(rep(NA, nb_simu))
logHR <- c(rep(NA, nb_simu))
classique <- c(rep(NA, nb_simu))


sous_base <- data[1:n1_planned,] 
d1 <- sum(sous_base$dc60)


LR_1 <- as.numeric(logrank.test(time = sous_base$tt,
                                event = sous_base$dc60,
                                group = sous_base$DXM,
                                alternative = "greater")$test[3])


(p1 <- 1 - pnorm(LR_1))
n1 <- dim(sous_base)[1]


if (p1 < early_stop) {
  early_stopping_D <- "efficacy"
  decision_D <- "reject H0"
  COMBI_D <- NA
  
  decision <- "reject H0"
  
  
  p_final <- NA
  decision_J <- "reject H0"
  
  decision_naive <- 1
  
  SS2[i] <- dim(sous_base)[1]
  adaptation[i] <- 0
  classique[i] <- 0
}  else{ # proceeds to stage 2
  early_stopping_D <- "NO"
  
  b_etoile <- qnorm(1-final_stop)
  CP_interim = 1 - pnorm( (b_etoile*sqrt(d/4) - LR_1*sqrt(d1/4) - ((d-d1)/4)*LR_1/sqrt(d1/4)) / sqrt((d-d1)/4) )
  
  d2_etoile <- (d1/LR_1^2) * (( (b_etoile*sqrt(d)-LR_1*sqrt(d1)) / sqrt(d-d1)) + qnorm(1-beta) )^2
  d_fin2 <- d2_etoile 
  lambda1_reestime <- sum(filter(sous_base,DXM=="B")$dc60)/sum(filter(sous_base,DXM=="B")$tt)
  lambda0_reestime <- sum(filter(sous_base,DXM=="A")$dc60)/sum(filter(sous_base,DXM=="A")$tt)
  P0 <- 1 - (1/(end_accrual_after_nb_month*lambda0_reestime))*(exp(-lambda0_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda0_reestime*end_accrual_after_nb_month)))
  P1 <- 1 - (1/(end_accrual_after_nb_month*lambda1_reestime))*(exp(-lambda1_reestime*tps_restant_apres_dernier_recrutement)*(1-exp(-lambda1_reestime*end_accrual_after_nb_month)))
  N_revised <- 2*d_fin2 / (P0+P1) 
  
  (N_revised <- ceiling(N_revised))
  n2 <- max(N_revised - n1,  0)#max(N_revised - n1,  N - n1)
  adaptation[i] <- 1
  

}
dim(data)
d_fin2
N_revised
n1
n2

