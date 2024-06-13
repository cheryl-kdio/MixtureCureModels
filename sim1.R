################################################################################
# Simulation d'un modèle de guerison

#' This routines simulates survival data with a cure fraction. The data is
#' simulated according to a mixture cure model with a logistic link for the
#' incidence part of the model and the latency part assumes a Weibull model with
#' no covariates.


set.seed(2022)
n = 500; wshape = 1.4; wscale = 0.05


#--- Incidence part (logistic model)
beta<- 0.7
censrate <- 0.1
tau0 <- 40

p <- exp(beta)/(1 + exp(beta)) 

#--- Generation of cure status
B <- stats::rbinom(n, size = 1, prob = p) # Cure status B = 1 --> cured ### prob=p or 0

#--- Generation of survival times for uncured (B=0) subject from Weibull

## Draw survival times from Weibull distribution
weibshape <- wshape # Weibull shape parameter > 0
weibscale <- wscale # Weibull scale parameter > 0
S0 <- function(t) exp(-weibscale * t ^ (weibshape)) # True baseline survival

# Generation of survival times for uncured population
U <- stats::runif(n) # Uniform draw
Tlat <- as.numeric((-log(U) / weibscale) ^
                     (1 / weibshape))
Tlat[which(B == 1)] <- 20000 # Large survival time for cured subjects  ### Delete this line to eliminate cure fraction
tobs <- Tlat

#--- Censoring follows exponential distribution with rate lambda

C <- stats::rexp(n,rate=0.15)
TgreatC <- which(Tlat > C)
tobs[TgreatC] <- C[TgreatC]
delta <- as.numeric(Tlat <= C)

#--- Summary statistics on cure and censoring rates
dataKM <- as.data.frame(cbind(tobs, delta))
fitKM <- survival::survfit(survival::Surv(tobs, delta) ~ 1, data = dataKM)
library(survminer)
ggsurvplot(fitKM)

plateau <-
  fitKM$time[utils::tail(which((diff(fitKM$surv) < 0) == TRUE), 1) + 1]
nobs_plateau <- sum(tobs > plateau)

cure_rate<-round(sum(B == 1) / n, 3) 

plot(tobs,S0(tobs))

# Extract variables
simdata <- data.frame(tobs, delta)
colnames(simdata) <- c("tobs","event")

#--- Fit a weibull/logit mixture cure model
library(cuRe)
cm<-fit.cure.model(Surv(tobs, delta) ~ 1, data = simdata, dist = "weibull", link = "logit", method = "L-BFGS-B")
summary(cm)

#--- Fit a Weibull model
library(flexsurv)
wei_reg<-flexsurvreg(Surv(tobs, delta) ~ 1, dist = "weibull", data = simdata)
summary(wei_reg)


#--- Plot the estimated survival curves vs the true one
plot(cm,type="surv",ci=F,col="blue",lwd=2)
lines(wei_reg,col="green",ci=F)
lines(fitKM, conf.int = F,col="red",lwd=2)
legend("topright",legend=c("Cure model","Weibull model","Kaplan-Meier"),col=c("blue","green","red"),lwd=2)

