####################################
# Data Analysis for S0106 and S1117
#
# Subodh Selukar
# 2022-03-24
####################################

library(survival)
library(flexsurvcure)

labSize <- 1.5

roundP <- function(num,dig=3){
  numRound <- round(num,dig)
  if (numRound==0) return(paste0("<",10**(-dig))) else return(numRound)
}


## Wrapper function for computing MLE
mleFun <- function(dat,dist="exp"){ # expects a df with columns Y,D
  if (dist=="exp"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="exp")
    return(c(pi=tmp$res[1,1],rate=tmp$res[2,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             rateLB=tmp$res[2,2],rateUB=tmp$res[2,3]))
  }
  if (dist=="wei"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="weibull")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],scale=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             scaleLB=tmp$res[3,2],scaleUB=tmp$res[3,3]))
  }
  if (dist=="llogis"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="llogis")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],scale=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             scaleLB=tmp$res[3,2],scaleUB=tmp$res[3,3]))
  }
  if (dist=="gam"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="gamma")
    return(c(pi=tmp$res[1,1],shape=tmp$res[2,1],rate=tmp$res[3,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             shapeLB=tmp$res[2,2],shapeUB=tmp$res[2,3],
             rateLB=tmp$res[3,2],rateUB=tmp$res[3,3]))
  }
  if (dist=="g.gam"){
    tmp <- flexsurvcure(Surv(Y,D)~1,data=dat,dist="gengamma")
    return(c(pi=tmp$res[1,1],mu=tmp$res[2,1],sigma=tmp$res[3,1],shape=tmp$res[4,1],
             piLB=tmp$res[1,2],piUB=tmp$res[1,3],
             muLB=tmp$res[2,2],muUB=tmp$res[2,3],
             sigmaLB=tmp$res[3,2],sigmaUB=tmp$res[3,3],
             shapeLB=tmp$res[4,2],shapeUB=tmp$res[4,3]))
  }
}

ratioTest <- function(dat,whichTau,dist="exp"){# expects a df with columns Y,D and administrative censoring time
  est <- mleFun(dat,dist)
  
  if (dist=="exp"){
    return(
      c(est,
        pexp(whichTau,rate=est[2],lower.tail=FALSE)/
          (est[1]+(1-est[1])*pexp(whichTau,rate=est[2],lower.tail=FALSE)))
    )
  }
  if (dist=="wei"){
    return(
      c(est,
        pweibull(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)/
          (est[1]+(1-est[1])*pweibull(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)))
    )
  }
  if (dist=="llogis"){
    return(
      c(est,
        pllogis(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)/
          (est[1]+(1-est[1])*pllogis(whichTau,shape=est[2],scale=est[3],lower.tail=FALSE)))
    )
  }
  if (dist=="gam"){
    return(
      c(est,
        pgamma(whichTau,shape=est[2],rate=est[3],lower.tail=FALSE)/
          (est[1]+(1-est[1])*pgamma(whichTau,shape=est[2],rate=est[3],lower.tail=FALSE)))
    )
  }
  if (dist=="g.gam"){
    return(
      c(est,
        pgengamma(whichTau,mu=est[2],sigma=est[3],Q=est[4],lower.tail=FALSE)/
          (est[1]+(1-est[1])*pgengamma(whichTau,mu=est[2],sigma=est[3],Q=est[4],lower.tail=FALSE)))
    )
  }
}

ratioOut <- function(dat,whichTau,dist="exp") {# expects a df with columns Y,D and administrative censoring time
  tmp <- ratioTest(dat,whichTau,dist)
  tmp[length(tmp)]
}

## M-Z test statistic Maller Zhou (1994)
mzTest <- function(dat){ # expects a df with columns Y,D
  maxE <- max(dat$Y[dat$D==1]) # last event time
  if (max(dat$Y) > maxE){
    plat <- max(dat$Y) - maxE
    numBefore <- sum(dat$D==1 & dat$Y > (maxE-plat)) # number of events that are plateau length before the last event
    return(
      (1-numBefore/nrow(dat))**nrow(dat) # alphahat
    )
  } else return(NA) # if the last observation is an event, then this test fails
}

## qn test statistic Maller Zhou (1996)
qnTest <- function(dat){ # expects a df with columns Y,D
  maxE <- max(dat$Y[dat$D==1]) # last event time
  maxAll <- max(dat$Y)
  if (maxAll > maxE){
    numPlat <- sum(dat$D==1 & dat$Y > (2*maxE-maxAll) & dat$Y <= maxE) # number of events that are between (2*maxE-maxAll) (exclusive) and maxE (inclusive)
    return(
      numPlat/nrow(dat) # qn
    )
  } else return(NA) # if the last observation is an event, then this test fails
}

## Shen 2000
shenTest <- function(dat){ 
  maxE <- max(dat$Y[dat$D==1]) # last event time
  maxAll <- max(dat$Y)
  if (maxAll > maxE){
    w <- (maxAll - maxE)/maxAll
    tauG <- w*maxE+(1-w)*maxAll
    numBefore <- sum(dat$D==1 & 
                       (dat$Y >= tauG*maxE/maxAll & dat$Y <= maxE)
    ) # number of events that are "plateau" length before the last event (plateau here is different than before)
    return(
      (1-numBefore/nrow(dat))**nrow(dat) # alphatilde
    )
  } else return(NA) # if the last observation is an event, then this test fails
}

## MZ 1996 Test for Immunes
immuneTest <- function(dat){
  kmfit <- survfit(Surv(dat$Y,dat$D)~1)
  
  lastObs <- max(dat$Y) 
  isCens <- (dat$D[which.max(dat$Y)]==0) # test if the last observation is a censoring, sanity check that phat < 1
  pHat <- 1-summary(kmfit,times=lastObs)$surv
  pCens <- 1-mean(dat$D) # proportion of censored
  
  return(c(pHat,pCens,lastObs,isCens))
}

s0106 <- read.csv("public_s0106.csv")
s1117 <- read.csv("public_s1117.csv")

### S0106: later plateau is not too distinguishable from early plateau 
## Data
dat11 <- s0106[,2:3]
dat18 <- s0106[,4:5]

colnames(dat11) <- colnames(dat18) <- c("Y","D")

dat11$Surv <- with(dat11,Surv(Y,D))
dat18$Surv <- with(dat18,Surv(Y,D))

## Analysis
pdf("S0106_Paper.pdf",width=8,height=5,onefile=FALSE)
plot(survfit(dat18$Surv~1),conf.int=FALSE,col="blue",mark.time=TRUE,lty=2,
     ylab="Survival Probability",main="Kaplan-Meier Estimate of S0106 Overall Survival",xlab="Years since Registration",
     cex.axis=labSize,cex.lab=labSize,cex.main=labSize)
lines(survfit(dat11$Surv~1),conf.int=FALSE,col="black",mark.time=TRUE)
legend("bottomleft",legend=c("2011 Data Release","8 Years of Follow-Up"),lty=1:2,col=c("black","blue"),
       bty="n",cex=labSize)
# dev.off()

dat11 <- data.frame(Y=dat11$Y,D=dat11$D)

dat11$Y[dat11$Y==0] <- 1e-6
dat11$SurvObj <- with(dat11,Surv(Y,D))
# flexsurvcure(SurvObj~1,data=dat11,dist="exp")$AIC
# flexsurvcure(SurvObj~1,data=dat11,dist="weibull")$AIC # best-fitting with AIC
# flexsurvcure(SurvObj~1,data=dat11,dist="gamma")$AIC
# flexsurvcure(SurvObj~1,data=dat11,dist="llogis")$AIC
# 
# flexsurvreg(SurvObj~1,data=dat11,dist="exp")$AIC
# flexsurvreg(SurvObj~1,data=dat11,dist="weibull")$AIC
# flexsurvreg(SurvObj~1,data=dat11,dist="gamma")$AIC
# flexsurvreg(SurvObj~1,data=dat11,dist="llogis")$AIC

ratioOut(dat11,max(dat11$Y),"wei")

mzTest(dat11)
immunOut <- immuneTest(dat11) 
immunOut[1]
round((1/immunOut[2])-1)
nrow(dat11)

shenTest(dat11)
qnTest(dat11)


s0106.weiFit <- flexsurvcure(SurvObj~1,data=dat11,dist="weibull")$res[,1]

curve(s0106.weiFit[1]+(1-s0106.weiFit[1])*
        pweibull(x,shape=s0106.weiFit[2],scale=s0106.weiFit[3],lower.tail=FALSE),
      lty=1,col="black",
      add=TRUE)

dev.off()



### S1117: early plateau but no plateau at all later
## Data

dat14 <- s1117[,2:3]
dat18 <- s1117[,4:5]

colnames(dat14) <- colnames(dat18) <- c("Y","D")

dat14$Surv <- with(dat14,Surv(Y,D))
dat18$Surv <- with(dat18,Surv(Y,D))

## Analysis
pdf("S1117_2014.pdf",width=8,height=5,onefile=FALSE)
plot(survfit(dat14$Surv~1),conf.int=FALSE,col="black",mark.time=TRUE,
     ylab="Survival Probability",main="Kaplan-Meier Estimate of S1117 Overall Survival",xlab="Years since Registration",
     cex.axis=labSize,cex.lab=labSize,cex.main=labSize)
legend("bottomleft",legend=c("2014 Data Release"),lty=1,col=c("black"),
       bty="n",cex=labSize)
dev.off()

pdf("S1117_Paper.pdf",width=8,height=5,onefile=FALSE)
plot(survfit(dat18$Surv~1),conf.int=FALSE,col="blue",mark.time=TRUE,lty=2,
     ylab="Survival Probability",main="Kaplan-Meier Estimate of S1117 Overall Survival",xlab="Years since Registration",
     cex.axis=labSize,cex.lab=labSize,cex.main=labSize)
lines(survfit(dat14$Surv~1),conf.int=FALSE,col="black",mark.time=TRUE)
legend("bottomleft",legend=c("2014 Data Release","5 Years of Follow-Up"),lty=1:2,col=c("black","blue"),
       bty="n",cex=labSize)


dat14 <- data.frame(Y=dat14$Y,D=dat14$D)

dat14$Y[dat14$Y==0] <- 1e-6 # flexsurv does not like exactly 0 times
dat14$SurvObj <- with(dat14,Surv(Y,D))
# flexsurvcure(SurvObj~1,data=dat14,dist="exp")$AIC
# flexsurvcure(SurvObj~1,data=dat14,dist="weibull")$AIC
# flexsurvcure(SurvObj~1,data=dat14,dist="gamma")$AIC 
# flexsurvcure(SurvObj~1,data=dat14,dist="llogis")$AIC
# 
# flexsurvreg(SurvObj~1,data=dat14,dist="exp")$AIC
# flexsurvreg(SurvObj~1,data=dat14,dist="weibull")$AIC
# flexsurvreg(SurvObj~1,data=dat14,dist="gamma")$AIC  
# flexsurvreg(SurvObj~1,data=dat14,dist="llogis")$AIC # best-fitting with AIC

# no need for RECeUS because pihat=0
mzTest(dat14)
immunOut2 <- immuneTest(dat14) 
immunOut2[1]
round((1/immunOut2[2])-1)
nrow(dat14)

shenTest(dat14)
qnTest(dat14)

dev.off()

dat18 <- data.frame(Y=dat18$Y,D=dat18$D)

dat18$Y[dat18$Y==0] <- 1e-6 # flexsurv does not like exactly 0 times
dat18$SurvObj <- with(dat18,Surv(Y,D))

# flexsurvcure(SurvObj~1,data=dat18,dist="exp")$AIC
# flexsurvcure(SurvObj~1,data=dat18,dist="weibull")$AIC
# flexsurvcure(SurvObj~1,data=dat18,dist="gamma")$AIC 
# flexsurvcure(SurvObj~1,data=dat18,dist="llogis")$AIC# best-fitting model
# 
# flexsurvreg(SurvObj~1,data=dat18,dist="exp")$AIC
# flexsurvreg(SurvObj~1,data=dat18,dist="weibull")$AIC
# flexsurvreg(SurvObj~1,data=dat18,dist="gamma")$AIC
# flexsurvreg(SurvObj~1,data=dat18,dist="llogis")$AIC # best-fitting non-cure AIC

ratioOut(dat18,max(dat18$Y),"llogis")

mzTest(dat18)
immunOut3 <- immuneTest(dat18) 
immunOut3[1]
round((1/immunOut3[2])-1)
nrow(dat18)

shenTest(dat18)
qnTest(dat18)
