rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
library(vars)
library(forecast)

# [8] "ESPbirthsRR.txt"              
# [9] "ESPexposRR.txt"               
# [10] "FRATNPasfrRR.txt"             
# [11] "FRATNPexposRR.txt"

## reading data
E0 <- read.table("ESPexposRR.txt", header=TRUE, skip=2)
Y0 <- read.table("ESPbirthsRR.txt", header=TRUE, skip=2)
E0 <- subset(E0, Year>=1960)
Y0 <- subset(Y0, Year>=1960)


## dimensions
x <- 12:55
m <- length(x)
t1 <- unique(E0$Year)
n1 <- length(t1)
t <- min(t1):2050
n <- length(t)
tF <- t[-c(1:n1)]
nF <- length(tF)

## colors
colon1 <- rainbow(n1)
colon <- rainbow(n)
colom <- rainbow(m)

## data in matrix
E1 <- matrix(E0$Exposure, m, n1)
Y1 <- matrix(Y0$Total, m, n1)
## in vector
e1 <- c(E1)
y1 <- c(Y1)

## actual ASFR
MX1 <- Y1/E1
mx1 <- c(MX1)
## in log
lMX1 <- log(MX1)
lmx1 <- c(lMX1)

## plotting data
whi <- floor(seq(1,n1,length=9))
rany <- range(MX1, 0)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(x, MX1[,whi[i]], pch=16, 
       col=colon1[whi[i]], main=paste(t1[whi[i]]),
       ylim=rany)
}
par(mfrow=c(1,1))

## in 3D
xx1 <- rep(x, n1)
tt1 <- rep(t1, each=m)
plot3d(xx1, tt1, mx1, type="p", col="red", 
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)

## compute observed TFR
TFR1 <- colSums(MX1)
plot(t1, TFR1)

## compute mean age at childbearing
MAB1 <- apply(MX1, 2, 
              function(mx) sum(mx*(x+0.5)) / sum(mx))
plot(t1, MAB1)

## compute variance age at childbearing
## compute variance age at childbearing
VAB1 <- numeric(n1)
for(i in 1:n1){
  mx.i <- MX1[,i]
  mab.i <- MAB1[i]
  VAB1[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}
plot(t1, VAB1)


## plotting summary measures
par(mfrow=c(1,3))
plot(t1, TFR1)
plot(t1, MAB1)
plot(t1, VAB1)
par(mfrow=c(1,1))


## augmenting data with arbitrary values for future years
Y <- matrix(10, m, n)
Y[1:m,1:n1] <- Y1
E <- matrix(10, m, n)
E[1:m,1:n1] <- E1
MX <- Y/E
MX[,(n1+1):n] <- NA

## set weights equal to zero for future years
WEI <- matrix(0, m, n)
WEI[,1:n1] <- 1
## set equal to zero when we have no exposures
WEI[E==0] <- 0
## unimportant in fertility, crucial in mortality



## B-splines
library(MortalitySmooth)
## over ages
Bx <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                       ndx=floor(m/4), deg=3)
nbx <- ncol(Bx)
## over all years (observed+forecast)
Bt <- diag(n)
nbt <- ncol(Bt)
## only for the observed years
Bt1 <- diag(n1)
nbt1 <- ncol(Bt1)
## only for the forecast years
BtF <- diag(nF)
nbtF <- ncol(BtF)


## overall basis
B <- kronecker(Bt, Bx)
nb <- ncol(B)
B1 <- kronecker(Bt1, Bx)
nb1 <- ncol(B1)
BF <- kronecker(BtF, Bx)
nbF <- ncol(BF)

## over ages
Bx1 <- kronecker(matrix(1, ncol = nbx, nrow = 1), Bx)
Bx2 <- kronecker(Bx, matrix(1, ncol = nbx, nrow = 1))
RTBx <- Bx1 * Bx2
## over all years (observed+forecast)
BBt1 <- kronecker(matrix(1, ncol = nbt, nrow = 1), Bt)
BBt2 <- kronecker(Bt, matrix(1, ncol = nbt, nrow = 1))
RTBt <- BBt1 * BBt2
## only for the observed years
BBt11 <- kronecker(matrix(1, ncol = nbt1, nrow = 1), Bt1)
BBt21 <- kronecker(Bt1, matrix(1, ncol = nbt1, nrow = 1))
RTBt1 <- BBt11 * BBt21
## only for the forecast years
BBt1F <- kronecker(matrix(1, ncol = nbtF, nrow = 1), BtF)
BBt2F <- kronecker(BtF, matrix(1, ncol = nbtF, nrow = 1))
RTBtF <- BBt1F * BBt2F

## penalty stuff
Dx <- diff(diag(nbx), diff=3)
tDDx <- t(Dx)%*%Dx
Dt <- diff(diag(nbt), diff=2)
tDDt <- t(Dt)%*%Dt
Px <- kronecker(diag(nbt), tDDx)
Pt <- kronecker(tDDt, diag(nbx))
## smoothing parameters
lambdas <- c(10^2, 10^5)
lambda.x <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.x * Px + lambda.t * Pt
PF <- P
# ## smoothing parameters
# v1 <- c(rep(1,n1-1), rep(0,nF-1))
# V1 <- diag(v1)
# tDDt1 <- t(Dt)%*%V1%*%Dt
# Pt1 <- kronecker(tDDt1, diag(nbx))
# vF <- c(rep(0,n1-1), rep(1,nF-1))
# VF <- diag(vF)
# tDDtF <- t(Dt)%*%VF%*%Dt
# PtF <- kronecker(tDDtF, diag(nbx))
# lambda.tF <- 10^-1
# PF <- lambda.x*Px + lambda.t*Pt1 + lambda.tF*PtF



## data in vector
y <- c(Y)
e <- c(E)
wei <- c(WEI)
mx <- c(MX)

## simple smoothing
ETA <- log((Y+1)/(E+1))
for(it in 1:20){
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBx, RTBt, nbx, nbt, WW)
  tBWBpP <- tBWB + P
  tBWz <- MortSmooth_BcoefB(t(Bx), t(Bt), (WW * z))
  betas <- solve(tBWBpP, c(tBWz))
  BETAS <- matrix(betas, nrow = nbx)
  old.ETA <- ETA
  ETA <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  dif.ETA <- max(abs(old.ETA - ETA))
  cat(it, dif.ETA, "\n")
  if(dif.ETA < 1e-4) break
}
BETAS.S <- BETAS
ETA.S <- ETA
MU.S <- exp(ETA)

## compute estimated TFR, MAB and VAB
TFR.S <- colSums(MU.S)
MAB.S <- apply(MU.S, 2, function(mx) sum(mx*(x+0.5)) / sum(mx))
VAB.S <- numeric(n)
for(i in 1:n){
  mx.i <- MU.S[,i]
  mab.i <- MAB.S[i]
  VAB.S[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}



# ranx <- range(t)
# for(i in 1:m){
#   rany <- range(MX[i,], MU.S[i,], na.rm = TRUE)
#   plot(t, MX[i,], xlim=ranx, ylim=rany)
#   lines(t, MU.S[i,], col=2, lwd=2, lty=2)
#   title(main=x[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }
# for(i in 1:n){
#   rany <- range(MX[,i], MU.S[,i], na.rm=TRUE)
#   plot(x, MX[,i], ylim=rany)
#   lines(x, MU.S[,i], col=2, lwd=2, lty=2)
#   title(main=t[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }

## plotting smooth estimates
whi <- floor(seq(1,n,length=9))
rany <- range(MX1, 0)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(1,1,t="n", xlim=range(x), ylim=rany, 
       main=paste(t[whi[i]], 
                  round(TFR.S[whi[i]],2),
                  round(MAB.S[whi[i]],2),
                  round(VAB.S[whi[i]],2)),
       xlab="age", ylab="asfr")
  points(x, MX[,whi[i]], pch=16, ylim=rany)
  lines(x, MU.S[,whi[i]], col=2, lwd=2)
}
par(mfrow=c(1,1))

plot3d(xx1, tt1, mx1, type="p", col="red",
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.S,
          back='line', front='line', col=2, lwd=1, alpha = 0.5)



## forecasting summary measures (somehow: UN, INSEE, time-series, etc)
TFR.ts <- ts(log(TFR1), start = t1[1])
ds <- 0
df.test.TFR <- ur.df(TFR.ts,lags=2,'trend')
if (df.test.TFR@teststat[1] > df.test.TFR@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
TFR.mod <- auto.arima(TFR.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
pred.TFR <- forecast(TFR.mod, h=nF,level=80)
plot(pred.TFR)


MAB.ts <- ts(MAB1, start = t1[1])
ds <- 0
df.test.MAB <- ur.df(MAB.ts,lags=2,'trend')
if (df.test.MAB@teststat[1] > df.test.MAB@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
MAB.mod <- auto.arima(MAB.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
pred.MAB <- forecast(MAB.mod, h=nF,level=80)
plot(pred.MAB)

VAB.ts <- ts(VAB1, start = t1[1])
ds <- 0
df.test.VAB <- ur.df(VAB.ts,lags=2,'trend')
if (df.test.VAB@teststat[1] > df.test.VAB@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
VAB.mod <- auto.arima(VAB.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
pred.VAB <- forecast(VAB.mod, h=nF,level=80)
plot(pred.VAB)


## define three bootstapping matrices
n.sim <- 100
boot.TFR <- matrix(NA,nrow=nF,ncol=n.sim)
boot.MAB <- boot.VAB <- boot.TFR
## bootsrap TFR, MAB and VAB
for(i in 1:n.sim){
  ## generate simulation with bootsrapping
  TFR.sim <- exp(simulate(TFR.mod, nsim=nF,
                          future=TRUE, bootstrap=TRUE))
  MAB.sim <- simulate(MAB.mod, nsim=nF,
                      future=TRUE, bootstrap=TRUE)
  VAB.sim <- simulate(VAB.mod, nsim=nF,
                      future=TRUE, bootstrap=TRUE)
  ## derive the bootsrap values
  boot.TFR[,i] <- TFR.sim
  boot.MAB[,i] <- MAB.sim
  boot.VAB[,i] <- VAB.sim
}

whi.sim <- sample(1:n.sim, 1)
TFR.tar <- as.vector(exp(pred.TFR$mean))#boot.TFR[,whi.sim]#
MAB.tar <- as.vector(pred.MAB$mean)#boot.MAB[,whi.sim]#
VAB.tar <- as.vector(pred.VAB$mean)#boot.VAB[,whi.sim]#

# par(mfrow=c(1,3))
# plot(tF, TFR.tar)
# TFR.tar <- loess(TFR.tar~tF, span = 0.5)$fitted
# lines(tF, TFR.tar, col=2)
# plot(tF, MAB.tar)
# MAB.tar <- loess(MAB.tar~tF, span = 0.5)$fitted
# lines(tF, MAB.tar, col=2)
# plot(tF, VAB.tar)
# VAB.tar <- loess(VAB.tar~tF, span=0.5)$fitted
# lines(tF, VAB.tar, col=2)


par(mfrow=c(1,3))
ranx <- range(t)
rany <- range(TFR1, boot.TFR, TFR.S)
plot(t1, TFR1, xlim=ranx, ylim=rany, t="l", lwd=3)
matlines(tF, boot.TFR, col="grey80", lty=1, lwd=1)
lines(tF, TFR.tar, col=3, lwd=3)
lines(t, TFR.S, col=2, lwd=3)
legend("top", c("Actual", 
                "Smooth+Forecast", 
                "Target",
                "Bootstrapped Instances"), col=c(1,2,3,"grey80"), 
       lty=c(1,1,1,1), lwd=c(3,3,3,1))

rany <- range(MAB1, boot.MAB, MAB.S)
plot(t1, MAB1, xlim=ranx, ylim=rany, t="l", lwd=3)
matlines(tF, boot.MAB, col="grey80", lty=1, lwd=1)
lines(tF, MAB.tar, col=3, lwd=3)
lines(t, MAB.S, col=2, lwd=3)
legend("top", c("Actual", 
                "Smooth+Forecast", 
                "Target",
                "Bootstrapped Instances"), col=c(1,2,3,"grey80"), 
       lty=c(1,1,1,1), lwd=c(3,3,3,1))


rany <- range(VAB1, boot.VAB, VAB.S)
plot(t1, VAB1, xlim=ranx, ylim=rany, t="l", lwd=3)
matlines(tF, boot.VAB, col="grey80", lty=1, lwd=1)
lines(tF, VAB.tar, col=3, lwd=3)
lines(t, VAB.S, col=2, lwd=3)
legend("top", c("Actual", 
                "Smooth+Forecast", 
                "Target",
                "Bootstrapped Instances"), col=c(1,2,3,"grey80"), 
       lty=c(1,1,1,1), lwd=c(3,3,3,1))


par(mfrow=c(1,1))



## CONSTRAINING ## CONSTRAINING ## CONSTRAINING ## CONSTRAINING  
## CONSTRAINING ## CONSTRAINING ## CONSTRAINING ## CONSTRAINING  

ONE <- matrix(0, nb1, nF)
## constraints on TFR over future years
hxFUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  TFR.F <- colSums(MU.F)
  hx <- TFR.F - TFR.tar
  return(hx)
}
hx1FUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  hx1 <- matrix(0, nbF, nF)
  for(j in 1:nF){
    MU.F.j <- rep(0, m*nF)
    MU.F.j[1:m+(j-1)*m] <- MU.F[,j]
    for(i in 1:nbF){
      hx1[i,j] <- sum(BF[,i] * MU.F.j)
    }
  }
  return(hx1)
}

## constrained MAB over time
gxFUN <- function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  x05 <- x+0.5
  MAB.F <- apply(MU.F, 2, function(mx)sum(mx*x05)/sum(mx))
  gx <- MAB.F - MAB.tar
  return(gx)
}

gx1FUN <- function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  x05 <- x+0.5

  den <- colSums(MU.F)
  den2 <- den^2
  MUx.F <- MU.F*x05
  mux.F <- c(MUx.F)
  sum.MUx.F <- colSums(MUx.F)

  gx1 <- matrix(0, nbF, nF)
  for(j in 1:nF){
    whi1 <- rep(0,nF*m)
    whi1[1:m+(j-1)*m] <- 1
    for(i in 1:nbF){
      Bmux.i <- BF[,i]*mux.F
      part1 <- sum(Bmux.i*whi1)/den[j]
      Bmu <- BF[,i]*c(MU.F)
      part2 <- ( sum.MUx.F[j] * sum(Bmu*whi1) ) / den2[j]
      gx1[i,j] <- part1 - part2
    }
  }
  return(gx1)
}

## constrained VAB over time
uxFUN <- function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  x05 <- x+0.5
  
  TFR.F <- colSums(MU.F)
  MAB.F <- apply(MU.F, 2, function(mx)sum(mx*x05)/sum(mx))
  VAB.F <- numeric(nF)
  for(i in 1:nF){
    VAB.F[i] <- sum(((x05 - MAB.F[i])^2)*MU.F[,i])/TFR.F[i]
  }
  ux <- VAB.F - VAB.tar
  return(ux)
}

ux1FUN <- function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  x05 <- x+0.5
  mu.F <- c(MU.F)

  TFR.F <- colSums(MU.F)
  TFR2.F <- TFR.F^2
  MAB.F <- apply(MU.F, 2, function(mx)sum(mx*x05)/sum(mx))

  x.mab <- sweep(matrix(x05,m,nF), 2, MAB.F)
  x.mab2 <- x.mab^2
  x.mab2.MU <- x.mab2*MU.F
  sum.x.mab2.mu <- colSums(x.mab2.MU)
  x.mab2.mu <- c(x.mab2.MU)

  ux1 <- matrix(0, nbF, nF)
  for(j in 1:nF){
    whi1 <- rep(0,nF*m)
    whi1[1:m+(j-1)*m] <- 1
    for(i in 1:nbF){
      B.x.mab2.mu <- BF[,i]*x.mab2.mu
      part1 <- sum(B.x.mab2.mu*whi1) / TFR.F[j]
      part2 <- (sum(BF[,i]*mu.F*whi1) * sum.x.mab2.mu[j])/TFR2.F[j]
      ux1[i,j] <- part1 - part2
    }
  }
  return(ux1)
}


BETAS <- BETAS.S
for(it in 1:50){
  betas <- c(BETAS)
  
  ETA <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBx, RTBt, nbx, nbt, WW)
  tBWBpP <- tBWB + PF
  tBWz <- MortSmooth_BcoefB(t(Bx), t(Bt), (WW * z))
  
  fx1 <- tBWBpP%*%betas - c(tBWz)
  fx2 <- tBWBpP
  
  BETASF <- BETAS[,(nbt1+1):nbt]
  hx <- hxFUN(BETASF)
  hx1 <- hx1FUN(BETASF)
  
  gx <- gxFUN(BETASF)
  gx1 <- gx1FUN(BETASF)
   
  ux <- uxFUN(BETASF)
  ux1 <- ux1FUN(BETASF)
  
  Am <- cbind(rbind(ONE, hx1), 
              rbind(ONE, gx1),
              rbind(ONE, ux1))
  bv <- -c(hx, gx, ux)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(ncol(Am))))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  D <- matrix(d, nbx, nbt)
  BETAS <- BETAS + D
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)

## compute estimated TFR, MAB and VAB
TFR.C <- colSums(MU.C)
MAB.C <- apply(MU.C, 2, function(mx) sum(mx*(x+0.5)) / sum(mx))
VAB.C <- numeric(n)
for(i in 1:n){
  mx.i <- MU.C[,i]
  mab.i <- MAB.C[i]
  VAB.C[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}

## estimated TFR, MAB and VAB
ranx <- range(t)
par(mfrow=c(1,3))
rany <- range(TFR1, TFR.tar, TFR.S, TFR.C)
plot(t1, TFR1, xlim=ranx, ylim=rany, pch=16, col=1)
lines(t, TFR.C, col=4, lwd=4)
lines(tF, TFR.tar, col=3, lwd=2, lty=3)
lines(t, TFR.S, col=2, lwd=2, lty=2)
points(t1, TFR1, pch=16, col=1)
abline(v=t[n1]+0.5)
legend("top", c("Actual",
                "Smooth+Forecast", 
                "Target Future TFR",
                "Smooth+Forecast+Constrained"), 
       col=c(1,2,3,4), pch=c(16,NA,NA,NA),
       lty=c(NA, 2, 3, 1), lwd=c(1,2,2,4))

rany <- range(MAB1, MAB.tar, MAB.S, MAB.C)
plot(t1, MAB1, xlim=ranx, ylim=rany, pch=16, col=1)
lines(t, MAB.C, col=4, lwd=4)
lines(tF, MAB.tar, col=3, lwd=2, lty=3)
lines(t, MAB.S, col=2, lwd=2, lty=2)
points(t1, MAB1, pch=16, col=1)
abline(v=t[n1]+0.5)
legend("top", c("Actual",
                "Smooth+Forecast", 
                "Target Future MAB",
                "Smooth+Forecast+Constrained"), 
       col=c(1,2,3,4), pch=c(16,NA,NA,NA),
       lty=c(NA, 2, 3, 1), lwd=c(1,2,2,4))

rany <- range(VAB1, VAB.tar, VAB.S, VAB.C)
plot(t1, VAB1, xlim=ranx, ylim=rany, pch=16, col=1)
lines(t, VAB.C, col=4, lwd=4)
lines(tF, VAB.tar, col=3, lwd=2, lty=3)
lines(t, VAB.S, col=2, lwd=2, lty=2)
points(t1, VAB1, pch=16, col=1)
abline(v=t[n1]+0.5)
legend("top", c("Actual",
                "Smooth+Forecast", 
                "Target Future VAB",
                "Smooth+Forecast+Constrained"), 
       col=c(1,2,3,4), pch=c(16,NA,NA,NA),
       lty=c(NA, 2, 3, 1), lwd=c(1,2,2,4))

par(mfrow=c(1,1))






# ranx <- range(t)
# for(i in 1:m){
#   rany <- range(MX[i,], MU.S[i,], MU.C[i,], na.rm = TRUE)
#   plot(t, MX[i,], xlim=ranx, ylim=rany)
#   lines(t, MU.S[i,], col=2, lwd=2, lty=2)
#   lines(t, MU.C[i,], col=4, lwd=3)
#   title(main=x[i])
#   locator(1)
#   Sys.sleep(0.1)
# }
# for(i in 1:n){
#   rany <- range(MX[,i], MU.S[,i], MU.C[,i], na.rm=TRUE)
#   plot(x, MX[,i], ylim=rany)
#   lines(x, MU.S[,i], col=2, lwd=2, lty=2)
#   lines(x, MU.C[,i], col=4, lwd=3)
#   title(main=t[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }
matplot(x, MU.C, t="l", lty=1, col=colon)
## plotting smooth estimates
whi <- floor(seq(1,n,length=9))
rany <- range(MX1, MU.S, MU.C, 0)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(1,1,t="n", xlim=range(x), ylim=rany, 
       main=paste(t[whi[i]], 
                  round(TFR.C[whi[i]],2),
                  round(MAB.C[whi[i]],2),
                  round(VAB.C[whi[i]],2)),
       xlab="age", ylab="asfr")
  points(x, MX[,whi[i]], pch=16, ylim=rany)
  lines(x, MU.S[,whi[i]], col=2, lwd=2, lty=2)
  lines(x, MU.C[,whi[i]], col=4, lwd=3)
}
par(mfrow=c(1,1))

whi <- floor(seq(1,m,length=9))
ranx <- range(t)
par(mfrow=c(3,3))
for(i in 1:9){
  rany <- range(MX1[whi[i],], 
                MU.S[whi[i],], 
                MU.C[whi[i],])
  plot(1,1,t="n", xlim=ranx, ylim=rany, 
       main=paste(x[whi[i]]),
       xlab="year", ylab="asfr")
  points(t, MX[whi[i],], pch=16)
  lines(t, MU.S[whi[i],], col=2, lwd=2, lty=2)
  lines(t, MU.C[whi[i],], col=4, lwd=3)
}
par(mfrow=c(1,1))

whi <- floor(seq(1,m,length=9))
ranx <- range(t)
par(mfrow=c(3,3))
for(i in 1:9){
  rany <- range(lMX1[whi[i],], 
                ETA.S[whi[i],], 
                ETA.C[whi[i],], finite=TRUE,
                na.rm = TRUE)
  plot(1,1,t="n", xlim=ranx, ylim=rany, 
       main=paste(x[whi[i]]),
       xlab="year", ylab="log-asfr")
  points(t1, lMX1[whi[i],], pch=16)
  lines(t, ETA.S[whi[i],], col=2, lwd=2, lty=2)
  lines(t, ETA.C[whi[i],], col=4, lwd=3)
}
par(mfrow=c(1,1))


plot3d(xx1, tt1, mx1, type="p", col="red",
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.S,
          back='line', front='line', col=2, lwd=1, alpha = 0.5)
surface3d(x, t, MU.C,
          back='line', front='line', col=4, lwd=2, alpha = 0.5)


save.image("ESPfert.RData")














# ## first trial: enforce the sum of the linear predictor
# ## for each year to be equal to a given vector
# 
# a <- round(colSums(ETA.T))
# A <- kronecker(diag(n), matrix(1,1,m))
# AB <- A%*%B
# 
# ## test constraints with starting values
# cbind(AB%*%betas.S, a)
# 
# eta <- eta.S
# for(it in 1:100){
#   mu <- e*exp(eta)
#   z <- (y - mu)/mu + eta
#   w <- c(mu)
#   tBWB <- t(B) %*% (w * B)
#   tBWBpP <- tBWB + P 
#   tBWz <- t(B) %*% (w * z)
#   
#   LHS <- rbind(cbind(tBWBpP, t(AB)),
#                cbind(AB, 0*diag(n)))
#   RHS <- c(tBWz, a)
#   
#   betas <- solve(LHS, RHS)[1:nb]
#   old.eta <- eta
#   eta <- B %*% betas
#   dif.eta <- max(abs(old.eta - eta))
#   cat(it, dif.eta, "\n")
#   if(dif.eta < 1e-4) break
# }
# betas.C <- betas
# BETAS.C <- matrix(betas.C, nbx, nbt)
# eta.C <- eta
# mu.C <- exp(eta.C)
# ETA.C <- matrix(eta.C, m, n)
# MU.C <- matrix(mu.C, m, n)  
# 
# cbind(AB%*%betas.S, AB%*%betas.C, a)
# 
# 
# TFR.C <- colSums(MU.C)
# plot(t, TFR.S)
# points(t, TFR, col=2, pch=3)
# points(t, TFR.C, col=3, pch=4)
# 
# 
# plot3d(xx, tt, mx, type="p", col="red", 
#        xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
# surface3d(x, t, MU.T, 
#           back='line', front='line', 
#           col=5, lwd=1, alpha = 0.5)
# surface3d(x, t, MU.S, 
#           back='line', front='line', 
#           col=2, lwd=1, alpha = 0.5)
# surface3d(x, t, MU.C, 
#           back='line', front='line', 
#           col=3, lwd=2, alpha = 0.5)
# 
# 
# # ranx <- range(t)
# # for(i in 1:m){
# #   rany <- range(MX[i,], MU.S[i,])
# #   plot(t, MX[i,], xlim=ranx, ylim=rany)
# #   lines(t, MU.S[i,], col=2, lwd=2)
# #   lines(t, MU.T[i,], col=4, lty=2)
# #   lines(t, MU.C[i,], col=3, lty=3, lwd=3)
# #   title(main=x[i])
# #   # locator(1)
# #   Sys.sleep(0.1)
# # }
# # for(i in 1:n){
# #   rany <- range(MX[,i], MU.S[,i])
# #   plot(x, MX[,i], ylim=rany)
# #   lines(x, MU.S[,i], col=2, lwd=2)
# #   lines(x, MU.T[,i], col=4, lty=2)
# #   lines(x, MU.C[,i], col=3, lty=3, lwd=3)
# #   title(main=t[i])
# #   # locator(1)
# #   Sys.sleep(0.1)
# # }


## second trial: enforce the sum of the expected values
## for each year to be equal to a given vector, TFR




  
  
##################################################


## Poisson-likelihood example with fertility law, two years,
## with B-splines
## Poisson-likelihood example with fertility law, two years,
## with B-splines
## Poisson-likelihood example with fertility law, two years,
## with B-splines
## Poisson-likelihood example with fertility law, two years,
## with B-splines

rm(list = ls())
options(device="X11")
library(quadprog)
library(magic)

Hadwiger <- function(a, b, c, x){
  ## function for computing asfr from
  ## Hadwiger (1940) model
  p1 <- a*b/c * (c/x)^(3/2)
  p2 <- -b^2 * (c/x + x/c -2)
  asfr <- p1*exp(p2)
  return(asfr)
}

xx <- 12:55
m <- length(xx)

# a, b and c are associated with total fertility, height of
# the curve and Mean Age at Childbearing.
## true parameters
a1T <- 1.75
b1T <- 3
c1T <- 28
a2T <- 1.7
b2T <- 3.1
c2T <- 29

mu1T <- Hadwiger(a1T, b1T, c1T, x=xx)
eta1T <- log(mu1T)
mu2T <- Hadwiger(a2T, b2T, c2T, x=xx)
eta2T <- log(mu2T)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e1 <- e/20
e2 <- e/19

y1T <- e1*mu1T
y1 <- rpois(m, y1T)
mx1 <- y1/e1
lmx1 <- log(mx1)
y2T <- e2*mu2T
y2 <- rpois(m, y2T)
mx2 <- y2/e2
lmx2 <- log(mx2)

plot(xx, mx1, t="p", col=1, 
     ylim=range(0,mx1), main="asfr")
lines(xx, mu1T, col=1, lwd=2)
points(xx, mx2, t="p", col=2)
lines(xx, mu2T, col=2, lwd=2)

## observed TFR
TFR1obs <- sum(mx1)
TFR2obs <- sum(mx2)
## target TFR
TFR1tar <- 3.9
TFR2tar <- 4

## B-splines
B1 <- MortSmooth_bbase(x=xx,xl=min(xx),xr=max(xx),ndx=floor(m/4),deg=3)
B2 <- B1
B <- adiag(B1, B2)
nb <- ncol(B)
## penalty stuff
D <- diff(diag(nb/2), diff=3)
tDD <- t(D)%*%D
lambda <- 10^2
P <- adiag(lambda*tDD, lambda*tDD)

## combining two years
y <- c(y1, y2)
e <- c(e1, e2)


eta <- log((y+1)/(e+1))
betas <- solve(t(B)%*%B, t(B)%*%eta)
for(it in 1:200){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  d <- solve(fx2, -fx1)
  betas <- betas + d
  if(max(abs(d))<1e-4) break
  cat(it, max(abs(d)), "\n")
}
## only smooth
betasS <- betas
etaS <- eta
muS <- mu
mx1S <- mu[1:m]/e1
mx2S <- mu[1:m+m]/e2

TFR1smo <- sum(mx1S)
TFR2smo <- sum(mx2S)

data.frame("years1"=c("smo"=TFR1smo, "obs"=TFR1obs, "tar"=TFR1tar),
           "years2"=c("smo"=TFR2smo, "obs"=TFR2obs, "tar"=TFR2tar))

hxFUN.1 <-  function(betas){
  betas.1 <- betas[1:(nb/2)]
  eta.1 <- B1%*%betas.1
  exp.eta.1 <- exp(eta.1)
  sum.exp.eta.1 <- sum(exp.eta.1)
  hx.1 <- sum.exp.eta.1 - TFR1tar # !!!
  return(hx.1)
}
hxFUN.2 <-  function(betas){
  betas.2 <- betas[1:(nb/2)+(nb/2)]
  eta.2 <- B2%*%betas.2
  exp.eta.2 <- exp(eta.2)
  sum.exp.eta.2 <- sum(exp.eta.2)
  hx.2 <- sum.exp.eta.2 - TFR2tar # !!!
  return(hx.2)
}

hxFUN.1(betas);TFR1smo-TFR1tar
hxFUN.2(betas);TFR2smo-TFR2tar

hx1FUN.1 <- function(betas){
  betas.1 <- betas[1:(nb/2)]
  exp.eta.1 <- exp(B1%*%betas.1)
  hx1.1 <- sum(B1[,1]*exp.eta.1)
  for(i in 2:(nb/2)){
    hx1.1 <- rbind(hx1.1, sum(B1[,i]*exp.eta.1))
  }
  return(hx1.1)
}
hx1FUN.1(betas)

hx1FUN.2 <- function(betas){
  betas.2 <- betas[1:(nb/2)+(nb/2)]
  exp.eta.2 <- exp(B2%*%betas.2)
  hx1.2 <- sum(B2[,1]*exp.eta.2)
  for(i in 2:(nb/2)){
    hx1.2 <- rbind(hx1.2, sum(B2[,i]*exp.eta.2))
  }
  return(hx1.2)
}
hx1FUN.2(betas)


betas <- betasS
for(it in 1:20){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  
  hx.1 <- hxFUN.1(betas)
  hx1.1 <- hx1FUN.1(betas)
  Am.1 <- hx1.1
  bv.1 <- -hx.1
  
  hx.2 <- hxFUN.2(betas)
  hx1.2 <- hx1FUN.2(betas)
  Am.2 <- hx1.2
  bv.2 <- -hx.2
  
  Am <- adiag(Am.1, Am.2)
  bv <- c(bv.1, bv.2)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), diag()))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(betas);hxFUN(betasS)

TFRsmo <- sum(exp(etaS))
TFRcon <- sum(exp(eta))
cbind(TFRobs, TFRsmo, TFRcon)

par(mfrow=c(2,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(etaS), col=3, lwd=2, lty=2)
lines(xx, e*exp(eta), col=4, lwd=2, lty=1)

plot(xx, y/e, t="h", lwd=2, 
     ylim=range(0, y/e, exp(etaS), exp(eta)), main="rates")
lines(xx, yT/e, col=2, lwd=2)
lines(xx, exp(etaS), col=3, lwd=2, lty=2)
lines(xx, exp(eta), col=4, lwd=2, lty=1)

plot(xx, log(y/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))
# yrsUN <- 2016:2100
# # ## Spain
# # TFRmed0 <- c(1.39,1.44,1.48,1.52,1.55,1.58,1.61,1.63,1.64,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.72)
# # TFRhigh0 <- c(1.64,1.84,1.98,2.02,2.05,2.08,2.11,2.13,2.14,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.22)
# # TFRlow0 <- c(1.14,1.04,0.98,1.02,1.05,1.08,1.11,1.13,1.14,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.22)
# # MABun0 <- c(31.92,32.17,32.35,32.48,32.55,32.60,32.61,32.61,32.61,32.61,32.61,32.61,32.61,32.61,32.61,32.61,32.61)
# ## France
# TFRmed0 <- c(1.97,1.97,1.96,1.96,1.95,1.95,1.95,1.95,1.95,1.95,1.94,1.94,1.94,1.94,1.94,1.94,1.94)
# TFRhigh0 <- c(2.22,2.37,2.46,2.46,2.45,2.45,2.45,2.45,2.45,2.45,2.44,2.44,2.44,2.44,2.44,2.44,2.44)
# TFRlow0 <- c(1.72,1.57,1.46,1.46,1.45,1.45,1.45,1.45,1.45,1.45,1.44,1.44,1.44,1.44,1.44,1.44,1.44)
# MABun0 <- c(30.34,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54,30.54)
# 
# TFRmed1 <- rep(TFRmed0, each=5)
# TFRhigh1 <- rep(TFRhigh0, each=5)
# TFRlow1 <- rep(TFRlow0, each=5)
# MABun1 <- rep(MABun0, each=5)
# 
# rany <- range(TFRmed1, TFRlow1, TFRhigh1)
# plot(yrsUN, TFRmed1, ylim=rany)
# points(yrsUN, TFRlow1, col=2)
# points(yrsUN, TFRhigh1, col=3)
