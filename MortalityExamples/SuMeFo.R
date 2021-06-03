rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(quadprog)
library(rgl)
library(MortalitySmooth)
library(HMDdata)
library(HMDHFDplus)
library(vars)
library(forecast)
library(colorspace)
library(RColorBrewer)
library(fields)

## ages
ages <- 0:104
m <- length(ages)
## observed years
years1 <- 1960:2018
n1 <- length(years1)
## all years: observed+future
years <- years1[1]:2050
n <- length(years)
## forecast years
yearsF <- (years1[n1]+1):2050
nF <- length(yearsF)
## domains for the estimation
a <- 1:m-1
t1 <- 1:n1-1
t <- 1:n-1
tF <- t[-c(1:n1)]

## colors
coln <- rainbow_hcl(n)
colm <- rainbow_hcl(m)
colnT <- adjustcolor(coln, 0.5)
colmT <- adjustcolor(colm, 0.5)
# ## USING readHMDweb() for the final version
# ## select population
# pop <- "USA"
# ## select sex
# sex <- "Male" ## alternative "Female"
# ## loading exposures from HMD 
# Ehmd <- readHMDweb(CNTRY = pop, 
#                    item = "Exposures_1x1", 
#                    fixup = TRUE)
# ## loading deaths from HMD 
# Yhmd <- readHMDweb(CNTRY = pop, 
#                    item = "Deaths_1x1",
#                    fixup = TRUE)
# ## in both two previous steps
# ## if arguments "username" and "password" are not
# ## provided you will be prompted about  
# ## selecting ages and years
# E0 <- subset(Ehmd, Year%in%years1 & Age%in%ages, sex)
# ## place in a mXn1 matrix
# E1 <- matrix(E0[,1], m, n1)
# ## selecting ages and years
# Y0 <- subset(Yhmd, Year%in%years1 & Age%in%ages, sex)
# ## place in a mXn1 matrix
# Y1 <- matrix(Y0[,1], m, n1)
# ## observed log-rates
# MX1 <- Y1/E1
# ETA1 <- log(MX1)

## USING HMDdata for the current version
## select population
pop <- "Scotland" ## alternative names(HMDdata)
## select sex
sex <- "M" ## alternative "F"
## extract deaths and exposures from HMDdata
Y1 <- matrix(selectHMDdata(pop, "Deaths", sex, ages, years1), m, n1)
E1 <- matrix(selectHMDdata(pop, "Exposures", sex, ages, years1), m, n1)
MX1 <- Y1/E1
ETA1 <- log(MX1)

## plotting
matplot(ages, ETA1, col=colnT[1:n1], t="l", lty=3, 
        xlab="ages", ylab="log-mortality")
whi <- floor(seq(1,n1,length=4))
for(i in whi) lines(ages, ETA1[,i], col=coln[i], lwd=4)
legend("topleft", inset=0.1, legend=years1[whi], col=coln[whi], lwd=4)

## observed e0
e01 <- apply(MX1, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
## plotting
plot(years1, e01, pch=16, xlab="years", ylab="e0")

## observed e-dagger
ed1 <-  apply(MX1, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
## plotting
plot(years1, ed1, pch=16, xlab="years", ylab="e-dagger")

## gini coefficients
gi1 <- apply(MX1, 2, function(mx) 1 - sum(exp(-cumsum(mx))^2)/(sum(exp(-cumsum(mx))))-.0006)
## plotting
plot(years1, gi1, pch=16, xlab="years", ylab="Gini coefficient")


# ## computing e0 using matrix algebra
# MX1vec <- c(MX1)
# In <- diag(n1)
# Im <- diag(m)
# C <- lower.tri(Im, diag = TRUE)
# C[C==1] <- -1
# CC <- kronecker(In, C)
# ones.m <- matrix(1,1,m)
# Iones <- kronecker(In, ones.m)
# e01a <- Iones %*% exp(CC%*%MX1vec) + 0.5


## matrices/vectors arbitrary values
Y <- matrix(10, m, n)
Y[1:m,1:n1] <- Y1
E <- matrix(10, m, n)
E[1:m,1:n1] <- E1
MX <- Y/E
MX[,(n1+1):n] <- NA
ETA <- log(MX)

## set weights equal to zero for future years
WEI <- matrix(0, m, n)
WEI[,1:n1] <- 1
## set equal to zero when we have no exposures
WEI[E1==0] <- 0


# ## PLAIN P-splines
# 
# ## B-splines
# ## over ages
# Ba <- MortSmooth_bbase(x=a, xl=min(a), xr=max(a), ndx=floor(m/5), deg=3)
# nba <- ncol(Ba)
# ## over all years (observed+forecast)
# Bt <- diag(n)
# nbt <- ncol(Bt)
# ## only for the observed years
# Bt1 <- diag(n1)
# nbt1 <- ncol(Bt1)
# ## only for the forecast years
# BtF <- diag(nF)
# nbtF <- ncol(BtF)
# 
# ## overall basis
# B <- kronecker(Bt, Ba)
# nb <- ncol(B)
# B1 <- kronecker(Bt1, Ba)
# nb1 <- ncol(B1)
# BF <- kronecker(BtF, Ba)
# nbF <- ncol(BF)
# 
# ## over ages
# Ba1 <- kronecker(matrix(1, ncol = nba, nrow = 1), Ba)
# Ba2 <- kronecker(Ba, matrix(1, ncol = nba, nrow = 1))
# RTBa <- Ba1 * Ba2
# ## over all years (observed+forecast)
# BBt1 <- kronecker(matrix(1, ncol = nbt, nrow = 1), Bt)
# BBt2 <- kronecker(Bt, matrix(1, ncol = nbt, nrow = 1))
# RTBt <- BBt1 * BBt2
# ## only for the observed years
# BBt11 <- kronecker(matrix(1, ncol = nbt1, nrow = 1), Bt1)
# BBt21 <- kronecker(Bt1, matrix(1, ncol = nbt1, nrow = 1))
# RTBt1 <- BBt11 * BBt21
# ## only for the forecast years
# BBt1F <- kronecker(matrix(1, ncol = nbtF, nrow = 1), BtF)
# BBt2F <- kronecker(BtF, matrix(1, ncol = nbtF, nrow = 1))
# RTBtF <- BBt1F * BBt2F
# 
# 
# ## penalty stuff
# Da <- diff(diag(nba), diff=2)
# tDDa <- t(Da)%*%Da
# Dt <- diff(diag(nbt), diff=2)
# tDDt <- t(Dt)%*%Dt
# Pa <- kronecker(diag(nbt), tDDa)
# Pt <- kronecker(tDDt, diag(nba))
# ## smoothing parameters
# lambdas <- c(10^0, 10^7)
# lambda.a <- lambdas[1]
# lambda.t <- lambdas[2]
# P <- lambda.a * Pa + lambda.t * Pt
# 
# ## simple smoothing
# ETA <- log((Y+1)/(E+1))
# for(it in 1:20){
#   MU <- E*exp(ETA)
#   W <- c(WEI*MU)
#   z <- ETA + (1/MU) * (Y - MU)
#   z[which(WEI == 0)] <- 0
#   WW <- WEI * W
#   tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
#   tBWBpP <- tBWB + P
#   tBWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
#   betas <- solve(tBWBpP, c(tBWz))
#   BETAS <- matrix(betas, nrow = nba)
#   old.ETA <- ETA
#   ETA <- MortSmooth_BcoefB(Ba, Bt, BETAS)
#   dif.ETA <- max(abs(old.ETA - ETA))
#   cat(it, dif.ETA, "\n")
#   if(dif.ETA < 1e-4) break
# }
# BETAS.S <- BETAS
# ETA.S <- ETA
# MU.S <- exp(ETA)
# 
# 
# ## plotting
# rany <- range(ETA1, ETA.S, finite=TRUE)
# matplot(ages, ETA1, col=colnT, t="l", lty=3, ylim=rany,
#         xlab="ages", ylab="log-mortality")
# whi <- floor(seq(1,n,length=6))
# for(i in whi) lines(ages, ETA.S[,i], col=coln[i], lwd=4)
# legend("topleft", inset=0.1, legend=years[whi], col=coln[whi], lwd=4)
# 
# par(mfrow=c(2,3))
# whi <- floor(seq(1,m,length=6))
# ranx <- range(years)
# for(i in whi){
#   rany <- range(ETA1[i,], ETA.S[i,], finite=TRUE)
#   plot(years, ETA.S[i,], t="l", lwd=4, col=colmT[i], main=paste(ages[i]),
#        ylim=rany, xlim=ranx,xlab="ages", ylab="log-mortality")
#   points(years1, ETA1[i,], col=colm[i], pch=16)
#   abline(v=years[n1]+0.5)
# }
# par(mfrow=c(1,1))
# 
# ## compute estimated+forecast e0
# e0.S <- apply(MU.S, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
# 
# ## plotting
# ranx <- range(years)
# rany <- range(e01, e0.S)
# plot(years, e0.S, t="l", col=2, lwd=4,
#      xlab="years", ylab="e0", xlim=ranx, ylim=rany)
# points(years1, e01, pch=16)
# abline(v=years[n1]+0.5)
# 
# ## compute estimated+forecast e-dagger
# ed.S <-  apply(MU.S, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
# 
# ## plotting
# ranx <- range(years)
# rany <- range(ed1, ed.S)
# plot(years, ed.S, t="l", col=2, lwd=4,
#      xlab="years", ylab="e-dagger", xlim=ranx, ylim=rany)
# points(years1, ed1, pch=16)
# abline(v=years[n1]+0.5)


## P-SPLINES + infant mortality coefficients
## over ages w/o age 0
a0 <- a[-1]
m0 <- m-1
Ba0 <- MortSmooth_bbase(x=a0, xl=min(a0), xr=max(a0), ndx=floor(m0/5), deg=3)
nba0 <- ncol(Ba0)
## adding infant-specific basis
Ba <- cbind(0, Ba0)
Ba <- rbind(c(1, rep(0,nba0)), Ba)
nba <- ncol(Ba)

matplot(a, Ba, t="l", lty=1, col=rainbow(nba))
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
B <- kronecker(Bt, Ba)
nb <- ncol(B)
B1 <- kronecker(Bt1, Ba)
nb1 <- ncol(B1)
BF <- kronecker(BtF, Ba)
nbF <- ncol(BF)

## over ages
Ba1 <- kronecker(matrix(1, ncol = nba, nrow = 1), Ba)
Ba2 <- kronecker(Ba, matrix(1, ncol = nba, nrow = 1))
RTBa <- Ba1 * Ba2
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
Da <- diff(diag(nba), diff=2)
Da[1,1] <- 0
tDDa <- t(Da)%*%Da
Dt <- diff(diag(nbt), diff=2)
tDDt <- t(Dt)%*%Dt
Pa <- kronecker(diag(nbt), tDDa)
Pt <- kronecker(tDDt, diag(nba))
## smoothing parameters
lambdas <- c(10^0, 10^7)
lambda.a <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.a * Pa + lambda.t * Pt

## simple smoothing
ETA <- log((Y+1)/(E+1))
for(it in 1:20){
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
  tBWBpP <- tBWB + P
  tBWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
  betas <- solve(tBWBpP, c(tBWz))
  BETAS <- matrix(betas, nrow = nba)
  old.ETA <- ETA
  ETA <- MortSmooth_BcoefB(Ba, Bt, BETAS)
  dif.ETA <- max(abs(old.ETA - ETA))
  cat(it, dif.ETA, "\n")
  if(dif.ETA < 1e-4) break
}
BETAS.SI <- BETAS
ETA.SI <- ETA
MU.SI <- exp(ETA)


## plotting
rany <- range(ETA1, ETA.SI, finite=TRUE)
matplot(ages, ETA1, col=colnT, t="l", lty=3, ylim=rany,
        xlab="ages", ylab="log-mortality")
whi <- floor(seq(1,n,length=6))
for(i in whi) lines(ages, ETA.SI[,i], col=coln[i], lwd=4)
legend("topleft", inset=0.1, legend=years[whi], col=coln[whi], lwd=4)


par(mfrow=c(2,5))
whi <- floor(seq(1,m,length=10))
ranx <- range(years)
for(i in whi){
  rany <- range(ETA1[i,], ETA.SI[i,], finite=TRUE)
  plot(years, ETA.SI[i,], t="l", lwd=4, col=colmT[i], main=paste(ages[i]),
       ylim=rany, xlim=ranx,xlab="ages", ylab="log-mortality")
  points(years1, ETA1[i,], col=colm[i], pch=16)
  abline(v=years[n1]+0.5)
}
par(mfrow=c(1,1))

## compute estimated+forecast e0
e0.SI <- apply(MU.SI, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)



## compute estimated+forecast e-dagger
ed.SI <-  apply(MU.SI, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))



## useful objects
ones.m <- matrix(1,m,1)
Im <- diag(m)
C <- lower.tri(Im, diag = TRUE)
C[C==1] <- -1

## derivative of e0 with respect to age for each year (observed+forecast)
der.e0.SI <- matrix(0, m, n)
for(i in 1:n){
  Ba.betas <- ETA.SI[,i]
  exp.Ba.betas <- exp(Ba.betas)
  Q <- C %*% exp.Ba.betas
  exp.Q <- exp(Q)
  d.exp.Q <- diag(c(exp.Q))
  d.exp.Ba.betas <- diag(c(exp.Ba.betas)) 
  deriv.coef <- t(ones.m) %*% d.exp.Q %*% C %*% d.exp.Ba.betas %*% Ba
  der.e0.SI[,i] <- Ba%*%c(deriv.coef)
}

## derivative of e-dagger with respect to age for each year (observed+forecast)
der.ed.SI <- matrix(0, m, n)
for(i in 1:n){
  Ba.betas <- ETA.SI[,i]
  exp.Ba.betas <- exp(Ba.betas)
  Q <- C %*% exp.Ba.betas
  exp.Q <- exp(Q)
  p1 <- t(Ba) %*% ( t(C) %*% (Q * exp.Q) * exp.Ba.betas)
  p2 <- t(Ba) %*% ( (t(C)%*%exp.Q) * exp.Ba.betas)
  deriv.coef <- -(p1 + p2)
  der.ed.SI[,i] <- Ba%*%c(deriv.coef)
}


## plotting
par(mfrow=c(2,2))
ranx <- range(years)
rany <- range(e01, e0.SI)
plot(years, e0.SI, t="l", col=2, lwd=4,
     xlab="years", ylab="e0", xlim=ranx, ylim=rany)
points(years1, e01, pch=16)
abline(v=years[n1]+0.5)
ranx <- range(years)
rany <- range(ed1, ed.SI)
plot(years, ed.SI, t="l", col=2, lwd=4,
     xlab="years", ylab="e-dagger", xlim=ranx, ylim=rany)
points(years1, ed1, pch=16)
abline(v=years[n1]+0.5)
matplot(ages, der.e0.SI, t="l", col=coln, lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
matplot(ages, der.ed.SI, t="l", col=coln, lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
par(mfrow=c(1,1))


par(mfrow=c(4,2))
whi <- floor(seq(1,m,length=4))
ranx <- range(years)
for(i in whi){
  rany <- range(ETA1[i,], ETA.SI[i,], finite=TRUE)
  plot(years, ETA.SI[i,], t="l", lwd=4, col=colmT[i], main=paste(ages[i]),
       ylim=rany, xlim=ranx,xlab="ages", ylab="log-mortality")
  points(years1, ETA1[i,], col=colm[i], pch=16)
  
  rany <- range(der.e0.SI[i,])
  plot(years, der.e0.SI[i,], t="l", lwd=4, col=colm[i], main=paste(ages[i]),
       ylim=rany, xlim=ranx,xlab="ages", ylab="d e0 / d mu")
  abline(v=years[n1]+0.5)
}
par(mfrow=c(1,1))

##
quantile(der.e0.SI)
brek <- seq(-2.25, 0, 0.25)#c(-2.1, -1.5, -1.25, -1, -0.75, -0.5, -0.25, -0.1, 0)
colo <- rev(brewer.pal(9, "Blues"))
image.plot(years, ages, t(der.e0.SI), breaks=brek, col=colo)
abline(v=years[n1]+0.5)

quantile(der.ed.SI)
brek <- seq(-1, 1.75, 0.25)
colo <- c(rev(brewer.pal(4, "Blues")), brewer.pal(7, "Reds"))
image.plot(years, ages, t(der.ed.SI), breaks=brek, col=colo)
abline(v=years[n1]+0.5)



## forecasting e0
e0.ts <- ts(e01, start = t1[1])
ds <- 0
df.test.e0 <- ur.df(e0.ts,lags=2,'trend')
if (df.test.e0@teststat[1] > df.test.e0@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
e0.mod <- auto.arima(e0.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
pred.e0 <- forecast(e0.mod, h=nF,level=80)
## target e0
e0.tar <- as.vector(pred.e0$mean)


## forecasting e-dagger
ed.ts <- ts(ed1, start = t1[1])
# ## option 1
# ds <- 0
# df.test.ed <- ur.df(ed.ts,lags=2,'trend')
# if (df.test.ed@teststat[1] > df.test.ed@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
# ed.mod <- auto.arima(ed.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
# pred.ed <- forecast(ed.mod, h=nF,level=80)
## option 2
ed.mod <- auto.arima(ed.ts, d=2, max.p=2, max.q=2)
pred.ed <- forecast(ed.mod, h=nF,level=95)
## target ed
ed.tar <- as.vector(pred.ed$mean)

par(mfrow=c(1,2))
plot(pred.e0)
lines(t, e0.SI, col=2, lwd=2, lty=2)
plot(pred.ed)
lines(t, ed.SI, col=2, lwd=2, lty=2)
par(mfrow=c(1,1))

## useful objects
ONE <- matrix(0, nb1, nF)
BETAS <- BETAS.SI
BETASF <- BETAS[,(nbt1+1):nbt]

## constraints on e0 over "future" years
Constr.e0 <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Ba, BtF, BETASF)
  MU.F <- exp(ETA.F)
  e0.F <- apply(MU.F, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
  delta.e0 <- e0.F - e0.tar
  deriv.e0 <- matrix(0, nbF, nF)
  for(i in 1:nF){
    Ba.betas <- ETA.F[,i]
    exp.Ba.betas <- exp(Ba.betas)
    Q <- C %*% exp.Ba.betas
    exp.Q <- exp(Q)
    d.exp.Q <- diag(c(exp.Q))
    d.exp.Ba.betas <- diag(c(exp.Ba.betas)) 
    wr <- 1:nba + (i-1)*nba
    deriv.e0[wr,i] <- t(ones.m) %*% d.exp.Q %*% C %*% d.exp.Ba.betas %*% Ba
  }
  out.e0 <- list(delta.e0=delta.e0,
                 deriv.e0=deriv.e0)
  return(out.e0)
}
Constr.ed <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Ba, BtF, BETASF)
  MU.F <- exp(ETA.F)
  ed.F <- apply(MU.F, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
  delta.ed <- ed.F - ed.tar
  deriv.ed <- matrix(0, nbF, nF)
  for(i in 1:nF){
    Ba.betas <- ETA.F[,i]
    exp.Ba.betas <- exp(Ba.betas)
    Q <- C %*% exp.Ba.betas
    exp.Q <- exp(Q)
    p1 <- t(Ba) %*% ( t(C) %*% (Q * exp.Q) * exp.Ba.betas)
    p2 <- t(Ba) %*% ( (t(C)%*%exp.Q) * exp.Ba.betas)
    wr <- 1:nba + (i-1)*nba
    deriv.ed[wr,i] <- -(p1 + p2)
  }
  out.ed <- list(delta.ed=delta.ed,
                 deriv.ed=deriv.ed)
  return(out.ed)
}


Pr <- 10^-4*diag(nba*n+nF+nF)

ONE <- matrix(0, nb1, nF)
BETAS <- BETAS.SI
PLOT <- TRUE
for(it in 1:10){
  betas <- c(BETAS)
  
  ETA <- MortSmooth_BcoefB(Ba, Bt, BETAS)
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
  tBWBpP <- tBWB + P
  tBWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
  
  fx1 <- tBWBpP%*%betas - c(tBWz)
  fx2 <- tBWBpP
  
  BETASF <- BETAS[,(nbt1+1):nbt]
  
  constr.e0 <- Constr.e0(BETASF)
  hx <- constr.e0$delta.e0
  hx1 <- constr.e0$deriv.e0
  
  constr.ed <- Constr.ed(BETASF)
  gx <- constr.ed$delta.ed
  gx1 <- constr.ed$deriv.ed
  
  Am <- cbind( rbind(ONE, hx1),
               rbind(ONE, gx1))
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(ncol(Am)))) + Pr
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  D <- matrix(d, nba, nbt)
  BETAS <- BETAS + D
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
  if(PLOT){
    par(mfrow=c(1,2))
    e0.C <- apply(exp(ETA), 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
    ranx <- range(t)
    rany <- range(e01, e0.tar, e0.SI)
    plot(t1, e01, xlim=ranx, ylim=rany)
    points(tF, e0.tar, pch=16, t="b", lwd=2)
    lines(t, e0.SI, col=2, lwd=2, lty=2)
    lines(t, e0.C, col=4, lwd=3)
    abline(v=t[n1]+0.5)
    
    ed.C <- apply(exp(ETA), 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
    ranx <- range(t)
    rany <- range(ed1, ed.tar, ed.SI)
    plot(t1, ed1, xlim=ranx, ylim=rany)
    points(tF, ed.tar, pch=16, t="b", lwd=2)
    lines(t, ed.SI, col=2, lwd=2, lty=2)
    lines(t, ed.C, col=4, lwd=3)
    abline(v=t[n1]+0.5)
    par(mfrow=c(1,1))
  }
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)


## plotting
rany <- range(ETA1, ETA.C, finite=TRUE)
matplot(ages, ETA1, col=colnT, t="l", lty=3, ylim=rany,
        xlab="ages", ylab="log-mortality")
whi <- floor(seq(1,n,length=6))
for(i in whi) lines(ages, ETA.C[,i], col=coln[i], lwd=4)
legend("topleft", inset=0.1, legend=years[whi], col=coln[whi], lwd=4)


par(mfrow=c(2,5))
whi <- floor(seq(1,m,length=10))
ranx <- range(years)
for(i in whi){
  rany <- range(ETA1[i,], ETA.C[i,], finite=TRUE)
  plot(years, ETA.C[i,], t="l", lwd=4, col=colmT[i], main=paste(ages[i]),
       ylim=rany, xlim=ranx,xlab="ages", ylab="log-mortality")
  points(years1, ETA1[i,], col=colm[i], pch=16)
  abline(v=years[n1]+0.5)
}
par(mfrow=c(1,1))






hxFUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  e0.F <- apply(MU.F, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
  hx <- e0.F - e0.tar
  return(hx)
}

hx1FUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  exp.sum.MU.F <- -exp(colSums(-MU.F))
  
  hx1 <- matrix(0, nbF, nF)
  for(k in 1:nF){
    for(i in 1:nbF){
      Bi <- BF[1:m+(k-1)*m,i]
      muk <- MU.F[,k]
      etak <- ETA.F[,k]
      p1 <- Bi[1]*exp(muk[2] + etak[1]) + Bi[1]*muk[1] + Bi[2]*muk[2]
      for(j in 3:m){
        p1 <- p1*exp(muk[j]) + sum(Bi[1:j]*muk[1:j])
      }
      hx1[i,k] <- p1*exp.sum.MU.F[k]
    }
  }
  return(hx1)
}


BETAS <- BETAS.S
PLOT <- TRUE
for(it in 1:50){
  betas <- c(BETAS)
  
  ETA <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBx, RTBt, nbx, nbt, WW)
  tBWBpP <- tBWB + P
  tBWz <- MortSmooth_BcoefB(t(Bx), t(Bt), (WW * z))
  
  fx1 <- tBWBpP%*%betas - c(tBWz)
  fx2 <- tBWBpP
  
  BETASF <- BETAS[,(nbt1+1):nbt]
  hx <- hxFUN(BETASF)
  hx1 <- hx1FUN(BETASF)
  
  Am <- rbind(ONE, hx1)
  bv <- -c(hx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(ncol(Am))))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  D <- matrix(d, nbx, nbt)
  BETAS <- BETAS + D
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
  if(PLOT){
    e0.C <- apply(exp(ETA), 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
    ranx <- range(t)
    rany <- range(e01, e0.tar, e0.S)
    plot(t1, e01, xlim=ranx, ylim=rany)
    lines(t, e0.S, col=2, lwd=2, lty=2)
    lines(t, e0.C, col=4, lwd=3)
    points(tF, e0.tar, pch=16, t="b", lwd=2)
    abline(v=t[n1]+0.5)
  }
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)


# ranx <- range(t)
# for(i in 1:m){
#   rany <- range(lMX[i,], ETA.S[i,], ETA.C[i,], na.rm = TRUE, finite=TRUE)
#   plot(t1, lMX1[i,], xlim=ranx, ylim=rany)
#   lines(t, ETA.S[i,], col=2, lwd=2)
#   lines(t, ETA.C[i,], col=4, lwd=3)
#   title(main=x[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }
# for(i in 1:n){
#   rany <- range(lMX[,i], ETA.S[,i], ETA.C, na.rm=TRUE, finite=TRUE)
#   plot(x, lMX[,i], ylim=rany)
#   lines(x, ETA.S[,i], col=2, lwd=2)
#   lines(x, ETA.C[,i], col=4, lwd=3)
#   title(main=years[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }



whi <- floor(seq(1,n,length=9))
rany <- range(lMX, ETA.S, ETA.C, 0, na.rm=TRUE, finite=TRUE)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(1,1,t="n", xlim=range(x), ylim=rany, 
       main=paste(years[whi[i]], round(e0.C[whi[i]],2)),
       xlab="age", ylab="log-mortality")
  points(x, lMX[,whi[i]], pch=16, ylim=rany)
  lines(x, ETA.S[,whi[i]], col=2, lwd=2, lty=2)
  lines(x, ETA.C[,whi[i]], col=4, lwd=3)
}
par(mfrow=c(1,1))



plot3d(xx, tt, lmxINF, type="p", col="red",
       xlab="age", ylab="year", zlab="log-mortality", 
       site=5, lwd=15)
surface3d(ages, years, ETA.S, back='line', front='line',
          col=2, lwd=1, alpha = 0.5)
surface3d(ages, years, ETA.C, back='line', front='line',
          col=4, lwd=1, alpha = 0.5)

## compute constrained e0
e0.C <- apply(MU.C, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(years)
rany <- range(e01, e0.tar, e0.S)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
lines(years, e0.C, col=4, lwd=3)
points(yearsF, e0.tar, pch=16, t="b", lwd=2)
abline(v=years[n1]+0.5)



round(cbind(e0.S, e0.C, e0.S-e0.C),2)




## save.image("JAPfem.RData")


## Final Figures for the IWSM abstract
rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
library(MortalitySmooth)
library(HMDdata)
library(vars)
library(forecast)
load("JAPfem.RData")

ranx <- range(years)
rany <- range(e01, e0.tar, e0.S)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
lines(years, e0.C, col=4, lwd=3)
points(yearsF, e0.tar, pch=3, t="b", lwd=2)
abline(v=years[n1]+0.5)












































## LC true parameters, only for setting up the simulation
AlphaTrue <- read.table("AlphaTrue.txt", header=TRUE)
BetaTrue <- read.table("BetaTrue.txt", header=TRUE)
KappaTrue <- read.table("KappaTrue.txt", header=TRUE)
aT <- AlphaTrue[,2]
bT <- BetaTrue[,2]
kT <- KappaTrue[,2]

## dimensions
x <- AlphaTrue[,1]
t <- KappaTrue[,1]
m <- length(x)
n <- length(t)

## exposures
E <- as.matrix(read.table("Exposures.txt", header=TRUE))

## simulating data
One <- matrix(1, nrow=n, ncol=1)  
ETA.T <- aT %*% t(One) + bT %*% t(kT)
MU.T <- exp(ETA.T)
Y.T <- MU.T * E
y <- rpois(m*n, lambda=c(Y.T))
Y <- matrix(y, m, n)
MX <- Y/E
lMX <- log(MX)
mx <- c(MX)
lmx <- c(lMX)

## plotting simulated and true log-mortality
lmxNA <- lmx
lmxNA[is.infinite(lmx)] <- NA
xx <- rep(x, n)
tt <- rep(t, each=m)
plot3d(xx, tt, lmxNA, type="p", col="red", 
       xlab="age", ylab="year", zlab="log-mortality", site=5, lwd=15)
surface3d(x, t, ETA.T, back='line', front='line', 
          col=5, lwd=1, alpha = 0.5)


## assuming we know only asfr
## up to max(t1) (e.g. 1990), but we have knowledge 
## of e0 up to 2006 (actually computed)
## and we wanna forecast all asfr up to
## 2006 by constraining each "future" age-pattern to 
## follow the previously mentioned summary measure (e0)

t1 <- t[1]:1980
n1 <- length(t1)
tF <- t[-c(1:n1)]
nF <- length(tF)

Y1 <- Y[,1:n1]
E1 <- E[,1:n1]
MX1 <- Y1/E1
lMX1 <- log(MX1)
e01 <- apply(MX1, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)

## keep aside all asfr and summary measures 
## we in theory know
Yobs <- Y
Eobs <- E
MXobs <- MX
lMXobs <- lMX
## e0 only for future years, target e0
## taken from the underlying force of mortality
e0.tar <- apply(MU.T[,1:nF+n1], 2, function(mx) sum(exp(-cumsum(mx)))+0.5)


## replace with arbitrary values all info
## after max(t1), what the last year we assume to know
## this step is need to estimate up to max(t1)
## and simoultaneously forecast afterwards
Y <- matrix(10, m, n)
Y[1:m,1:n1] <- Y1
E <- matrix(10, m, n)
E[1:m,1:n1] <- E1
MX <- Y/E
MX[,(n1+1):n] <- NA

## set weights equal to zero
## to years max(t1):2016
## 0/1 weight for forecasts
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
# ## smoothing parameters
# v1 <- c(rep(1,n1-1), rep(0,nF-1))
# V1 <- diag(v1)
# tDDt1 <- t(Dt)%*%V1%*%Dt
# Pt1 <- kronecker(tDDt1, diag(nbx))
# vF <- c(rep(0,n1-1), rep(1,nF-1))
# VF <- diag(vF)
# tDDtF <- t(Dt)%*%VF%*%Dt
# PtF <- kronecker(tDDtF, diag(nbx))
# lambda.tF <- 10^2
# PF <- lambda.x*Px + lambda.t*Pt1 + lambda.tF*PtF


## data in vector
y <- c(Y)
e <- c(E)
wei <- c(WEI)
mx <- c(MX)
lmx <- log(mx)

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


# ranx <- range(t)
# for(i in 1:m){
#   rany <- range(lMX[i,], ETA.S[i,], na.rm = TRUE, finite=TRUE)
#   plot(t1, lMX1[i,], xlim=ranx, ylim=rany)
#   lines(t, ETA.S[i,], col=2, lwd=2)
#   title(main=x[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }
# for(i in 1:n){
#   rany <- range(MX[,i], MU.S[,i], na.rm=TRUE)
#   plot(x, MX[,i], ylim=rany)
#   lines(x, MU.S[,i], col=2, lwd=2)
#   title(main=t[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }


plot3d(xx, tt, lmx, type="p", col="red",
       xlab="age", ylab="year", zlab="log-mortality", 
       site=5, lwd=15)
surface3d(x, t, ETA.T, back='line', front='line',
          col=5, lwd=1, alpha = 0.5)
surface3d(x, t, ETA.S, back='line', front='line',
          col=2, lwd=1, alpha = 0.5)

## compute estimated e0
e0.S <- apply(MU.S, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(t)
rany <- range(e01, e0.tar, e0.S)
plot(t1, e01, xlim=ranx, ylim=rany)
lines(t, e0.S, col=2, lwd=2, lty=2)
points(tF, e0.tar, pch=16, t="b", lwd=2)
abline(v=t[n1]+0.5)


ONE <- matrix(0, nb1, nF)
BETAS <- BETAS.S
BETASF <- BETAS[,(nbt1+1):nbt]
## constraints on e0 over "future" years
hxFUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  e0.F <- apply(MU.F, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
  hx <- e0.F - e0.tar
  return(hx)
}

#exp(-exp(b11 a1 + b12 a2 + b13 a3)) + exp(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3))+ exp(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3)-exp(b31 a1 + b32 a2 + b33 a3))+ exp(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3)-exp(b31 a1 + b32 a2 + b33 a3)-exp(b41 a1 + b42 a2 + b43 a3)) + 0.5
hx1FUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  exp.sum.MU.F <- -exp(colSums(-MU.F))
  
  hx1 <- matrix(0, nbF, nF)
  for(k in 1:nF){
    for(i in 1:nbF){
      Bi <- BF[1:m+(k-1)*m,i]
      muk <- MU.F[,k]
      etak <- ETA.F[,k]
      p1 <- Bi[1]*exp(muk[2] + etak[1]) + Bi[1]*muk[1] + Bi[2]*muk[2]
      for(j in 3:m){
        p1 <- p1*exp(muk[j]) + sum(Bi[1:j]*muk[1:j])
      }
      hx1[i,k] <- p1*exp.sum.MU.F[k]
    }
  }
  return(hx1)
}


BETAS <- BETAS.S
PLOT <- TRUE
for(it in 1:50){
  betas <- c(BETAS)
  
  ETA <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU <- E*exp(ETA)
  W <- c(WEI*MU)
  z <- ETA + (1/MU) * (Y - MU)
  z[which(WEI == 0)] <- 0
  WW <- WEI * W
  tBWB <- MortSmooth_BWB(RTBx, RTBt, nbx, nbt, WW)
  tBWBpP <- tBWB + P
  tBWz <- MortSmooth_BcoefB(t(Bx), t(Bt), (WW * z))
  
  fx1 <- tBWBpP%*%betas - c(tBWz)
  fx2 <- tBWBpP
  
  BETASF <- BETAS[,(nbt1+1):nbt]
  hx <- hxFUN(BETASF)
  hx1 <- hx1FUN(BETASF)
  
  Am <- rbind(ONE, hx1)
  bv <- -c(hx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(ncol(Am))))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  D <- matrix(d, nbx, nbt)
  BETAS <- BETAS + D
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
  if(PLOT){
    e0.C <- apply(exp(ETA), 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
    ranx <- range(t)
    rany <- range(e01, e0.tar, e0.S)
    plot(t1, e01, xlim=ranx, ylim=rany)
    lines(t, e0.S, col=2, lwd=2, lty=2)
    lines(t, e0.C, col=4, lwd=3)
    points(tF, e0.tar, pch=16, t="b", lwd=2)
    abline(v=t[n1]+0.5)
  }
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)


ranx <- range(t)
for(i in 1:m){
  rany <- range(lMX[i,], ETA.S[i,], ETA.C[i,], na.rm = TRUE, finite=TRUE)
  plot(t1, lMX1[i,], xlim=ranx, ylim=rany)
  lines(t, ETA.S[i,], col=2, lwd=2)
  lines(t, ETA.C[i,], col=4, lwd=3)
  lines(t, ETA.T[i,], col=3, lwd=2, lty=3)
  title(main=x[i])
  # locator(1)
  Sys.sleep(0.1)
}
# for(i in 1:n){
#   rany <- range(MX[,i], MU.S[,i], na.rm=TRUE)
#   plot(x, MX[,i], ylim=rany)
#   lines(x, MU.S[,i], col=2, lwd=2)
#   title(main=t[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }


plot3d(xx, tt, lmx, type="p", col="red",
       xlab="age", ylab="year", zlab="log-mortality", 
       site=5, lwd=15)
# surface3d(x, t, ETA.T, back='line', front='line',
#           col=5, lwd=1, alpha = 0.5)
surface3d(x, t, ETA.C, back='line', front='line',
          col=2, lwd=1, alpha = 0.5)

## compute constrained e0
e0.C <- apply(MU.C, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(t)
rany <- range(e01, e0.tar, e0.S)
plot(t1, e01, xlim=ranx, ylim=rany)
lines(t, e0.S, col=2, lwd=2, lty=2)
lines(t, e0.C, col=4, lwd=3)
points(tF, e0.tar, pch=16, t="b", lwd=2)
abline(v=t[n1]+0.5)

























## mortality 1D example ## mortality 1D example ## mortality 1D example 
## mortality 1D example ## mortality 1D example ## mortality 1D example 
## mortality 1D example ## mortality 1D example ## mortality 1D example 
rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(MortalitySmooth)
library(rgl)

e <- c(32649,32508,32376,32252,32133,32018,31905,31792,31678,31560,31438,31343,31302,31300,31327,31370,31416,31454,31470,31452,31388,31267,31128,31012,30910,30810,30702,30575,30419,30224,29981,29680,29313,28902,28474,28027,27559,27067,26550,26008,25439,24843,24220,23569,22933,22340,21768,21199,20617,20005,19352,18648,17883,17055,16160,15303,14552,13861,13191,12511,11797,11031,10203,9311,8364,7376,6408,5512,4692,3950,3288,2705,2198,1763,1396,1090,839,637,476,351,254,181,127,87,59,39,26)
eta.T <- c(-7.903,-7.716,-7.546,-7.396,-7.263,-7.148,-7.053,-6.978,-6.916,-6.881,-6.855,-6.82,-6.788,-6.75,-6.705,-6.657,-6.608,-6.555,-6.501,-6.449,-6.398,-6.344,-6.287,-6.228,-6.168,-6.105,-6.042,-5.977,-5.91,-5.842,-5.773,-5.702,-5.628,-5.552,-5.473,-5.394,-5.317,-5.24,-5.162,-5.085,-5.008,-4.93,-4.851,-4.77,-4.687,-4.602,-4.514,-4.423,-4.329,-4.235,-4.139,-4.041,-3.942,-3.843,-3.742,-3.638,-3.532,-3.425,-3.317,-3.208,-3.099,-2.991,-2.883,-2.775,-2.667,-2.561,-2.455,-2.349,-2.244,-2.141,-2.039,-1.939,-1.84,-1.744,-1.649,-1.554,-1.462,-1.372,-1.284,-1.201,-1.122,-1.047,-0.975,-0.907,-0.843,-0.782,-0.724)
m <- length(eta.T)
x <- 12:98
mu.T <- exp(eta.T)
y.T <- e*mu.T
y <- rpois(m, y.T)
mx <- y/e
lmx <- log(mx)

## simulated e0
e0.obs <- sum(exp(-cumsum(mx)))+0.5

plot(x, lmx)
lines(x, eta.T, col=2)

## B-splines
B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                      ndx=floor(m/4), deg=3)
nb <- ncol(B)
## penalty stuff
lambda <- 10^1
D <- diff(diag(nb), diff=2)
tDD <- t(D)%*%D
P <- lambda * tDD

eta <- log((y+1)/(e+1))
for(it in 1:20){
  mu <- e*exp(eta)
  z <- (y-mu)/mu + eta
  w <- as.vector(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  betas <- solve(tBWBpP, tBWz)
  old.eta <- eta
  eta <- B%*%betas
  dif.eta <- max(abs(old.eta - eta))
  cat(it, dif.eta, "\n")
  if(dif.eta < 1e-4) break
}
eta.S <- eta
mu.S <- exp(eta.S)
betas.S <- betas

plot(x, lmx)
lines(x, eta.T, col=2)
lines(x, eta.S, col=3, lwd=2, lty=2)

## estimated e0
e0.S <- sum(exp(-cumsum(mu.S)))+0.5

## target e0
e0.tar <- 65




hxFUN <-  function(betas){
  eta.F <- B%*%betas
  mu.F <- exp(eta.F)
  e0.F <- sum(exp(-cumsum(mu.F)))+0.5
  hx <- e0.F - e0.tar
  return(hx)
}

hx1FUN <- function(betas){
  eta.F <- B%*%betas
  mu.F <- exp(eta.F)
  sum.exp.mu <- -exp(sum(-mu.F))
  hx1 <- matrix(0, nb, 1)
  for(i in 1:nb){
    Bi <- B[,i]
    p1 <- Bi[1]*exp(mu.F[2] + eta.F[1]) + Bi[1]*mu.F[1] + Bi[2]*mu.F[2]
    for(j in 3:m){
       p1 <- p1*exp(mu.F[j]) + sum(Bi[1:j]*mu.F[1:j])
    }
    hx1[i,1] <- p1*sum.exp.mu
  }
  return(hx1)
}

betas <- betas.S
for(it in 1:50){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y-mu)/mu + eta
  w <- as.vector(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - c(tBWz)
  fx2 <- tBWBpP
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  Am <- hx1
  bv <- -c(hx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(ncol(Am))))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
}  

eta.C <- eta
mu.C <- exp(eta.C)
betas.C <- betas

## constrained e0
e0.C <- sum(exp(-cumsum(mu.C)))+0.5

cbind(e0.tar, e0.C, e0.S)

## plotting
plot(x, lmx)
lines(x, eta.T, col=2)
lines(x, eta.S, col=3, lwd=2, lty=2)
lines(x, eta.C, col=4, lwd=3)

