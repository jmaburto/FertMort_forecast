rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
library(MortalitySmooth)
library(HMDdata)
library(vars)
library(forecast)

cou <- "Italy"
sex <- "F"
ages <- 0:105
m <- length(ages)
years1 <- 1960:2014
n1 <- length(years1)
Y1 <- matrix(selectHMDdata(cou, "Deaths", sex, ages, years1), m, n1)
E1 <- matrix(selectHMDdata(cou, "Exposures", sex, ages, years1), m, n1)
MX1 <- Y1/E1
lMX1 <- log(MX1)
mx1 <- c(MX1)
lmx1 <- c(lMX1)
## domains for the estimation
x <- 1:m-1
t1 <- 1:n1-1

## plotting in 3D
lmx1NA <- lmx1
lmx1NA[is.infinite(lmx1)] <- NA
xx1 <- rep(ages, n1)
tt1 <- rep(years1, each=m)
plot3d(xx1, tt1, lmx1NA, type="p", col="red", 
       xlab="age", ylab="year", zlab="log-mortality", site=5, lwd=15)

years <- min(years1):2050
n <- length(years)
t <- 1:n-1
tF <- t[-c(1:n1)]
nF <- length(tF)
yearsF <- years[-c(1:n1)]
MX1na <- MX1
MX1na[is.nan(MX1)] <- 0
## observed e0
e01 <- apply(MX1na, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
## observed e-dagger
ed1 <-  apply(MX1na, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
## observed Gini
G1  <-  apply(MX1na, 2, function(mx) 1 - sum(exp(-cumsum(mx))^2)/(sum(exp(-cumsum(mx))))-.0006)

# MX1vec <- c(MX1)
# In <- diag(n1)
# Im <- diag(m)
# C <- lower.tri(Im, diag = TRUE)
# C[C==1] <- -1
# CC <- kronecker(In, C)
# ones.m <- matrix(1,1,m)
# Iones <- kronecker(In, ones.m)
# e01 <-  Iones %*%   exp(CC%*%MX1vec) + 0.5
# ed1 <- - Iones %*% (exp(CC%*%MX1vec) * CC%*%MX1vec)
# 
# mx1 <- MX1[,1]
# 
# ones.m %*% exp(C%*%mx1) + 0.5
# e01[1]
# -t(exp(C%*%mx1)) %*% (C%*%mx1)
# ed1[1]

## matrices/vectors arbitrary values
Y <- matrix(10, m, n)
Y[1:m,1:n1] <- Y1
E <- matrix(10, m, n)
E[1:m,1:n1] <- E1
MX <- Y/E
MX[,(n1+1):n] <- NA
lMX <- log(MX)
mx <- c(MX)
lmx <- c(lMX)

## set weights equal to zero for future years
WEI <- matrix(0, m, n)
WEI[,1:n1] <- 1
## set equal to zero when we have no exposures
WEI[E1==0] <- 0

## using MortalitySmooth optimized lambdas
fitPS <- Mort2Dsmooth(x=x, y=t, Z=Y, offset=log(E), W=WEI,
                      overdispersion = TRUE)
summary(fitPS)
plot(fitPS)

ETA.PS <- fitPS$logmortality
MU.PS <- exp(ETA.PS)
e0.PS <- apply(MU.PS, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(years)
rany <- range(e01, e0.PS)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.PS, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
## compute estimated e-dagger
ed.PS <-  apply(MU.PS, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
ranx <- range(years)
rany <- range(ed1, ed.PS)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)

## B-splines

## over ages
## w/o age 0
x0 <- x[-1]
m0 <- m-1
Bx0 <- MortSmooth_bbase(x=x0, xl=min(x0), xr=max(x0),
                        ndx=floor(m0/4), deg=3)
nbx0 <- ncol(Bx0)
Bx <- cbind(0, Bx0)
Bx <- rbind(c(1, rep(0,nbx0)), Bx)
# Bx <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
#                        ndx=floor(m/4), deg=3)
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


## overall basis dimensions
nb <- nbx*n
nb1 <- nbx*n1
nbF <- nbtF*nbx

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
Dx <- diff(diag(nbx), diff=2)
tDDx <- t(Dx)%*%Dx
Dt <- diff(diag(nbt), diff=2)
tDDt <- t(Dt)%*%Dt
Px <- kronecker(diag(nbt), tDDx)
Pt <- kronecker(tDDt, diag(nbx))
## smoothing parameters
lambdas <- c(10^0, 10^5)
lambda.x <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.x * Px + lambda.t * Pt

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

## useful objects
ones.m <- matrix(1,m,1)
Im <- diag(m)
C <- lower.tri(Im, diag = TRUE)
C[C==1] <- -1

## e0
deriv.e0.S <- matrix(0, m, n)
for(i in 1:n){
  Ba <- ETA.S[,i]
  exp.Ba <- exp(Ba)
  Q <- C %*% exp.Ba
  exp.Q <- exp(Q)
  d.exp.Q <- diag(c(exp.Q))
  d.exp.Ba <- diag(c(exp.Ba)) 
  deriv.coef <- t(ones.m) %*% d.exp.Q %*% C %*% d.exp.Ba %*% Bx
  deriv.e0.S[,i] <- Bx%*%c(deriv.coef)
}
matplot(ages, deriv.e0.S, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
## e-dagger
deriv.ed.S <- matrix(0, m, n)
for(i in 1:n){
  Ba <- ETA.S[,i]
  exp.Ba <- exp(Ba)
  Q <- C %*% exp.Ba
  exp.Q <- exp(Q)
  p1 <- t(Bx) %*% ( t(C) %*% (Q * exp.Q) * exp.Ba)
  p2 <- t(Bx) %*% ( (t(C)%*%exp.Q) * exp.Ba)
  deriv.coef <- -(p1 + p2)
  deriv.ed.S[,i] <- Bx%*%c(deriv.coef)
}
matplot(ages, deriv.ed.S, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)




# ranx <- range(years)
# for(i in 1:m){
#   rany <- range(lMX[i,], ETA.S[i,], na.rm = TRUE, finite=TRUE)
#   plot(years1, lMX1[i,], xlim=ranx, ylim=rany)
#   lines(years, ETA.S[i,], col=2, lwd=2)
#   title(main=ages[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }
# for(i in 1:n){
#   rany <- range(lMX[,i], ETA.S[,i], na.rm=TRUE, finite=TRUE)
#   plot(x, lMX[,i], ylim=rany)
#   lines(x, ETA.S[,i], col=2, lwd=2)
#   title(main=years[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }


xx <- rep(ages, n)
tt <- rep(years, each=m)
lmxINF <- lmx
lmxINF[is.infinite(lmx)] <- NA
plot3d(xx, tt, lmxINF, type="p", col="red",
       xlab="age", ylab="year", zlab="log-mortality", 
       site=5, lwd=15)
surface3d(ages, years, ETA.S, back='line', front='line',
          col=2, lwd=1, alpha = 0.5)

## compute estimated e0
e0.S <- apply(MU.S, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(years)
rany <- range(e01, e0.S)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
## compute estimated e-dagger
ed.S <-  apply(MU.S, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
ranx <- range(years)
rany <- range(ed1, ed.S)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)


## forecasting e0
# e0.ts <- ts(e01, start = t1[1])
# ds <- 0
# df.test.e0 <- ur.df(e0.ts,lags=2,'trend')
# if (df.test.e0@teststat[1] > df.test.e0@cval[1,][1]) ds <- 1 ## if TRUE -> d = 1
# e0.mod <- auto.arima(e0.ts,d=ds,max.p=3,max.q=3,trace=TRUE)
# pred.e0 <- forecast(e0.mod, h=nF,level=95)
# 
# ranx <- range(years)
# rany <- range(e01, e0.S, pred.e0$lower, pred.e0$upper)
# plot(years1, e01, xlim=ranx, ylim=rany)
# lines(years, e0.S, col=2, lwd=2, lty=2)
# ## future mean for e0
# e0mu <- c(e01[n1], as.vector(pred.e0$mean))
# ## future standard deviation for e0
# e0sig <- c(0.001, as.vector((pred.e0$upper- pred.e0$lower)/4))
# ## product of each combination of elements of the sd vector
# e0var <- e0sig %*% t(e0sig)
# ## variance-covariance matrix
# e0Cor <- matrix(0.9999, nF+1, nF+1)
# diag(e0Cor) <- 1
# ## scalar multiply e0var by the variance-covariance matrix 
# ## obtaining covariance matrix of the variables
# e0Sig <- e0var * e0Cor
# ## generate target e0
# e0.tar <- mvrnorm(n = 1, mu=e0mu, Sigma=e0Sig)
# e0.tar <- e0.tar[-1]
# lines(yearsF, e0.tar, col=4, lwd=2, lty=2)

## from UN WPP 
yrsUN1 <- seq(2017.5, 2097.5, 5)
## general UN e0 data
UNfile <- paste("e0", sex, ".csv", sep="")
eUNall <- read.table(UNfile, sep=",", header=TRUE)
eUN <- t(as.matrix(eUNall[which(eUNall[,1]==cou), 5:21]))
eUN <- eUN[,3]
## only for Japan, females
# ## age 0
# eUN <- c(87.18,87.89,88.56,89.24,89.89,90.54,91.17,91.80,92.42,93.02,93.63,94.20,94.80,95.40,95.98,96.56,97.14)
# ## age 1
# eUN <- c(86.34,87.03,87.68,88.35,88.98,89.61,90.24,90.87,91.47,92.07,92.68,93.25,93.84,94.43,95.02,95.58,96.16)
# ## age 5
# eUN <- c(82.40,83.08,83.72,84.39,85.01,85.64,86.27,86.89,87.49,88.08,88.69,89.26,89.85,90.44,91.03,91.59,92.17)
# ## age 15
# eUN <- c(72.45,73.13,73.77,74.43,75.05,75.68,76.29,76.92,77.52,78.11,78.71,79.28,79.87,80.45,81.04,81.60,82.18)
# ## age 10
# eUN <- c(77.43,78.11,78.75,79.41,80.03,80.66,81.28,81.90,82.51,83.10,83.70,84.27,84.86,85.45,86.04,86.60,87.17)
fitUN <- lm(eUN~yrsUN1)
e0.tar <- predict(fitUN, list(yrsUN1=yearsF))
ranx <- range(years)
rany <- range(e01, e0.S)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
points(yearsF, e0.tar, col=5)







## rough forecast of e-dagger
ed.ts <- ts(ed1, start = t1[1])
ed.mod <- auto.arima(ed.ts, d=2, max.p=2, max.q=2)
pred.ed <- forecast(ed.mod, h=nF,level=95)
plot(pred.ed)
ranx <- range(years)
rany <- range(ed1, ed.S, pred.ed$lower, pred.ed$upper)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)

# ## future mean for e-dagger
# edmu <- c(ed1[n1], as.vector(pred.ed$mean))
# ## future standard deviation for e-dagger
# edsig <- c(0.001, as.vector((pred.ed$upper- pred.ed$lower)/4))
# ## product of each combination of elements of the sd vector
# edvar <- edsig %*% t(edsig)
# ## variance-covariance matrix
# edCor <- matrix(0.9999, nF+1, nF+1)
# diag(edCor) <- 1
# ## scalar multiply e0var by the variance-covariance matrix 
# ## obtaining covariance matrix of the variables
# edSig <- edvar * edCor
# is.positive.definite(edSig)
# ## generate target e0
# ed.tar <- mvrnorm(n = 1, mu=edmu, Sigma=edSig)
# ed.tar <- ed.tar[-1]
ed.tar <- as.vector(pred.ed$mean)
lines(yearsF, ed.tar, col=4, lwd=2, lty=2)

par(mfrow=c(1,2))
## plotting
ranx <- range(t)
rany <- range(e01, e0.tar, e0.S)
plot(t1, e01, xlim=ranx, ylim=rany)
lines(t, e0.S, col=2, lwd=2, lty=2)
points(tF, e0.tar, pch=16, t="b", lwd=2)
abline(v=t[n1]+0.5)

ranx <- range(t)
rany <- range(ed1, ed.tar, ed.S)
plot(t1, ed1, xlim=ranx, ylim=rany)
lines(t, ed.S, col=2, lwd=2, lty=2)
points(tF, ed.tar, pch=16, t="b", lwd=2)
abline(v=t[n1]+0.5)
par(mfrow=c(1,1))






## constraints on e0 over "future" years
Constr.e0 <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  e0.F <- apply(MU.F, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
  delta.e0 <- e0.F - e0.tar
  deriv.e0 <- matrix(0, nbF, nF)
  for(i in 1:nF){
    Ba <- ETA.F[,i]
    exp.Ba <- exp(Ba)
    Q <- C %*% exp.Ba
    exp.Q <- exp(Q)
    d.exp.Q <- diag(c(exp.Q))
    d.exp.Ba <- diag(c(exp.Ba)) 
    wr <- 1:nbx + (i-1)*nbx
    deriv.e0[wr,i] <- t(ones.m) %*% d.exp.Q %*% C %*% d.exp.Ba %*% Bx
  }
  out.e0 <- list(delta.e0=delta.e0,
                 deriv.e0=deriv.e0)
  return(out.e0)
}

Constr.ed <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  ed.F <- apply(MU.F, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
  delta.ed <- ed.F - ed.tar
  deriv.ed <- matrix(0, nbF, nF)
  for(i in 1:nF){
    Ba <- ETA.F[,i]
    exp.Ba <- exp(Ba)
    Q <- C %*% exp.Ba
    exp.Q <- exp(Q)
    p1 <- t(Bx) %*% ( t(C) %*% (Q * exp.Q) * exp.Ba)
    p2 <- t(Bx) %*% ( (t(C)%*%exp.Q) * exp.Ba)
    wr <- 1:nbx + (i-1)*nbx
    deriv.ed[wr,i] <- -(p1 + p2)
  }
  out.ed <- list(delta.ed=delta.ed,
                 deriv.ed=deriv.ed)
  return(out.ed)
}

ONE <- matrix(0, nb1, nF)
BETAS <- BETAS.S
PLOT <- TRUE
for(it in 1:10){
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
               cbind(t(Am), 0*diag(ncol(Am))))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  D <- matrix(d, nbx, nbt)
  BETAS <- BETAS + D
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
  if(PLOT){
    par(mfrow=c(1,2))
    e0.C <- apply(exp(ETA), 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
    ranx <- range(t)
    rany <- range(e01, e0.tar, e0.S)
    plot(t1, e01, xlim=ranx, ylim=rany)
    points(tF, e0.tar, pch=16, t="b", lwd=2)
    lines(t, e0.S, col=2, lwd=2, lty=2)
    lines(t, e0.C, col=4, lwd=3)
    abline(v=t[n1]+0.5)
    
    ed.C <- apply(exp(ETA), 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
    ranx <- range(t)
    rany <- range(ed1, ed.tar, ed.S)
    plot(t1, ed1, xlim=ranx, ylim=rany)
    points(tF, ed.tar, pch=16, t="b", lwd=2)
    lines(t, ed.S, col=2, lwd=2, lty=2)
    lines(t, ed.C, col=4, lwd=3)
    abline(v=t[n1]+0.5)
    par(mfrow=c(1,1))
  }
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)

S.S <- apply(MU.S, 2, function(mx) exp(-cumsum(mx)))
F.S <- S.S*MU.S
S.C <- apply(MU.C, 2, function(mx) exp(-cumsum(mx)))
F.C <- S.C*MU.C


par(mfrow=c(1,2))
matplot(F.S, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
matplot(F.C, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
par(mfrow=c(1,1))
## e0
deriv.e0.C <- matrix(0, m, n)
for(i in 1:n){
  Ba <- ETA.C[,i]
  exp.Ba <- exp(Ba)
  Q <- C %*% exp.Ba
  exp.Q <- exp(Q)
  d.exp.Q <- diag(c(exp.Q))
  d.exp.Ba <- diag(c(exp.Ba)) 
  deriv.coef <- t(ones.m) %*% d.exp.Q %*% C %*% d.exp.Ba %*% Bx
  deriv.e0.C[,i] <- Bx %*%c(deriv.coef)
}
matplot(ages, deriv.e0.C, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)

## e-dagger
deriv.ed.C <- matrix(0, m, n)
for(i in 1:n){
  Ba <- ETA.C[,i]
  exp.Ba <- exp(Ba)
  Q <- C %*% exp.Ba
  exp.Q <- exp(Q)
  p1 <- t(Bx) %*% ( t(C) %*% (Q * exp.Q) * exp.Ba)
  p2 <- t(Bx) %*% ( (t(C)%*%exp.Q) * exp.Ba)
  deriv.coef <- -(p1 + p2)
  deriv.ed.C[,i] <- Bx%*%c(deriv.coef)
}
matplot(ages, deriv.ed.C, t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)


par(mfrow=c(2,2))
matplot(ages, deriv.e0.S[,1:nF+n1], t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
matplot(ages, deriv.ed.S[,1:nF+n1], t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
matplot(ages, deriv.e0.C[,1:nF+n1], t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
matplot(ages, deriv.ed.C[,1:nF+n1], t="l", col=rainbow(n), lty=c(rep(1,n1), rep(3,nF)))
abline(h=0, col=8, lwd=2, lty=2)
par(mfrow=c(1,1))





ranx <- range(t)
for(i in 1:m){
  rany <- range(lMX[i,], ETA.S[i,], ETA.C[i,], na.rm = TRUE, finite=TRUE)
  plot(t1, lMX1[i,], xlim=ranx, ylim=rany)
  lines(t, ETA.S[i,], col=2, lwd=2)
  lines(t, ETA.C[i,], col=4, lwd=3)
  title(main=x[i])
  #locator(1)
  Sys.sleep(0.1)
}


rany <- range(lMX, ETA.S, ETA.C, na.rm=TRUE, finite=TRUE)
for(i in 1:n){
  matplot(x, ETA.C, col=8, t="l", lty=1)
  points(x, lMX[,i], pch=16, cex=0.8)
  lines(x, ETA.S[,i], col=2, lwd=2)
  lines(x, ETA.C[,i], col=4, lwd=3)
  title(main=years[i])
  # locator(1)
  Sys.sleep(0.1)
}



whi <- floor(seq(1,m,length=9))
par(mfrow=c(3,3))
for(i in 1:9){
  rany <- range(lMX[whi[i],], ETA.S[whi[i],], ETA.C[whi[i],], na.rm=TRUE, finite=TRUE)
  plot(1,1,t="n", xlim=range(years), ylim=rany, 
       main=paste(ages[whi[i]]),
       xlab="years", ylab="log-mortality")
  points(years, lMX[whi[i],], pch=16, ylim=rany)
  lines(years, ETA.S[whi[i],], col=2, lwd=2, lty=2)
  lines(years, ETA.PS[whi[i],], col=5, lwd=2, lty=3)
  lines(years, ETA.C[whi[i],], col=4, lwd=3)
}
par(mfrow=c(1,1))



whi <- floor(seq(1,n,length=9))
rany <- range(lMX, ETA.S, ETA.C, 0, na.rm=TRUE, finite=TRUE)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(1,1,t="n", xlim=range(x), ylim=rany, 
       main=paste(years[whi[i]], round(e0.C[whi[i]],2)),
       xlab="age", ylab="log-mortality")
  points(x, lMX[,whi[i]], pch=16, ylim=rany)
  lines(x, ETA.S[,whi[i]], col=2, lwd=2, lty=2)
  lines(x, ETA.PS[,whi[i]], col=5, lwd=2, lty=3)
  lines(x, ETA.C[,whi[i]], col=4, lwd=3)
}
par(mfrow=c(1,1))



plot3d(xx, tt, lmxINF, type="p", col="red",
       xlab="age", ylab="year", zlab="log-mortality", 
       site=5, lwd=15)
# surface3d(ages, years, ETA.S, back='line', front='line',
#           col=2, lwd=1, alpha = 0.5)
surface3d(ages, years, ETA.C, back='line', front='line',
          col=4, lwd=1, alpha = 0.5)
surface3d(ages, years, ETA.PS, back='line', front='line',
          col=5, lwd=1, alpha = 0.5)

par(mfrow=c(1,2))
## compute constrained e0
e0.C <- apply(MU.C, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(years)
rany <- range(e01, e0.tar, e0.PS)
plot(years1, e01, xlim=ranx, ylim=rany)
points(yearsF, e0.tar, pch=16, t="b", lwd=2)
lines(years, e0.C, col=4, lwd=3)
lines(years, e0.PS, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
## compute constrained e-dagger
ed.C <- apply(MU.C, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
ranx <- range(years)
rany <- range(ed1, ed.tar, ed.PS)
plot(years1, ed1, xlim=ranx, ylim=rany)
points(yearsF, ed.tar, pch=16, t="b", lwd=2)
lines(years, ed.C, col=4, lwd=3)
lines(years, ed.PS, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
par(mfrow=c(1,1))


round(cbind(years, e0.S, e0.C, e0.S-e0.C),2)
round(cbind(years, ed.S, ed.C, ed.S-ed.C),2)


save.image("ITAfemMort.RData")




































## e-dagger from UN life tables, only for japan, femelas
JPNmxUN <- read.table("JPNmxUN.csv", sep=",", header=TRUE)
len <- JPNmxUN$length
xUN <- JPNmxUN$age
MXun <- JPNmxUN[,-c(1,2)]
MXLENun <- matrix(0, sum(len), ncol(MXun))
for(i in 1:ncol(MXLENun)){
  MXLENun[,i] <- rep(MXun[,i], len)
}
e0.UN <- apply(MXLENun, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ed.UN <- apply(MXLENun, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))

ranx <- range(years)
rany <- range(ed1, ed.S)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
points(yrsUN1, ed.UN, col=4)

fitUN <- lm(ed.UN~yrsUN1)
ed.tar <- predict(fitUN, list(yrsUN1=yearsF))
ranx <- range(years)
rany <- range(ed1, ed.S)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
points(yearsF, ed.tar, col=5)
