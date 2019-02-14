rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
source("LifeTableFUN.R")
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

## compute life expectancy
e0 <- apply(MX, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
plot(t, e0)

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

## keep aside all asfr and summary measures 
## we in theory know
Yobs <- Y
Eobs <- E
MXobs <- MX
lMXobs <- lMX
## e0 only for future years, target e0
e0.tar <- e0[(n1+1):n]

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

## smoothing parameters
v1 <- c(rep(1,n1-1), rep(0,nF-1))
V1 <- diag(v1)
tDDt1 <- t(Dt)%*%V1%*%Dt
Pt1 <- kronecker(tDDt1, diag(nbx))
vF <- c(rep(0,n1-1), rep(1,nF-1))
VF <- diag(vF)
tDDtF <- t(Dt)%*%VF%*%Dt
PtF <- kronecker(tDDtF, diag(nbx))
lambda.tF <- 10^-1
PF <- lambda.x*Px + lambda.t*Pt1 + lambda.tF*PtF

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
#   rany <- range(lMX[,i], ETA.S[,i], na.rm=TRUE, finite=TRUE)
#   plot(x, lMX[,i], ylim=rany,  t="n")
#   if(i<=n1) points(x, lMX1[,i])
#   lines(x, ETA.S[,i], col=2, lwd=2)
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
plot(t, e0)
points(t, e0.S, col=2)
abline(v=t[n1]+01)


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


hx1FUN <-  function(BETASF){
  ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
  MU.F <- exp(ETA.F)
  exp.sum.MU.F <- -exp(colSums(-MU.F))

  hx1 <- matrix(0, nbF, nF)
  for(k in 1:nF){
    for(i in 1:nbF){
      Bi <- BF[1:m+(k-1)*m,i]
      mui <- MU.F[,k]
      etai <- ETA.F[,k]
      p1 <- Bi[1]*exp(mui[2] + etai[1]) + Bi[1]*mui[1] + Bi[2]*mui[2]
      for(j in 3:m){
        p1 <- p1*exp(mui[j]) + sum(Bi[1:j]*mui[1:j])
      }
      hx1[i,k] <- p1*exp.sum.MU.F[k]
    }
  }
  return(hx1)
}

BETAS <- BETAS.S
for(it in 1:100){
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
}  
BETAS.C <- BETAS
ETA.C <- ETA
MU.C <- exp(ETA)















































Hadwiger <- function(a, b, c, x){
  ## function for computing asfr from
  ## Hadwiger (1940) model
  p1 <- a*b/c * (c/x)^(3/2)
  p2 <- -b^2 * (c/x + x/c -2)
  asfr <- p1*exp(p2)
  return(asfr)
}

persp3d(x, t, ETA.T, col=5, alpha=0.8)






x <- 0:110
m <- length(x)
t <- 1946:2016
n <- length(t)
colon <- rainbow(n)
colom <- rainbow(m)

# a, b and c are associated with total fertility, height of
# the curve and Mean Age at Childbearing.

## true series of parameters
PAR <- read.table("PAR.Hadwiger.txt", header=TRUE)

## true underlying fertility over age and time
MU.T <- matrix(0, m, n)
for(i in 1:n){
  MU.T[,i] <- Hadwiger(PAR[i,1], PAR[i,2], PAR[i,3], x=x)
}

matplot(x, MU.T, t="l", col=colon)

persp3d(x, t, MU.T, 
        col=5, alpha=0.8)


# true linear predictor
ETA.T <- log(MU.T)

## exposures
E0 <- read.table("FRATNPexposRR.txt", header=TRUE, skip=2)
E <- matrix(E0$Exposure, m, n)
E <- E/20
e <- c(E)
## Poisson expected values: true number of births
Y.T <- E * MU.T

## simulating births
y <- rpois(m*n, c(Y.T))
Y <- matrix(y, m, n)

## simulated ASFR
MX <- Y/E
mx <- c(MX)
## in log
lMX <- log(Y/E)
lmx <- c(lMX)

whi <- floor(seq(1,n,length=4))
rany <- range(MX, MU.T)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(x, MX[,whi[i]], pch=16, 
       col=colon[whi[i]], main=paste(t[whi[i]]),
       ylim=rany)
  lines(x, MU.T[,whi[i]], col=colon[whi[i]], lwd=2)
}
par(mfrow=c(1,1))
xx <- rep(x, n)
tt <- rep(t, each=m)

plot3d(xx, tt, mx, type="p", col="red", 
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.T, 
          back='line', front='line', 
          col=5, lwd=1, alpha = 0.5)


## compute observed TFR
TFR <- colSums(MX)
plot(t, TFR)

## compute mean age at childbearing
MAB <- apply(MX, 2, 
             function(mx) sum(mx*(x+0.5)) / sum(mx))
plot(t, MAB)

## compute variance age at childbearing
## compute variance age at childbearing
VAB <- numeric(n)
for(i in 1:n){
  mx.i <- MX[,i]
  mab.i <- MAB[i]
  VAB[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}
plot(t, VAB)


## plotting summary measures
par(mfrow=c(1,3))
plot(t, TFR)
plot(t, MAB)
plot(t, VAB)
par(mfrow=c(1,1))

## assuming we know only asfr
## up to max(t1) (e.g. 1996), but we have knowledge 
## of TFR, MAB and VAB up to 2016 (actually computed)
## and we wanna forecast all asfr up to
## 2016 by constraining each "future" age-pattern to 
## follow the previously mentioned summary measures

t1 <- t[1]:1996
n1 <- length(t1)

Y1 <- Y[,1:n1]
E1 <- E[,1:n1]
MX1 <- Y1/E1
lMX1 <- log(MX1)

## keep aside all asfr and summary measures 
## we in theory know
Yobs <- Y
Eobs <- E
MXobs <- MX
lMXobs <- lMX
TFR.tar <- TFR
MAB.tar <- MAB
VAB.tar <- VAB

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
Bt <- diag(n)#MortSmooth_bbase(x=t, xl=min(t), xr=max(t),
            #           ndx=floor(n/4), deg=3)
nbt <- ncol(Bt)

## overall basis
B <- kronecker(Bt, Bx)
nb <- ncol(B)

Bx1 <- kronecker(matrix(1, ncol = nbx, nrow = 1), Bx)
Bx2 <- kronecker(Bx, matrix(1, ncol = nbx, nrow = 1))
RTBx <- Bx1 * Bx2
Bt1 <- kronecker(matrix(1, ncol = nbt, nrow = 1), Bt)
Bt2 <- kronecker(Bt, matrix(1, ncol = nbt, nrow = 1))
RTBt <- Bt1 * Bt2

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
# eta <- log((y+1)/(e+1))
# for(it in 1:20){
#   mu <- e*exp(eta)
#   z <- (y - mu)/mu + eta
#   w <- c(wei*mu)
#   tBWB <- t(B) %*% (w * B)
#   tBWBpP <- tBWB + P
#   tBWz <- t(B) %*% (w * z)
#   betas <- solve(tBWBpP, tBWz)
#   old.eta <- eta
#   eta <- B %*% betas
#   dif.eta <- max(abs(old.eta - eta))
#   cat(it, dif.eta, "\n")
#   if(dif.eta < 1e-4) break
# }
BETAS.S <- BETAS
ETA.S <- ETA
MU.S <- exp(ETA)


# ranx <- range(t)
# for(i in 1:m){
#   rany <- range(MX[i,], MU.S[i,], na.rm = TRUE)
#   plot(t, MX[i,], xlim=ranx, ylim=rany)
#   lines(t, MU.S[i,], col=2, lwd=2)
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


plot3d(xx, tt, mx, type="p", col="red",
       xlab="age", ylab="year", zlab="asfr", 
       site=5, lwd=15)
# surface3d(x, t, MU.T,
#           back='line', front='line',
#           col=5, lwd=1, alpha = 0.5)
surface3d(x, t, MU.S,
          back='line', front='line',
          col=2, lwd=1, alpha = 0.5)


## compute estimated TFR, MAB and VAB
TFR.S <- colSums(MU.S)
MAB.S <- apply(MU.S, 2, function(mx) sum(mx*(x+0.5)) / sum(mx))
VAB.S <- numeric(n)
for(i in 1:n){
  mx.i <- MU.S[,i]
  mab.i <- MAB.S[i]
  VAB.S[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}


## estimated TFR, MAB and VAB
par(mfrow=c(1,3))
plot(t, TFR.S)
points(t, TFR.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", "Target"), col=1:2, pch=c(1,16))

plot(t, MAB.S)
points(t, MAB.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", "Target"), col=1:2, pch=c(1,16))

plot(t, VAB.S)
points(t, VAB.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", "Target"), col=1:2, pch=c(1,16))
par(mfrow=c(1,1))

## constraints on TFR over time
hxFUN <-  function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  TFR.it <- colSums(MU.it)
  hx <- TFR.it - TFR.tar
  return(hx)
}
hx1FUN <-  function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  hx1 <- matrix(0, nb, n)
  for(j in 1:n){
    MU.it.j <- rep(0, m*n)
    MU.it.j[1:m+(j-1)*m] <- MU.it[,j]
    for(i in 1:nb){
      hx1[i,j] <- sum(B[,i] * MU.it.j)
    }
  }
  return(hx1)
}

## constrained MAB over time
gxFUN <- function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  x05 <- x+0.5
  MAB.t <- apply(MU.it, 2, function(mx)sum(mx*x05)/sum(mx))
  gx <- MAB.t - MAB.tar
  return(gx)
}

gx1FUN <- function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  x05 <- x+0.5
  
  den <- colSums(MU.it)
  den2 <- den^2
  MUx.it <- MU.it*x05
  mux.it <- c(MUx.it)
  sum.MUx.it <- colSums(MUx.it)
  
  gx1 <- matrix(0, nb, n)
  for(j in 1:n){
    whi1 <- rep(0,n*m)
    whi1[1:m+(j-1)*m] <- 1
    for(i in 1:nb){
      Bmux.i <- B[,i]*mux.it
      part1 <- sum(Bmux.i*whi1)/den[j]
      Bmu <- B[,i]*c(MU.it)
      part2 <- ( sum.MUx.it[j] * sum(Bmu*whi1) ) / den2[j]
      gx1[i,j] <- part1 - part2  
    }
  }
  return(gx1)
}

## constrained VAB over time
uxFUN <- function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  x05 <- x+0.5
  TFR.it <- colSums(MU.it)
  MAB.it <- apply(MU.it, 2, function(mx)sum(mx*x05)/sum(mx))
  VAB.t <- numeric(n)
  for(i in 1:n){
    VAB.t[i] <- sum(((x05 - MAB.it[i])^2)*MU.it[,i])/TFR.it[i]
  }
  ux <- VAB.t - VAB.tar
  return(ux)
}

ux1FUN <- function(BETAS){
  ETA.it <- MortSmooth_BcoefB(Bx, Bt, BETAS)
  MU.it <- exp(ETA.it)
  x05 <- x+0.5
  mu.it <- c(MU.it)
  
  TFR.it <- colSums(MU.it)
  TFR2.it <- TFR.it^2
  MAB.it <- apply(MU.it, 2, function(mx)sum(mx*x05)/sum(mx))
  
  x.mab <- sweep(matrix(x05,m,n), 2, MAB.it)
  x.mab2 <- x.mab^2
  x.mab2.MU <- x.mab2*MU.it
  sum.x.mab2.mu <- colSums(x.mab2.MU)
  x.mab2.mu <- c(x.mab2.MU)
  
  ux1 <- matrix(0, nb, n)
  for(j in 1:n){
    whi1 <- rep(0,n*m)
    whi1[1:m+(j-1)*m] <- 1
    for(i in 1:nb){
      B.x.mab2.mu <- B[,i]*x.mab2.mu
      part1 <- sum(B.x.mab2.mu*whi1) / TFR.it[j]
      part2 <- (sum(B[,i]*mu.it*whi1) * sum.x.mab2.mu[j])/TFR2.it[j]
      ux1[i,j] <- part1 - part2  
    }
  }
  return(ux1)
}


BETAS <- BETAS.S
for(it in 1:100){
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
  
  hx <- hxFUN(BETAS)
  hx1 <- hx1FUN(BETAS)
  
  gx <- gxFUN(BETAS)
  gx1 <- gx1FUN(BETAS)
  
  #ux <- uxFUN(BETAS)
  #ux1 <- ux1FUN(BETAS)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(2*n)))
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



## estimated, TFR, MAB and VAB
par(mfrow=c(1,3))
plot(t, TFR.S)
points(t, TFR.C, col=4, pch=3, lwd=2, cex=1.2)
points(t, TFR.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", 
                "Smooth+Constrained+Forecast", 
                "Target"), col=c(1,4,2), pch=c(1,3,16),
       pt.cex=c(1,1.2,1), lwd=c(1,2,1), lty=NA)

plot(t, MAB.S)
points(t, MAB.C, col=4, pch=3, lwd=2, cex=1.2)
points(t, MAB.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", 
                "Smooth+Constrained+Forecast", 
                "Target"), col=c(1,4,2), pch=c(1,3,16),
       pt.cex=c(1,1.2,1), lwd=c(1,2,1), lty=NA)

plot(t, VAB.S)
points(t, VAB.C, col=4, pch=3, lwd=2, cex=1.2)
points(t, VAB.tar, col=2, pch=c(rep(16,n1), rep(17,n-n1)))
abline(v=t[n1]+0.5)
legend("top", c("Smooth+Forecast", 
                "Smooth+Constrained+Forecast", 
                "Target"), col=c(1,4,2), pch=c(1,3,16),
       pt.cex=c(1,1.2,1), lwd=c(1,2,1), lty=NA)
par(mfrow=c(1,1))


ranx <- range(t)
for(i in 1:m){
  rany <- range(MX[i,], MU.S[i,], na.rm = TRUE)
  plot(t, MX[i,], xlim=ranx, ylim=rany)
  lines(t, MU.S[i,], col=2, lwd=2)
  lines(t, MU.T[i,], col=4, lty=2)
  lines(t, MU.C[i,], col=3, lty=3, lwd=3)
  title(main=x[i])
  locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  rany <- range(MX[,i], MU.S[,i], MU.C, na.rm = TRUE)
  plot(x, MX[,i], ylim=rany)
  lines(x, MU.S[,i], col=2, lwd=2)
  lines(x, MU.T[,i], col=4, lty=2)
  lines(x, MU.C[,i], col=3, lty=3, lwd=3)
  title(main=t[i])
  # locator(1)
  Sys.sleep(0.1)
}
  
whi <- floor(seq(1,n,length=9))
rany <- range(MX, MU.T, MU.S, MU.C, 0, na.rm = TRUE)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(x, MX[,whi[i]], pch=16, 
       col=1, main=paste(t[whi[i]]),
       ylim=rany)
  lines(x, MU.T[,whi[i]], col=2, lwd=2)
  lines(x, MU.S[,whi[i]], col=3, lwd=2)
  lines(x, MU.C[,whi[i]], col=4, lwd=2)
  legend("top", c("Simul", "True", "Smooth", "Constrained"), 
         col=1:4, lwd=2, lty=1, pch=c(16,NA,NA,NA))
}
par(mfrow=c(1,1))




plot3d(xx, tt, mx, type="p", col="red",
       xlab="age", ylab="year", zlab="asfr", 
       site=5, lwd=15)
surface3d(x, t, MU.S,
          back='line', front='line',
          col=5, lwd=1, alpha = 0.5)
surface3d(x, t, MU.C,
          back='line', front='line',
          col=2, lwd=1, alpha = 0.5)

















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

