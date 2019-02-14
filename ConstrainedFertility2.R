rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)


Hadwiger <- function(a, b, c, x){
  ## function for computing asfr from
  ## Hadwiger (1940) model
  p1 <- a*b/c * (c/x)^(3/2)
  p2 <- -b^2 * (c/x + x/c -2)
  asfr <- p1*exp(p2)
  return(asfr)
}

x <- 12:55
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

## eventual weights
WEI <- matrix(0,m,n)
WEI[,1:61] <- 1

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
TFR.tar <- TFR#+0.2#runif(n, -0.5, 0.5)
plot(t, TFR, ylim=range(TFR, TFR.tar))
points(t, TFR.tar, col=2, pch=2)

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


## penalty stuff
Dx <- diff(diag(nbx), diff=2)
tDDx <- t(Dx)%*%Dx
Dt <- diff(diag(nbt), diff=2)
tDDt <- t(Dt)%*%Dt
Px <- kronecker(diag(nbt), tDDx)
Pt <- kronecker(tDDt, diag(nbx))
## smoothing parameters
lambdas <- c(10^1, 10^3.5)
lambda.x <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.x * Px + lambda.t * Pt

## data in vector
y <- c(Y)
e <- c(E)
wei <- c(WEI)
## simple smoothing
eta <- log((y+1)/(e+1))
for(it in 1:20){
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(wei*mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  betas <- solve(tBWBpP, tBWz)
  old.eta <- eta
  eta <- B %*% betas
  dif.eta <- max(abs(old.eta - eta))
  cat(it, dif.eta, "\n")
  if(dif.eta < 1e-4) break
}
betas.S <- betas
BETAS.S <- matrix(betas.S, nbx, nbt)
eta.S <- eta
mu.S <- exp(eta.S)
ETA.S <- matrix(eta.S, m, n)
MU.S <- matrix(mu.S, m, n)


ranx <- range(t)
for(i in 1:m){
  rany <- range(MX[i,], MU.S[i,])
  plot(t, MX[i,], xlim=ranx, ylim=rany)
  lines(t, MU.S[i,], col=2, lwd=2)
  lines(t, MU.T[i,], col=4, lty=2)
  title(main=x[i])
  # locator(1)
  Sys.sleep(0.1)
}
# for(i in 1:n){
#   rany <- range(MX[,i], MU.S[,i])
#   plot(x, MX[,i], ylim=rany)
#   lines(x, MU.S[,i], col=2, lwd=2)
#   lines(x, MU.T[,i], col=4, lty=2)
#   title(main=t[i])
#   # locator(1)
#   Sys.sleep(0.1)
# }


plot3d(xx, tt, mx, type="p", col="red",
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.T,
          back='line', front='line',
          col=5, lwd=1, alpha = 0.5)
surface3d(x, t, MU.S,
          back='line', front='line',
          col=2, lwd=1, alpha = 0.5)


## compute and plot estimated TFR
TFR.S <- colSums(MU.S)
plot(t, TFR, col=1, pch=1)
points(t, TFR.tar, col=2, pch=2)
points(t, TFR.S, col=3, pch=3)
legend("top", c("True", "Target", "Smooth"), col=1:3, pch=1:3)

## constraints on TFR ver time
hxFUN <-  function(betas){
  eta.it <- B%*%betas
  mu.it <- exp(eta.it)
  MU.it <- matrix(mu.it, m, n)
  TFR.it <- colSums(MU.it)
  hx <- TFR.it - TFR.tar
  return(hx)
}

# cbind(hxFUN(betas.S), TFR-TFR.S)

hx1FUN <-  function(betas){
  eta.it <- B%*%betas
  mu.it <- exp(eta)
  MU.it <- matrix(mu.it, m, n)
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

eta <- log((y+1)/(e+1))
betas <- solve(t(B)%*%B, t(B)%*%eta)
for(it in 1:100){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(wei*mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(n)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
}
betas.C <- betas
BETAS.C <- matrix(betas.C, nbx, nbt)
eta.C <- eta
mu.C <- exp(eta.C)
ETA.C <- matrix(eta.C, m, n)
MU.C <- matrix(mu.C, m, n)

TFR.C <- colSums(MU.C)
plot(t, TFR, col=1, pch=1, ylim=range(TFR, TFR.tar))
points(t, TFR.tar, col=2, pch=2)
points(t, TFR.S, col=3, pch=3)
points(t, TFR.C, col=4, pch=4)
legend("top", c("True", "Target", "Smooth", "Constrained"), 
       col=1:4, pch=1:4)

  
  
ranx <- range(t)
for(i in 1:m){
  rany <- range(MX[i,], MU.S[i,])
  plot(t, MX[i,], xlim=ranx, ylim=rany)
  lines(t, MU.S[i,], col=2, lwd=2)
  lines(t, MU.T[i,], col=4, lty=2)
  lines(t, MU.C[i,], col=3, lty=3, lwd=3)
  title(main=x[i])
  locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  rany <- range(MX[,i], MU.S[,i])
  plot(x, MX[,i], ylim=rany)
  lines(x, MU.S[,i], col=2, lwd=2)
  lines(x, MU.T[,i], col=4, lty=2)
  lines(x, MU.C[,i], col=3, lty=3, lwd=3)
  title(main=t[i])
  # locator(1)
  Sys.sleep(0.1)
}
  
  




















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
