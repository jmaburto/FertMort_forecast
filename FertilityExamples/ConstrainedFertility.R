rm(list = ls())
options(device="X11")
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
VAB <- numeric(n)
for(i in 1:n){
  mx.i <- MX[,i]
  mab.i <- MAB[i]
  VAB[i] <- sum((((x+0.5) - mab.i)^2)*mx.i)/sum(mx.i)
}

## plotting summary measures
par(mfrow=c(1,3))
plot(t, TFR)
plot(t, MAB)
plot(t, VAB)
par(mfrow=c(1,1))

# plot(VAB, TFR, t="n")
# text(x=VAB, y=TFR, labels=t-1900, col=colon)


## assuming we know only asfr
## up to 1996, but we have knowledge 
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

## keep aside all asfr we in theory know
Yobs <- Y
Eobs <- E
MXobs <- MX
lMXobs <- lMX
TFRobs <- TFR
MABobs <- MAB
VABobs <- VAB


## replace with arbitrary values all info
## after 1996, whatthe last year we assume to know
## this step is need to estimate up to 1996
## and simoultaneously forecast afterwards
Y <- matrix(10, m, n)
Y[1:m,1:n1] <- Y1
E <- matrix(10, m, n)
E[1:m,1:n1] <- E1
MX <- Y/E
MX[,(n1+1):n] <- NA

## set weights equal to zero
## to years 1997:2016
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
## select only for the observed years
Bt1 <- Bt[1:n1, 1:(nbt-3-1)]

## plotting(checking) B-splines bases over years
matplot(t, Bt, t="l", lwd=3, lty=3,
        col=rainbow(nbt))
matlines(t1, Bt1, t="l", lty=1, lwd=2,
         col=rainbow(nbt)[1:16])
abline(v=t[n1])

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
lambdas <- c(10^0, 10^5)
lambda.x <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.x * Px + lambda.t * Pt

## data in vector

y <- c(Y)
e <- c(E)
wei <- c(WEI)
mx <- c(MX)
## simple data-driven forecast
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
  rany <- range(MX1[i,], MU.S[i,])
  plot(t1, MX1[i,], xlim=ranx, ylim=rany)
  lines(t, MU.S[i,], col=2, lwd=2)
  lines(t, MU.T[i,], col=4, lty=2)
  title(main=x[i])
  # locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  if(i<=n1){
    rany <- range(MX1[,i], MU.S[,i],0 )
    plot(x, MX1[,i], ylim=rany)
    lines(x, MU.S[,i], col=2, lwd=2)
    lines(x, MU.T[,i], col=4, lty=2)
  }else{
    rany <- range(MU.S[,i], MU.T[,i],0 )
    plot(x, MU.S[,i], t="l", col=2, lwd=2,
         ylim=rany)
    lines(x, MU.T[,i], col=4, lty=2)
  }
  title(main=t[i])
  # locator(1)
  Sys.sleep(0.1)
}


plot3d(xx, tt, mx, type="p", col="red", 
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.T, 
          back='line', front='line', 
          col=5, lwd=1, alpha = 0.5)
surface3d(x, t, MU.S, 
          back='line', front='line', 
          col=2, lwd=1, alpha = 0.5)


## compute estimated and forecast summary measures
TFR.S <- colSums(MU.S)
plot(t, TFR.S)
points(t, TFRobs, col=c(rep(2,n1), rep(3,n-n1)))

MAB.S <- apply(MU.S, 2, 
               function(mx) sum(mx*(x+0.5)) / sum(mx))
plot(t, MAB.S)
points(t, MABobs, col=c(rep(2,n1), rep(3,n-n1)))

## compute variance age at childbearing
VAB.S <- numeric(n)
for(i in 1:n){
  mx <- MU.S[,i]
  mab <- MAB.S[i]
  VAB.S[i] <- sum((((x+0.5) - mab)^2)*mx)/sum(mx)
}
plot(t, VAB.S)
points(t, VABobs, col=c(rep(2,n1), rep(3,n-n1)))

## constraints on TFR ver time
hxFUN <-  function(betas){
  eta.it <- B%*%betas
  mu.it <- exp(eta.it)
  MU.it <- matrix(mu.it, m, n)
  TFR.it <- colSums(MU.it)
  hx <- TFR.it - TFRobs
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


ranx <- range(t)
for(i in 1:m){
  rany <- range(MX1[i,], MU.S[i,], MU.C[i,])
  plot(t1, MX1[i,], xlim=ranx, ylim=rany)
  lines(t, MU.S[i,], col=2, lwd=2)
  lines(t, MU.T[i,], col=4, lty=2)
  lines(t, MU.C[i,], col=3, lty=3, lwd=2)
  title(main=x[i])
  locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  if(i<=n1){
    rany <- range(MX1[,i], MU.S[,i],0 )
    plot(x, MX1[,i], ylim=rany)
    lines(x, MU.S[,i], col=2, lwd=2)
    lines(x, MU.T[,i], col=4, lty=2)
  }else{
    rany <- range(MU.S[,i], MU.T[,i],0 )
    plot(x, MU.S[,i], t="l", col=2, lwd=2,
         ylim=rany)
    lines(x, MU.T[,i], col=4, lty=2)
  }
  title(main=t[i])
  # locator(1)
  Sys.sleep(0.1)
}


plot3d(xx, tt, mx, type="p", col="red", 
       xlab="age", ylab="year", zlab="asfr", site=5, lwd=15)
surface3d(x, t, MU.T, 
          back='line', front='line', 
          col=5, lwd=1, alpha = 0.5)
surface3d(x, t, MU.S, 
          back='line', front='line', 
          col=2, lwd=1, alpha = 0.5)



















a <- round(colSums(ETA.T))
A <- kronecker(diag(n), matrix(1,1,m))
AB <- A%*%B

cbind(AB%*%betas.S, a)

eta <- log((y+1)/(e+1))
for(it in 1:100){
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(wei*mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P 
  tBWz <- t(B) %*% (w * z)
  
  LHS <- rbind(cbind(tBWBpP, t(AB)),
               cbind(AB, 0*diag(n)))
  RHS <- c(tBWz, a)
  
  betas <- solve(LHS, RHS)[1:nb]
  old.eta <- eta
  eta <- B %*% betas
  dif.eta <- max(abs(old.eta - eta))
  cat(it, dif.eta, "\n")
  if(dif.eta < 1e-4) break
}
  
cbind(AB%*%betas, a)

myfun <- function(betas){
  return(AB%*%betas - a)
}

eta <- log((y+1)/(e+1))
betas <- solve(t(B)%*%B, t(B)%*%eta)
for(it in 1:100){
  old.betas <- betas
  
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(wei*mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P 
  tBWz <- t(B) %*% (w * z)
  
  # LHS <- rbind(cbind(tBWBpP, t(AB)),
  #              cbind(AB, 0*diag(n)))
  # RHS <- c(tBWz, a)
  # betas <- solve(LHS, RHS)[1:nb]

  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP

  hx <- myfun(betas)
  
  LHS <- rbind(cbind(fx2, t(AB)),
               cbind(AB, 0*diag(n)))
  RHS <- c(-fx1, a)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  
  
  
  dif.betas <- max(abs(old.betas - betas))
  cat(it, dif.betas, "\n")
  if(dif.betas < 1e-4) break
}

cbind(AB%*%betas, a)





  
  betas.old <- betas
  betas <- solve(LHS, RHS)[1:nb]
  
  
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  
  Am <- hx1
  bv <- -hx
  
  # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
  #                 Amat=Am, bvec=bv, meq=n)
  # d <- opt$solution
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(n)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}







betas <- betas.S
## constraints on TFR ver time
hxFUN <-  function(betas){
  eta <- B%*%betas
  exp.eta <- exp(eta)
  exp.ETA <- matrix(exp.eta, m, n)
  sum.exp.ETA <- colSums(exp.ETA)
  hx <- sum.exp.ETA - TFRobs # !!!
  return(hx)
}


betas <- betas.S
hx1FUN <-  function(betas){
  eta <- B%*%betas
  exp.eta <- exp(eta)
  exp.ETA <- matrix(exp.eta, m, n)
  hx1 <- matrix(0, nb, n)
  for(j in 1:n){
    bla <- rep(0, m*n)
    bla[1:m+(j-1)*m] <- exp.ETA[,j]
    for(i in 1:nb){
      hx1[i,j] <- sum(B[,i] * bla)
    }
  }
  return(hx1)
}



betas <- betas.S
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
  
  # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
  #                 Amat=Am, bvec=bv, meq=n)
  # d <- opt$solution

  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(n)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
betas.C <- betas
eta.C <- eta
mu.C <- exp(eta.C)
ETA.C <- matrix(eta.C, m, n)
MU.C <- matrix(mu.C, m, n)

hxFUN(betas.C);hxFUN(betas.S)

TFR.C <- colSums(MU.C)
plot(t, TFR.S, cex=1.2)
points(t, TFRobs, col=c(rep(2,n1), rep(3,n-n1)))
points(t, TFR.C, col=4, pch=16, cex=1)









## simple unidimensional fit of asfr 
## with constraint on the TFR
## which is the sum of the asfr

rm(list = ls())
options(device="X11")

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

# a, b and c are associated with total fertility, height of
# the curve and Mean Age at Childbearing.
mu.T <- Hadwiger(1.25, 3.5, 28, x=x)

# true linear predictor
eta.T <- log(mu.T)

## exposures
e <- c(421660,415695,420043,415994,418983,
       413376,414893,414386,416835,413434,
       416274,415330,421366,429169,432178,
       427665,419672,361242,309296,308672,
       294653,274481,252437,290595,294925,
       298035,302802,305092,314129,318449,
       324701,332193,338598,329218,327132,
       328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/100
## true expected values
y.T <- e*mu.T

## simulating counts
y <- rpois(m, y.T)
## simulating rates
mx <- y/e
## simulating log-rates
lmx <- log(mx)

plot(x, lmx)
lines(x, eta.T)

## simple discrete smothing
Im <- diag(m)
D <- diff(Im, diff=2)
tDD <- t(D)%*%D



rm(list = ls())
options(device="X11")

library(MortalitySmooth)

x <- 12:55
m <- length(x)
B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                       ndx=floor(m/4), deg=3)
nb <- ncol(B)

betas <- c(-13.113,-6.686,-3.576,-2.123,-1.618,
           -1.564,-1.914,-2.445,-3.168,-3.96,
           -4.782,-5.807,-6.904,-8.297)

plot(x, exp(B%*%betas))
sum(exp(B%*%betas)) / sum(exp(betas))

sum(B)




