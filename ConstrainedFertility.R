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

# true linear predictor
ETA.T <- log(MU.T)

## exposures
E0 <- read.table("FRATNPexposRR.txt", header=TRUE, skip=2)
E <- matrix(E0$Exposure, m, n)
# E <- E/100

## Poisson expected values: true number of births
Y.T <- E * MU.T

## simulating births
Y <- matrix(rpois(m*n, c(Y.T)), m, n)

## simulated ASFR
MX <- Y/E
## in log
lMX <- log(Y/E)

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
  mx <- MX[,i]
  mab <- MAB[i]
  VAB[i] <- sum((((x+0.5) - mab)^2)*mx)/sum(mx)
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
Bt <- MortSmooth_bbase(x=t, xl=min(t), xr=max(t),
                       ndx=floor(n/4), deg=3)
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
lambdas <- c(10^2, 10^5)
lambda.x <- lambdas[1]
lambda.t <- lambdas[2]
P <- lambda.x * Px + lambda.t * Pt

## data in vector

y <- c(Y)
e <- c(E)
wei <- c(WEI)

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
eta.hat <- eta
mu.hat <- exp(eta)
ETA.hat <- matrix(eta.hat, m, n)
MU.hat <- matrix(mu.hat, m, n)


ranx <- range(t)
for(i in 1:m){
  rany <- range(MX1[i,], MU.hat[i,])
  plot(t1, MX1[i,], xlim=ranx, ylim=rany)
  lines(t, MU.hat[i,], col=2, lwd=2)
  title(main=x[i])
  # locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  if(i<=n1){
    rany <- range(MX1[,i], MU.hat[,i])
    plot(x, MX1[,i], ylim=rany)
    lines(x, MU.hat[,i], col=2, lwd=2)
  }else{
    plot(x, MU.hat[,i], t="l", col=2, lwd=2)
  }
  title(main=t[i])
  # locator(1)
  Sys.sleep(0.1)
}

## compute estimated and forecast summary measures
TFR <- colSums(MU.hat)
plot(t, TFR)
points(t, TFRobs, col=c(rep(2,n1), rep(3,n-n1)))

MAB <- apply(MU.hat, 2, 
             function(mx) sum(mx*(x+0.5)) / sum(mx))
plot(t, MAB)
points(t, MABobs, col=c(rep(2,n1), rep(3,n-n1)))

## compute variance age at childbearing
VAB <- numeric(n)
for(i in 1:n){
  mx <- MU.hat[,i]
  mab <- MAB[i]
  VAB[i] <- sum((((x+0.5) - mab)^2)*mx)/sum(mx)
}
plot(t, VAB)
points(t, VABobs, col=c(rep(2,n1), rep(3,n-n1)))


## constrining only on TFR











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
tDD <- t(D)%*D
