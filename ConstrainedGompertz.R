rm(list = ls())
options(device="X11")

## true Gompertz parameters
aT <- 0.00001
bT <- 0.1

## ages and dimension
m <- 70
x <- 1:m+30

## sex for the construction of the life-table (basically for the first a_x)
sexLT <- "F"

## take given exposures
library(MortalitySmooth)
e <- selectHMDdata("Sweden", "Exposures", "Females", ages=x, years=2000)
e <- as.vector(e)
## reduce/increase exposure to play on the variability
e <- e/10

## true linear predictor
etaT <- log(aT) + bT*x
## true force of mortality
muT <- exp(etaT)
## true expected values
dT <- e*muT
## simulated death counts
d <- rpois(m, dT)
## simulated rates
mx <- d/e
## simulated log-rates
lmx <- log(mx)

## plotting true and simulated log-mortality
plot(x, lmx)
lines(x, etaT, col=2)


## estimating gompertz with Poisson log-likelihood
lnLLKgom <- function(theta, death, expo, ages){
  a <- theta[1]
  b <- theta[2]
  bages <- b*ages
  mu <- a*exp(bages)
  LLK <- -sum(death*log(a) + death*bages - expo*mu)
  return(LLK)
}

opt <- optim(par=c(aT,bT), fn=lnLLKgom, death=d, expo=e, ages=x)

# ## alternative: estimating Gompertz as Poisson-GLM
# ## however lnLLKgom is needed afterwards
# fit <- glm(d~x, offset=log(e), family=poisson)

## fitted gompertz parameters
a.hat <- opt$par[1]
b.hat <- opt$par[2]
## fitted linear predictor and hazard
eta.hat <- log(a.hat) + b.hat*x
mu.hat <- exp(eta.hat)

## plotting simulated, true and fitted log-mortality
plot(x, lmx)
lines(x, etaT, col=2)
lines(x, eta.hat, col=3)

## functions for life-table stuff
source("LifeTableFUN.R")
## true life-table
LT.T <- lifetable.mx(x=x, mx=muT, sex=sexLT)
## true life expectancy
e0.T <- LT.T$ex[1]
## true e-dagger
ed.T <- eDagger(LT.T)[1]
## observed life-table
LT.obs <- lifetable.mx(x=x, mx=mx, sex=sexLT)
## observed life expectancy
e0.obs <- LT.obs$ex[1]
## observed e-dagger
ed.obs <- eDagger(LT.obs)[1]
## fitted life-table
LT.hat <- lifetable.mx(x=x, mx=mu.hat, sex=sexLT)
## fitted life expectancy
e0.hat <- LT.hat$ex[1]
## fitted e-dagger
ed.hat <- eDagger(LT.hat)[1]

## comparing
comp <- data.frame(true=c(e0.T, ed.T),
                   obs=c(e0.obs, ed.obs),
                   hat=c(e0.hat, ed.hat))
rownames(comp) <- c("e0", "ed")
comp


## library for optimization with non-linear constraints
library(Rsolnp)

## constrain function
constr <- function(theta, death, expo, ages){
  a <- theta[1]
  b <- theta[2]
  bages <- b*ages
  mu <- a*exp(bages)
  LT <- lifetable.mx(x=ages, mx=mu, sex=sexLT)
  e0 <- LT$ex[1]
  ed <- eDagger(LT)[1]
  return(c(e0,ed))
}

## both objective and constrain functions must have the same arguments
args(lnLLKgom)
args(constr)

## constraints vector
a <- c(e0.obs,ed.obs)

## starting value, unconstrained parameters
st.par <- c(opt$par)

## testing functions
lnLLKgom(st.par, death=d, expo=e, ages=x)
constr(st.par, death=d, expo=e, ages=x)
a

## consrained optimization
optC <- solnp(pars=st.par, fun=lnLLKgom, 
              eqfun=constr, eqB=a,
              death=d, expo=e, age=x)

## comparison
cbind(opt$par, optC$par)

## checking constraints
constr(opt$par, death=d, expo=e, ages=x)
constr(optC$par, death=d, expo=e, ages=x)
a

## fitted gompertz parameters
a.hatC <- optC$par[1]
b.hatC <- optC$par[2]
## fitted linear predictor and hazard
eta.hatC <- log(a.hatC) + b.hatC*x
mu.hatC <- exp(eta.hatC)

## plotting simulated, true and fitted log-mortality
plot(x, lmx)
lines(x, etaT, col=2)
lines(x, eta.hat, col=3)
lines(x, eta.hatC, col=4)

comp$hatC <- constr(optC$par, death=d, expo=e, ages=x)
comp







## smooth hazard

rm(list = ls())
options(device="X11")

## ages and dimension
m <- 70
x <- 1:m+30

## sex for the construction of the life-table (basically for the first a_x)
sexLT <- "F"

## take given exposures
library(MortalitySmooth)
e <- selectHMDdata("Sweden", "Exposures", "Females", ages=x, years=2000)
e <- as.vector(e)
## reduce/increase exposure to play on the variability
e <- e/10

## true linear predictor
aT <- 0.00001
bT <- 0.1
etaT <- log(aT) + bT*x + cos(x/15)
## true force of mortality
muT <- exp(etaT)
## true expected values
dT <- e*muT
## simulated death counts
d <- rpois(m, dT)
## simulated rates
mx <- d/e
## simulated log-rates
lmx <- log(mx)

## plotting true and simulated log-mortality
plot(x, lmx)
lines(x, etaT, col=2)

## estimating smooth curve with penalized Poisson log-likelihood
library(MortalitySmooth)
B <- MortSmooth_bbase(x=x, min(x), max(x), 10, 3)
nb <- ncol(B)
D <- diff(diag(nb), diff=2)

## estimating smooth hazard with penalized Poisson log-likelihood
lnPLLKnp <- function(theta, B, death, expo, D, lmbd){
  eta <- B%*%theta
  mu <- exp(eta)
  P <- lmbd * sum((D%*%theta)^2)
  LLK <- sum(death*eta - expo*mu)
  PLLK <- - (LLK - P)
  return(PLLK)
}

st.val <- solve(t(B)%*%B, t(B)%*%etaT)
lambda <- 10
opt <- optim(par=st.val, fn=lnPLLKnp, B=B, death=d, expo=e, D=D, lmbd=lambda)

## fitted B-splines coeff
theta.hat <- opt$par
## fitted linear predictor and hazard
eta.hat <- B%*%theta.hat
mu.hat <- exp(eta.hat)

## plotting simulated, true and fitted log-mortality
plot(x, lmx)
lines(x, etaT, col=2)
lines(x, eta.hat, col=3)

## functions for life-table stuff
source("LifeTableFUN.R")
## true life-table
LT.T <- lifetable.mx(x=x, mx=muT, sex=sexLT)
## true life expectancy
e0.T <- LT.T$ex[1]
## true e-dagger
ed.T <- eDagger(LT.T)[1]
## observed life-table
LT.obs <- lifetable.mx(x=x, mx=mx, sex=sexLT)
## observed life expectancy
e0.obs <- LT.obs$ex[1]
## observed e-dagger
ed.obs <- eDagger(LT.obs)[1]
## fitted life-table
LT.hat <- lifetable.mx(x=x, mx=mu.hat, sex=sexLT)
## fitted life expectancy
e0.hat <- LT.hat$ex[1]
## fitted e-dagger
ed.hat <- eDagger(LT.hat)[1]

## comparing
comp <- data.frame(true=c(e0.T, ed.T),
                   obs=c(e0.obs, ed.obs),
                   hat=c(e0.hat, ed.hat))
rownames(comp) <- c("e0", "ed")
comp


## library for optimization with non-linear constraints
library(Rsolnp)

## constraining function
constr <- function(theta, B, death, expo, D, lmbd){
  eta <- B%*%theta
  mu <- exp(eta)
  LT <- lifetable.mx(x=1:length(mu), mx=mu, sex=sexLT)
  e0 <- LT$ex[1]
  ed <- eDagger(LT)[1]
  return(c(e0,ed))
}

## both objective and constrain functions must have the same arguments
args(lnPLLKnp)
args(constr)

## constraints vector
a <- c(e0.obs,ed.obs)

## starting value, unconstrained parameters
st.par <- c(opt$par)

## testing functions
lnPLLKnp(st.par, B=B, death=d, expo=e, D=D, lmbd=lambda)
constr(st.par, B=B, death=d, expo=e)
a

## consrained optimization
optC <- solnp(pars=st.par, fun=lnPLLKnp, 
              eqfun=constr, eqB=a,
              B=B, death=d, expo=e, D=D, lmbd=lambda)

## comparison
cbind(opt$par, optC$par)

## checking constraints
constr(opt$par, B=B, death=d, expo=e, D=D, lmbd=lambda)
constr(optC$par, B=B, death=d, expo=e)
a

## fitted B-splines coeff
theta.hatC <- optC$par
## fitted linear predictor and hazard
eta.hatC <- B%*%theta.hatC
mu.hatC <- exp(eta.hatC)

## plotting simulated, true and fitted log-mortality
plot(x, lmx)
lines(x, etaT, col=2)
lines(x, eta.hat, col=3)
lines(x, eta.hatC, col=4)

comp$hatC <- constr(optC$par, B=B, death=d, expo=e)
comp














## smooth hazard and actual data

rm(list = ls())
options(device="X11")

## ages and dimension
x <- 1:100
m <- length(x)

## sex for the construction of the life-table (basically for the first a_x)
sexLT <- "F"

## take actual data
library(MortalitySmooth)
e <- selectHMDdata("Sweden", "Exposures", "Females", ages=x, years=2000)
e <- as.vector(e)
d <- selectHMDdata("Sweden", "Deaths", "Females", ages=x, years=2000)
d <- as.vector(d)

## actual rates
mx <- d/e
## actual log-rates
lmx <- log(mx)

## plotting actual log-mortality
plot(x, lmx)

## estimating smooth curve with penalized Poisson log-likelihood
library(MortalitySmooth)
B <- MortSmooth_bbase(x=x, min(x), max(x), floor(m/4), 3)
nb <- ncol(B)
D <- diff(diag(nb), diff=2)

## estimating smooth hazard with penalized Poisson log-likelihood
lnPLLKnp <- function(theta, B, death, expo, D, lmbd){
  eta <- B%*%theta
  mu <- exp(eta)
  P <- lmbd * sum((D%*%theta)^2)
  LLK <- sum(death*eta - expo*mu)
  PLLK <- - (LLK - P)
  return(PLLK)
}

## starting from a log-linear model
# fit0 <- lm(lmx ~ x)
# eta0 <- fit0$fitted
# or from plain smoothing
fit0 <- Mort1Dsmooth(x=x, y=d, offset=log(e))
eta0 <- fit0$logmortality
st.val <- solve(t(B)%*%B, t(B)%*%eta0)
lambda <- fit0$lambda
opt <- optim(par=st.val, fn=lnPLLKnp, B=B, death=d, expo=e, D=D, lmbd=lambda)

## fitted B-splines coeff
theta.hat <- opt$par
## fitted linear predictor and hazard
eta.hat <- B%*%theta.hat
mu.hat <- exp(eta.hat)

## plotting actual and fitted log-mortality
plot(x, lmx)
lines(x, eta.hat, col=2)

## functions for life-table stuff
source("LifeTableFUN.R")
## observed life-table
LT.obs <- lifetable.mx(x=x, mx=mx, sex=sexLT)
## observed life expectancy
e0.obs <- LT.obs$ex[1]
## observed e-dagger
ed.obs <- eDagger(LT.obs)[1]
## fitted life-table
LT.hat <- lifetable.mx(x=x, mx=mu.hat, sex=sexLT)
## fitted life expectancy
e0.hat <- LT.hat$ex[1]
## fitted e-dagger
ed.hat <- eDagger(LT.hat)[1]

## comparing
comp <- data.frame(obs=c(e0.obs, ed.obs),
                   hat=c(e0.hat, ed.hat))
rownames(comp) <- c("e0", "ed")
comp


## library for optimization with non-linear constraints
library(Rsolnp)

## constraining function
constr <- function(theta, B, death, expo, D, lmbd){
  eta <- B%*%theta
  mu <- exp(eta)
  LT <- lifetable.mx(x=1:length(mu), mx=mu, sex=sexLT)
  e0 <- LT$ex[1]
  ed <- eDagger(LT)[1]
  return(c(e0,ed))
}

## both objective and constrain functions must have the same arguments
args(lnPLLKnp)
args(constr)

## constraints vector
a <- c(e0.obs,ed.obs)

## starting value, unconstrained parameters
st.par <- c(opt$par)

## testing functions
lnPLLKnp(st.par, B=B, death=d, expo=e, D=D, lmbd=lambda)
constr(st.par, B=B, death=d, expo=e)
a

## consrained optimization
optC <- solnp(pars=st.par, fun=lnPLLKnp, 
              eqfun=constr, eqB=a,
              B=B, death=d, expo=e, D=D, lmbd=lambda)

## comparison
cbind(opt$par, optC$par)

## checking constraints
constr(opt$par, B=B, death=d, expo=e, D=D, lmbd=lambda)
constr(optC$par, B=B, death=d, expo=e)
a

## fitted B-splines coeff
theta.hatC <- optC$par
## fitted linear predictor and hazard
eta.hatC <- B%*%theta.hatC
mu.hatC <- exp(eta.hatC)

## plotting actual and fitted log-mortality, w/ and w/o constraints
plot(x, lmx)
lines(x, eta.hat, col=2, lwd=2)
lines(x, eta.hatC, col=3, lwd=2)

## constraintts
comp$hatC <- constr(optC$par, B=B, death=d, expo=e)
comp








## try with the Lee-Carter
rm(list = ls())
options(device="X11")


## domains and dimension
x <- 0:100
m <- length(x)
t <- 1950:2010
n <- length(t)

## sex for the construction of the life-table (basically for the first a_x)
sexLT <- "F"

## take given exposures
library(MortalitySmooth)
E <- selectHMDdata("Sweden", "Exposures", "Females", ages=x, years=t)
E <- as.matrix(E)
Y <- selectHMDdata("Sweden", "Deaths", "Females", ages=x, years=t)
Y <- as.matrix(Y)

## rates
MX <- Y/E
## log-rates
lMX <- log(MX)

matplot(x, lMX, t="l")


Da <- diff(diag(m), diff=2)
Dk <- diff(diag(n), diff=2)
lambdaa <- 10^6
lambdab <- 10^4
lambdak <- 10^0


## estimating LC hazard with penalized Poisson log-likelihood
lnPLLKlc <- function(theta, Y, E, 
                     Da, lambdaa, lambdab, Dk, lambdak){
  m <- nrow(Y)
  n <- ncol(Y)
  a <- theta[1:m]
  b <- theta[1:m+m]
  k <- theta[1:n+2*m]
  
  ETA <- matrix(a,m,n) + b%*%t(k)
  MU <- exp(ETA)
  LLK <- sum(Y*ETA - E*MU)
  
  Pa <- lambdaa * sum((Da%*%a)^2)
  Pb <- lambdab * sum((Da%*%b)^2)
  Pk <- lambdak * sum((Dk%*%k)^2)
  
  PLLK <- - (LLK - Pa - Pb - Pk)
  return(PLLK)
}

source("FUNCTIONSleecarter.r")
fit0 <- LCpoi(Dth=Y, Exp=E)
st.val <- c(fit0$Alpha, fit0$Beta, fit0$Kappa)

opt <- optim(par=st.val, fn=lnPLLKlc, Y=Y, E=E, 
             Da=Da, lambdaa=lambdaa,
             lambdab=lambdab,
             Dk=Dk, lambdak=lambdak,
             control=list(maxit=10^5))

a.hat <- opt$par[1:m]
b.hat <- opt$par[1:m+m]
k.hat <- opt$par[1:n+2*m]

par(mfrow=c(1,3))
plot(x, a.hat)
points(x, fit0$Alpha, col=2, pch=3)
plot(x, b.hat)
points(x, fit0$Beta, col=2, pch=3)
plot(t, k.hat)
points(t, fit0$Kappa, col=2, pch=3)




