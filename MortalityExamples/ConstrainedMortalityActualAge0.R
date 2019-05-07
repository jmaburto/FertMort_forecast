rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
library(MortalitySmooth)
library(HMDdata)
library(vars)
library(forecast)
library(MASS)
library(matrixcalc)


## select data
cou <- "Japan"
sex <- "F"
ages <- 0:105
m <- length(ages)
years1 <- 1960:2016
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
## w/o age 0
x0 <- x[-1]
m0 <- m-1

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
## observed e0
e01 <- apply(MX1, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
## observed e-dagger
ed1 <-  apply(MX1, 2, function(mx)  - sum(exp(-cumsum(mx))*(-cumsum(mx))))
## observed Gini
G1  <-  apply(MX1, 2, function(mx) 1 - sum(exp(-cumsum(mx))^2)/(sum(exp(-cumsum(mx))))-.0006)
# MX1vec <- c(MX1)
# In <- diag(n1)
# Im <- diag(m)
# C <- lower.tri(Im, diag = TRUE)
# C[C==1] <- -1 
# CC <- kronecker(In, C)
# ones.m <- matrix(1,1,m)
# Iones <- kronecker(In, ones.m)
# e01 <- Iones %*% exp(CC%*%MX1vec) + 0.5



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


## B-splines
library(MortalitySmooth)
## over ages, considering age 0
Bx0 <- MortSmooth_bbase(x=x0, xl=min(x0), xr=max(x0),
                        ndx=floor(m0/4), deg=3)
nbx0 <- ncol(Bx0)
Bx <- cbind(0, Bx0)
Bx <- rbind(c(1, rep(0,nbx0)), Bx)
# Bx <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
#                         ndx=floor(m/4), deg=3)
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


## overall basis, used only to compute summary measures
## as function of coefficients
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
Dx <- diff(diag(nbx), diff=2)
# Dx[1,1] <- 0
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

matplot(years, t(ETA.S), col=8, t="l", lty=1)
matpoints(years, t(ETA.S[c(1:10),]), col=rainbow(10))

ranx <- range(years)
for(i in 1:m){
  rany <- range(lMX[i,], ETA.S[i,], na.rm = TRUE, finite=TRUE)
  plot(years1, lMX1[i,], xlim=ranx, ylim=rany)
  lines(years, ETA.S[i,], col=2, lwd=2)
  title(main=ages[i])
  # locator(1)
  Sys.sleep(0.1)
}
matplot(x, ETA.S, col=8, lty=1, t="l")
for(i in 1:n){
  points(x, lMX[,i], pch=16, cex=0.8)
  lines(x, ETA.S[,i], col=2, lwd=2)
  title(main=years[i])
  # locator(1)
  Sys.sleep(0.1)
  points(x, lMX[,i], pch=16, cex=0.8, col="white")
  lines(x, ETA.S[,i], col="white", lwd=2)
  title(main=years[i], col.main="white")
}


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

## from UN WPP
yrsUN1 <- seq(2017.5, 2097.5, 5)
UNfile <- paste("e0", sex, ".csv", sep="")
e0UNall <- read.table(UNfile, sep=",", header=TRUE)
e0UN <- t(as.matrix(e0UNall[which(e0UNall[,1]==cou), 5:21]))
ranx <- range(years)
rany <- range(e01, e0.S, e0UN)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
abline(v=years[n1]+0.5)
matpoints(yrsUN1, e0UN, pch=16, col=c("green", "green4", "black", "blue4", "blue"))

## interpolate for each year
e0s.tar <- matrix(0, nF, 3)
for(i in 1:ncol(e0s.tar)){
  ww <- c(2,3,4)[i]
  bla <- c(years1,yrsUN1)
  fitUN <- loess(c(e01, e0UN[,ww]) ~ bla, span = 0.1)
  e0s.tar[,i] <- predict(fitUN, newdata=data.frame(bla=yearsF))
}
matlines(yearsF, e0s.tar, pch=17, col=c("green", "green4", "black", "blue4", "blue"))

# ## using time-series
# e0.ts <- ts(e01, start = t1[1])
# e0.mod <- auto.arima(e0.ts, d=1, max.p=1, max.q=1)
# pred.e0 <- forecast(e0.mod, h=nF, level=80)
# e0s.tar <- cbind(as.vector(pred.e0$lower), 
#                  as.vector(pred.e0$mean), 
#                  as.vector(pred.e0$upper))
# e0.tar <- e0s.tar[,2]
ranx <- range(years)
rany <- range(e01, e0.S, e0UN)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)


# ## future mean for e0
# e0mu <- c(e01[n1], e0s.tar[,3])
# ## future standard deviation for e0
# e0sig <- c(0.001, (e0s.tar[,5] - e0s.tar[,1])/4)
# ## product of each combination of elements of the sd vector
# e0var <- e0sig %*% t(e0sig)
# ## variance-covariance matrix
# e0Cor <- matrix(0.9999, nF+1, nF+1)
# diag(e0Cor) <- 1
# ## scalar multiply e0var by the variance-covariance matrix 
# ## obtaining covariance matrix of the variables
# e0Sig <- e0var * e0Cor
# is.positive.definite(e0Sig)
# ## generate target e0
# e0.tar <- mvrnorm(n = 1, mu=e0mu, Sigma=e0Sig)
# e0.tar <- e0.tar[-1]
e0.tar <- e0s.tar[,2]
lines(yearsF, e0.tar, col=4, lwd=2, lty=2)
matlines(yearsF, e0s.tar, col=3, lwd=2, lty=3)


## rough forecast of e-dagger
ed.ts <- ts(ed1, start = t1[1])
ed.mod <- auto.arima(ed.ts, d=1, max.p=1, max.q=1)
pred.ed <- forecast(ed.mod, h=nF, level=50)

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
eds.tar <- cbind(as.vector(pred.ed$lower), 
                 as.vector(pred.ed$mean), 
                 as.vector(pred.ed$upper))
ed.tar <- eds.tar[,2]
lines(yearsF, ed.tar, col=4, lwd=2, lty=2)
matpoints(yearsF, eds.tar, pch=17, col=c("green", "green4", "black", "blue4", "blue"))


## for combination of targeted upper and lower e0 and e-dagger

MUs.C <- ETAs.C <- array(0, dim=c(m,n,4))
BETASs.C <- array(0, dim=c(nbx,nbt,4))

for(l in 1:4){
  w0 <- c(1,3,1,3)
  wd <- c(1,1,3,3)
  e0.tar <- e0s.tar[,w0[l]]
  ed.tar <- eds.tar[,wd[l]]


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
  ## constraints on e-dagger over "future" years
  gxFUN <-  function(BETASF){
    ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
    MU.F <- exp(ETA.F)
    ed.F <- apply(MU.F, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
    gx <- ed.F - ed.tar
    return(gx)
  }
  gx1FUN <- function(BETASF){
    ETA.F <- MortSmooth_BcoefB(Bx, BtF, BETASF)
    MU.F <- exp(ETA.F)
    
    gx1 <- matrix(0, nbF, nF)
    k=i=1
    for(k in 1:nF){
      for(i in 1:nbF){
        Bi <- BF[1:m+(k-1)*m,i]
        muk <- MU.F[,k]
        p1 <- -(-muk[1]) * (-Bi[1]*muk[1])*exp(-muk[1]) - (-Bi[1]*muk[1])*exp(-muk[1])
        for(j in 2:m){
          p1j <- -(-sum(muk[1:j])) * (-sum(Bi[1:j]*muk[1:j]))*exp(-sum(muk[1:j])) -  (-sum(Bi[1:j]*muk[1:j]))*exp(-sum(muk[1:j]))
          p1 <- p1 + p1j
        }   
        gx1[i,k] <- p1
      }
    }
    return(gx1)
  }
  
  
  # 
  # mx <- c(0.0002748482, 0.0001851505, 0.0001337158, 0.0001028918)#MU.F[,1]#
  # - sum(exp(-cumsum(mx))*(-cumsum(mx)))
  # 
  # - (  exp(-mx[1])*(-mx[1]) + 
  #      exp(-mx[1]-mx[2])*(-mx[1]-mx[2]) + 
  #      exp(-mx[1]-mx[2]-mx[3])*(-mx[1]-mx[2]-mx[3]) + 
  #      exp(-mx[1]-mx[2]-mx[3]-mx[4])*(-mx[1]-mx[2]-mx[3]-mx[4]) )
  # 
  # 
  # 
  # mx[1] <- exp(b11 a1 + b12 a2 + b13 a3)
  # mx[2] <- exp(b21 a1 + b22 a2 + b23 a3)
  # mx[3] <- exp(b31 a1 + b32 a2 + b33 a3)
  # 
  # - (  exp(-exp(b11 a1 + b12 a2 + b13 a3))*(-exp(b11 a1 + b12 a2 + b13 a3)) + 
  #      exp(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3))*(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3)) + 
  #      exp(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3)-exp(b31 a1 + b32 a2 + b33 a3))*(-exp(b11 a1 + b12 a2 + b13 a3)-exp(b21 a1 + b22 a2 + b23 a3)-exp(b31 a1 + b32 a2 + b33 a3)) )
  # 
  # - (  exp(-exp(b11 a1 + b12 a2))*(-exp(b11 a1 + b12 a2)) + 
  #      exp(-exp(b11 a1 + b12 a2)-exp(b21 a1 + b22 a2))*(-exp(b11 a1 + b12 a2)-exp(b21 a1 + b22 a2)) + 
  #      exp(-exp(b11 a1 + b12 a2)-exp(b21 a1 + b22 a2)-exp(b31 a1 + b32 a2))*(-exp(b11 a1 + b12 a2)-exp(b21 a1 + b22 a2)-exp(b31 a1 + b32 a2)) )
  
  
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
    gx <- gxFUN(BETASF)
    gx1 <- gx1FUN(BETASF)
    
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
    if(max(abs(d))<1e-3) break
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

  MUs.C[,,l] <- MU.C
  ETAs.C[,,l] <- ETA.C
  BETASs.C[,,l] <- BETAS.C
}
  





ranx <- range(t)
for(i in 1:m){
  rany <- range(lMX[i,], ETA.S[i,], ETA.C[i,], na.rm = TRUE, finite=TRUE)
  plot(t1, lMX1[i,], xlim=ranx, ylim=rany)
  lines(t, ETA.S[i,], col=2, lwd=2)
  lines(t, ETA.C[i,], col=4, lwd=3)
  title(main=x[i])
  # locator(1)
  Sys.sleep(0.1)
}
for(i in 1:n){
  rany <- range(lMX[,i], ETA.S[,i], ETA.C, na.rm=TRUE, finite=TRUE)
  plot(x, lMX[,i], ylim=rany)
  lines(x, ETA.S[,i], col=2, lwd=2)
  lines(x, ETA.C[,i], col=4, lwd=3)
  title(main=years[i])
  # locator(1)
  Sys.sleep(0.1)
}



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
par(mfrow=c(1,2))
e0.C <- apply(MU.C, 2, function(mx) sum(exp(-cumsum(mx)))+0.5)
ranx <- range(years)
rany <- range(e01, e0.tar, e0.S)
plot(years1, e01, xlim=ranx, ylim=rany)
lines(years, e0.S, col=2, lwd=2, lty=2)
lines(years, e0.C, col=4, lwd=3)
points(yearsF, e0.tar, pch=16, t="b", lwd=2)
abline(v=years[n1]+0.5)
## compute constrained e-dagger
ed.C <- apply(MU.C, 2, function(mx) - sum(exp(-cumsum(mx))*(-cumsum(mx))))
ranx <- range(years)
rany <- range(ed1, ed.tar, ed.S)
plot(years1, ed1, xlim=ranx, ylim=rany)
lines(years, ed.S, col=2, lwd=2, lty=2)
lines(years, ed.C, col=4, lwd=3)
points(yearsF, ed.tar, pch=16, t="b", lwd=2)
abline(v=years[n1]+0.5)
par(mfrow=c(1,1))


round(cbind(years, e0.S, e0.C, e0.S-e0.C),2)
round(cbind(years, ed.S, ed.C, ed.S-ed.C),2)




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

