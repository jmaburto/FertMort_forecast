rm(list = ls())
options(device="X11")

x <- 1:4
X <- cbind(1, x)
y <- c(6,5,7,10)

tXX <- t(X)%*%X
betas <- solve(tXX, t(X)%*%y)

exp(betas[1]) + exp(betas[2])

betas1 <- betas



LLKfun <- function(betas, X, y){
  y.hat <- X%*%betas
  RSS <- sum( (y-y.hat)^2 )
  return(RSS)
}


st.val <- c(3, 1)
opt <- optim(st.val, LLKfun, X=X, y=y)
opt$par

CONfun <- function(betas, X, y){
  A <- c(exp(betas[1]) + exp(betas[2]), 
         (1/2)*(betas[1]-betas[2])^2)
  return(A)
}

a <- c(40,2.3)

library(Rsolnp)

optC <- solnp(opt$par, fun=LLKfun, 
              eqfun=CONfun, eqB=a, 
              X=X, y=y)
optC$par
opt$par


CONfun(optC$pars, X=X, y=y)
CONfun(opt$par, X=X, y=y)

y.hat <- X%*%opt$par
y.hatC <- X%*%optC$pars
plot(x, y)
lines(x, y.hat, col=2)
lines(x, y.hatC, col=3)

## manual

## about OLS
# betas <- c(4,1)
# for(it in 1:10){
#   y.hat <- X%*%betas
#   r <- y.hat - y
#   oldbetas <- betas
#   betas <- betas - (solve(t(X)%*%X) %*% t(X)%*%r)
#   conv <- max(abs(oldbetas-betas))
#   cat(betas, conv, "\n")
#   if(conv< 10^-4) break
# }
# betas1


## about constraints
JACfun <- function(betas){
  gr1 <- c(exp(betas[1]), betas[1]-betas[2])
  gr2 <- c(exp(betas[2]), betas[2]-betas[1])
  GR <- cbind(gr1, gr2)
  return(GR)
}
Sfun <- function(betas){
  A <- c(exp(betas[1]) + exp(betas[2]), 
         (1/2)*(betas[1]-betas[2])^2)
  return(A)
}
# OBJfun <- function(betas, a){
#   s <- Sfun(betas)
#   sums <- sum((s-a)^2)
#   return(sums)
# }
# betas <- c(3, 1)
# 
# opt <- optim(betas, OBJfun, a=a)
# 
# for(it in 1:10){
#   A.cur <- Sfun(betas)
#   s <- A.cur - a
#   J <- JACfun(betas)
#   oldbetas <- betas
#   betas <- betas - (solve(t(J)%*%J) %*% t(J)%*%s)
#   conv <- max(abs(oldbetas-betas))
#   cat(betas, conv, "\n")
#   if(conv< 10^-4) break
# }
# betas
# opt$par


## both OLS and constraints
betas <- opt$par
for(it in 1:10){
  cat(it, "\n")
  y.hat <- X%*%betas
  r <- y.hat - y
  oldbetas <- betas
  betas <- betas - (solve(t(X)%*%X) %*% t(X)%*%r)
  cat("OLS:", betas, "\n")
  A.cur <- Sfun(betas)
  s <- A.cur - a
  J <- JACfun(betas)
  oldbetas <- betas
  betas <- betas - (solve(t(J)%*%J) %*% t(J)%*%s)
  cat("CON:", betas, "\n\n\n")
}
betas
optC$pars
CONfun(optC$pars, X=X, y=y)




## another example## another example## another example
## another example## another example## another example
## another example## another example## another example

rm(list = ls())
options(device="X11")


# f(x) = 6*(x1/x2) + x2/x1
# s.t.
# h(x) = x1*x2 - 2 = 0
# g(x) = 1 - x1 - x2 <= 0

fxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx <- 6*(x1/x2) + x2/(x1^2)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gr1 <- 6/x2 - 2*(x2/(x1^3))
  gr2 <- -6*(x1/(x2^2)) + 1/(x1^2)
  fx1 <- rbind(gr1,gr2)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gr11 <- 6*(x2/(x1^4))
  gr12 <- -6/(x2^2) - 2/(x1^3)
  gr21 <- -6/(x2^2) - 2/(x1^3)
  gr22 <- 12*x1/(x2^3)
  fx2 <- rbind(cbind(gr11,gr12),
               cbind(gr21,gr22))
  return(fx2)
}

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- x1*x2-2
  return(hx)
}
  
hx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(x2,x1)
  return(hx1)
}

hx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx2 <- rbind(cbind(0,1),
               cbind(1,0))
}

gxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx <- 1-x1-x2
  return(gx)
}

gx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx1 <- rbind(-1,-1)
  return(gx1)
}

gx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx2 <- rbind(cbind(0,0),
               cbind(0,0))
  return(gx2)
}

## numerical constrained minimization
library(Rsolnp)

st.val <- c(1, 1)
optnum <- solnp(pars=st.val, fun=fxFUN, 
                eqfun=hxFUN, eqB=0,
                ineqfun=gxFUN, 
                ineqUB=0, 
                ineqLB=-10^8)
xhat <- optnum$pars
xhat
fxFUN(xhat)
hxFUN(xhat)
gxFUN(xhat)


## graphical view
x1s <- seq(-2, 2, length=1000)
x1s <- x1s[x1s!=0]
l1 <- length(x1s)
x2s <- seq(-1, 3, length=999)
x2s <- x2s[x2s!=0]
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}


## delete below ineq constr
whiNA <- matrix(NA, l1, l2)
for(i in 1:l1){
  for(j in 1:l2){
    whiNA[i,j] <- ifelse(1-x1s[i]-x2s[j]<=0, 1, NA)
  }
}

FXna <- FX*whiNA

brek <- quantile(FXna, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FXna, breaks=brek, col=colo)

pmin <- which(FXna==min(FXna, na.rm=TRUE), arr.ind=TRUE)
xhat1 <- c(x1s[pmin[1]], x2s[pmin[2]])

hxFUN(xhat1)

## eq constr
x1s.eq <- 2/x2s
points(x1s.eq, x2s, col=4, pch=16, cex=0.5)
## ineq constr
x1s.ineq <- 1-x2s
lines(x1s.ineq, x2s, col=3, lwd=4)
points(xhat[1], xhat[2], pch=16, col=1, cex=2)

## optimization by SQP

## starting x
library(quadprog)
## starting x
x <- c(1,1)
for(it in 1:10){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  hx2 <- hx2FUN(x)
  gx <- gxFUN(x)
  gx1 <- gx1FUN(x)
  gx2 <- gx2FUN(x)

  A <- t(rbind(t(hx1), t(gx1)))
  a <- c(-hx, gx)
  opt <- solve.QP(Dmat=fx2, dvec=-fx1,
                  Amat=A, bvec=a, meq=1)
  d <- opt$solution

  x <- x + d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}




## STILL another example## STILL another example
## STILL another example## STILL another example
## STILL another example## STILL another example
## STILL another example## STILL another example
## which doesn't really work, with SQL

rm(list = ls())
options(device="X11")


# f(x) = sin(x1)cos(x2) + cos(x1)sin(x2)
# s.t.
# h(x) : x1 - x2^3 = 0
# g(x) : 0 <= x1 + x2 <= pi

# x1=0.688
# x2=0.883

fxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx <- - sin(x1)*cos(x2) - cos(x1)*sin(x2)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gr1 <- sin(x2)*sin(x1) - cos(x2)*cos(x1)
  gr2 <- sin(x1)*sin(x2) - cos(x1)*cos(x2)
  fx1 <- rbind(gr1,gr2)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gr11 <- cos(x2)*sin(x1) + sin(x2)*cos(x1)
  gr12 <- cos(x1)*sin(x2) + sin(x1)*cos(x2)
  gr21 <- cos(x2)*sin(x1) + sin(x2)*cos(x1)
  gr22 <- cos(x1)*sin(x2) + sin(x1)*cos(x2)
  fx2 <- rbind(cbind(gr11,gr12),
               cbind(gr21,gr22))
  return(fx2)
}

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- x1 - x2^3
  return(hx)
}

hx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(1,-3*x2^2)
  return(hx1)
}

hx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx2 <- rbind(cbind(0,0),
               cbind(0,-6*x2))
}

gxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx <- x1 + x2
  return(gx)
}

gx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx1 <- rbind(1,1)
  return(gx1)
}

gx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  gx2 <- rbind(cbind(0,0),
               cbind(0,0))
  return(gx2)
}


## numerical constrained minimization
library(Rsolnp)

# x1=0.688
# x2=0.883
st.val <- c(0.5, 1)
optnum <- solnp(pars=st.val, fun=fxFUN, 
                eqfun=hxFUN, eqB=0,
                ineqfun=gxFUN, 
                ineqUB=pi, 
                ineqLB=0)
xhat <- optnum$pars
xhat
fxFUN(xhat)
hxFUN(xhat)
gxFUN(xhat)


## graphical view
x1s <- seq(-2, 2, length=1000)
l1 <- length(x1s)
x2s <- seq(-1, 3, length=999)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}


## delete below ineq constr
whiNA <- matrix(NA, l1, l2)
for(i in 1:l1){
  for(j in 1:l2){
    whiNA[i,j] <- ifelse(x1s[i]+x2s[j]>=0 & x1s[i]+x2s[j]<=pi, 1, NA)
  }
}

FXna <- FX*whiNA

brek <- quantile(FXna, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FXna, breaks=brek, col=colo)

pmin <- which(FXna==min(FXna, na.rm=TRUE), arr.ind=TRUE)
xhat1 <- c(x1s[pmin[1]], x2s[pmin[2]])

hxFUN(xhat1)

## eq constr
x1s.eq <- x2s^3
points(x1s.eq, x2s, col=4, pch=16, cex=0.5)
## ineq constr
x1s.ineqL <- -x2s
lines(x1s.ineqL, x2s, col=3, lwd=4)
x1s.ineqU <- pi-x2s
lines(x1s.ineqU, x2s, col=3, lwd=4)
points(xhat[1], xhat[2], pch=16, col=1, cex=2)


## optimization by SQP

## starting x
library(quadprog)
library(Matrix)
## starting x
x <- c(0.6, 0.8, 0.1)
x1 <- x[1]
x2 <- x[2]
Z <- sin(x1)*cos(x2) + cos(x1)*sin(x2)
lambda <- x[3]
D2 <- rbind(c(-Z, -Z, 1),
            c(-Z, -(Z+6*lambda*x2), -3*(x2^2)),
            c(1, -3*(x2^2), 0))
D <- rbind(cos(x1)*cos(x2) - sin(x1)*sin(x2) + lambda,
           cos(x1)*cos(x2) - sin(x1)*sin(x2) + 2*lambda*(x2^2),
           x1-x2^3)

d <-solve(D2, D)

for(it in 1:10){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)+1e-2*diag(2)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  hx2 <- hx2FUN(x)
  gx <- gxFUN(x)
  gx1 <- gx1FUN(x)
  gx2 <- gx2FUN(x)
  
  A <- t(rbind(t(hx1), t(gx1)))
  a <- c(-hx, gx)
  opt <- solve.QP(Dmat=fx2, dvec=-fx1,
                  Amat=A, bvec=a, meq=1)
  x <- opt$solution
  
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}







###########################################
## OLS + nonlinear constraint 
## 11111111111111111111111111111111111111111

rm(list = ls())
options(device="X11")
library(quadprog)

xx <- sort(runif(100))
A <- cbind(1,xx)
xUC.T <- c(1,2)
y <- A%*%xUC.T + rnorm(length(xx),0,0.1)
plot(xx, y)
lines(xx, A%*%xUC.T, col=2)

fxFUN <- function(x){
  y.hat <- A%*%x
  #fx <- sum( (y-y.hat)^2 )
  fx <- t(y-y.hat)%*%(y-y.hat)
  fx <- as.numeric(fx)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx1 <- 2 * (t(A)%*%A%*%x) - 2*(t(A)%*%y)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx2 <- 2*(t(A)%*%A)
  return(fx2)
}

exp(xUC.T[1]) + exp(xUC.T[2])

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- exp(x1) + exp(x2) - 9 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(exp(x1), exp(x2))
  return(hx1)
}

x0 <- solve(t(A)%*%A, t(A)%*%y)
lines(xx, A%*%xUC.T, col=2, lty=2)
lines(xx, A%*%x0, col=4)


x <- x0
for(it in 1:20){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
 
  Am <- hx1
  bv <- -hx
  # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
  #                 Amat=Am, bvec=bv, meq=1)
  # d <- opt$solution
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:2]
  
  x <- x + d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(x);hxFUN(xUC.T)
fxFUN(x);fxFUN(xUC.T)

plot(xx, y)
lines(xx, A%*%x, col=2, lwd=2, lty=2)
lines(xx, A%*%x0, col=3, lwd=2, lty=1)



## graphical view
x1s <- seq(0.5, 1.6, length=100)
l1 <- length(x1s)
x2s <- seq(1, 2.19, length=99)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}


brek <- quantile(FX, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FX, breaks=brek, col=colo)
## eq constr
x2s.eq <- seq(1, 2.19, length=1000)
x1s.eq <- log(9 - exp(x2s.eq)) #### !!!!!!!!!!
lines(x1s.eq, x2s.eq, col=4, pch=16, cex=0.5)

## SQL fit
points(x[1], x[2], pch=16, col=1, cex=2)
## uncostrained fit
points(x0[1], x0[2], pch=17, col=3, cex=2)

cbind(x0,x)


## numerical constrained minimization
library(Rsolnp)

st.val <- solve(t(A)%*%A, t(A)%*%y)
optnum <- solnp(pars=st.val, fun=fxFUN, 
                eqfun=hxFUN, eqB=0)
xhat <- optnum$pars
cbind(x0,x,xhat)
fxFUN(xhat)
hxFUN(xhat)
points(xhat[1], xhat[2], pch=18, col=5, cex=1)






###########################################
## OLS + nonlinear constraint
## 222222222222222222222222222222222222222

rm(list = ls())
options(device="X11")
library(quadprog)

xx <- 30:100
A <- cbind(1,xx)
xUC.T <- c(-10,0.1)
y <- A%*%xUC.T + rnorm(length(xx),0,0.1)
plot(xx, y)
lines(xx, A%*%xUC.T, col=2)

fxFUN <- function(x){
  y.hat <- A%*%x
  #fx <- sum( (y-y.hat)^2 )
  fx <- t(y-y.hat)%*%(y-y.hat)
  fx <- as.numeric(fx)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx1 <- 2 * (t(A)%*%A%*%x) - 2*(t(A)%*%y)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx2 <- 2*(t(A)%*%A)
  return(fx2)
}

exp(xUC.T[1]) + exp(xUC.T[2])

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- exp(x1) + exp(x2) - 1.1 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(exp(x1), exp(x2))
  return(hx1)
}

x0 <- solve(t(A)%*%A, t(A)%*%y)
lines(xx, A%*%xUC.T, col=2, lty=2)
lines(xx, A%*%x0, col=4)

x <- c(-10.5, 0.11)
for(it in 1:20){
  fx1 <- t(A)%*%A%*%x - t(A)%*%y
  fx2 <- t(A)%*%A
  d <- solve(fx2, -fx1)
  x <- x + d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
x
x0



x <- x0
for(it in 1:20){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  
  Am <- hx1
  bv <- -hx
  # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
  #                 Amat=Am, bvec=bv, meq=1)
  # d <- opt$solution
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:2]
  
  x <- x + d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}

plot(xx, y)
lines(xx, A%*%x, col=2, lwd=2, lty=2)
lines(xx, A%*%x0, col=3, lwd=2, lty=1)



## graphical view
x1s <- seq(-11, -9, length=100)
l1 <- length(x1s)
x2s <- seq(0.09, 0.11, length=99)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}


brek <- quantile(FX, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FX, breaks=brek, col=colo)
## eq constr
x2s.eq <- seq(0.09, 0.11, length=1000)
x1s.eq <- log(1.1 - exp(x2s.eq))
lines(x1s.eq, x2s.eq, col=4, pch=16, cex=0.5)

## SQL fit
points(x[1], x[2], pch=16, col=1, cex=2)
## uncostrained fit
points(x0[1], x0[2], pch=17, col=3, cex=2)

cbind(x0,x)


## numerical constrained minimization
library(Rsolnp)

st.val <- solve(t(A)%*%A, t(A)%*%y)
optnum <- solnp(pars=st.val, fun=fxFUN, 
                eqfun=hxFUN, eqB=0)
xhat <- optnum$pars
cbind(x0,x,xhat)
fxFUN(xhat)
hxFUN(xhat)
points(xhat[1], xhat[2], pch=18, col=5, cex=1)









###########################################
## OLS + nonlinear constraint on fitted values
## 33333333333333333333333333333333333333333

rm(list = ls())
options(device="X11")
library(quadprog)

xx <- sort(runif(100))
m <- length(xx)
A <- cbind(1,xx)
xUC.T <- c(1,2)
yT <- A%*%xUC.T 
y <- yT + rnorm(length(xx),0,0.1)
plot(xx, y)
lines(xx, yT, col=2)

x0 <- solve(t(A)%*%A, t(A)%*%y)
y0 <- A%*%x0
lines(xx, yT, col=2, lty=2)
lines(xx, y0, col=4)


fxFUN <- function(x){
  y.hat <- A%*%x
  #fx <- sum( (y-y.hat)^2 )
  fx <- t(y-y.hat)%*%(y-y.hat)
  fx <- as.numeric(fx)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx1 <- 2 * (t(A)%*%A%*%x) - 2*(t(A)%*%y)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx2 <- 2*(t(A)%*%A)
  return(fx2)
}

sum(exp(A%*%xUC.T))
sum(exp(A%*%x0))
sum(exp(y))

x <- x0
hxFUN <-  function(x){
  yy <- A%*%x
  hx <- sum(exp(yy)) - 800 ## !!!!!
  return(hx)
}

hx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  exp.y <- exp(A%*%x)
  A1.exp.y <- A[,1] * exp.y
  A2.exp.y <- A[,2] * exp.y
  sum.A1.exp.y <- sum(A1.exp.y)
  sum.A2.exp.y <- sum(A2.exp.y)
  hx1 <- rbind(sum.A1.exp.y, sum.A2.exp.y)
  return(hx1)
}




x <- x0
for(it in 1:20){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)

  Am <- hx1
  bv <- -hx
  # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
  #                 Amat=Am, bvec=bv, meq=1)
  # d <- opt$solution
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:2]

  x <- x + d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
yC <- A%*%x
cbind(x,x0)
hxFUN(x);hxFUN(x0);hxFUN(xUC.T)

plot(xx, y)
lines(xx, yC, col=2, lwd=2, lty=2)
lines(xx, y0, col=3, lwd=2, lty=1)

## graphical view
x1s <- seq(0, 2, length=1000)
l1 <- length(x1s)
x2s <- seq(0.5, 3, length=999)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}


brek <- quantile(FX,
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FX, breaks=brek, col=colo)
## eq constr
x2s.eq <- seq(0.5, 3, length=1000)
x1s.eq <- numeric(length(x2s.eq))
for(j in 1:length(x2s.eq)){
  bla <- numeric(l1)
  for(i in 1:l1){
    bla[i] <- hxFUN(x=c(x1s[i], x2s.eq[j]))
  }
  x1s.eq[j] <- x1s[which.min(abs(bla))]
  # cat(hxFUN(c(x1s.eq[j], x2s.eq[j])), "\n")
}
lines(x1s.eq, x2s.eq, col=4, pch=16, cex=0.5)
## SQL fit
points(x[1], x[2], pch=16, col=1, cex=2)
## uncostrained fit
points(x0[1], x0[2], pch=17, col=3, cex=2)

cbind(x0,x)


## numerical constrained minimization
library(Rsolnp)

st.val <- solve(t(A)%*%A, t(A)%*%y)
optnum <- solnp(pars=st.val, fun=fxFUN,
                eqfun=hxFUN, eqB=0)
xhat <- optnum$pars
cbind(x0,x,xhat)
fxFUN(xhat)
hxFUN(xhat)
points(xhat[1], xhat[2], pch=18, col=5, cex=1)










## Poisson-likelihood example with two variables
## Poisson-likelihood example with two variables
## Poisson-likelihood example with two variables
## Poisson-likelihood example with two variables

rm(list = ls())
options(device="X11")
library(quadprog)

m <- 100
xx <- sort(runif(m))
A <- cbind(1, xx)
xUC.T <- c(3, 2)
etaT <- A%*%xUC.T
muT <- exp(etaT)
range(muT)
y <- rpois(m, lambda=muT)
range(y)

par(mfrow=c(1,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)))
lines(xx, muT, col=2, lwd=2)
plot(xx, log(y), t="p", pch=16)
lines(xx, etaT, col=2, lwd=2)
par(mfrow=c(1,1))


exp(xUC.T[1]) + exp(xUC.T[2])

# x0 <- lm(log(y+0.1)~xx)$coef
# x <- x0
# for(it in 1:20){
#   eta <-  A%*%x
#   mu <- exp(eta)
#   z <- (y-mu)/mu + eta
#   W <- diag(c(mu))
#   fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
#   fx2 <- t(A)%*%W%*%A
#   d <- solve(fx2, -fx1)
#   x <- x+d
#   cat(x, "\n")
#   if(max(abs(d))<1e-6) break
# }
# cbind(x0, glm(y~xx, family=poisson)$coef, x)

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- exp(x1) + exp(x2) - 30 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(exp(x1), exp(x2))
  return(hx1)
}

x0 <- as.vector(glm(y~xx, family=poisson)$coef)
eta <-  A%*%x0
mu <- exp(eta)
z <- (y-mu)/mu + eta
W <- diag(c(mu))
tAWA <- t(A)%*%W%*%A
tAWz <- t(A)%*%W%*%z
hx <- hxFUN(x0)
hx1 <- hx1FUN(x0)
## approach 1
fx1 <- tAWA%*%x0 - tAWz
fx2 <- tAWA
LHS <- rbind(cbind(fx2, hx1),
             cbind(t(hx1), 0))
RHS <- rbind(fx1, hx)
d <- solve(LHS, RHS)[1:2]
x1 <- x0-d
cbind(x0,x1)
## approach 2
bla <- matrix(hx127, 2,1)
LHS <- rbind(cbind(tAWA, bla),
             cbind(t(bla), 0))
RHS <- rbind(tAWz, hx+30)
x1 <- solve(LHS, RHS)[1:2]
cbind(x0,x1)


for(it in 1:20){
  eta <-  A%*%x
  mu <- exp(eta)
  z <- (y-mu)/mu + eta
  W <- diag(c(mu))

  # tAWA <- t(A)%*%W%*%A
  # tAWz <- t(A)%*%W%*%z
  # 
  # 
  # Am <- hx1
  # bv <- -hx
  # 
  # LHS <- rbind(cbind(tAWA, Am),
  #              cbind(t(Am), 0))
  # RHS <- rbind(tAWz, bv)
  # d <- solve(LHS, RHS)[1:2]
  # 
  
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
  fx2 <- t(A)%*%W%*%A
  LHS <- rbind(cbind(fx2, hx1),
               cbind(t(hx1), 0))
  RHS <- rbind(-fx1, -hx)
  d <- solve(LHS, RHS)[1:2]
  x <- x+d
  
  
  
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(x);hxFUN(x0);hxFUN(xUC.T)







x0 <- as.vector(glm(y~xx, family=poisson)$coef)
x <- x0
for(it in 1:20){
  eta <-  A%*%x
  mu <- exp(eta)
  z <- (y-mu)/mu + eta
  W <- diag(c(mu))
  fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
  fx2 <- t(A)%*%W%*%A
  
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:2]
  x <- x+d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(x);hxFUN(x0);hxFUN(xUC.T)



par(mfrow=c(1,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)))
lines(xx, muT, col=2, lwd=2)
lines(xx, exp(A%*%x0), col=3, lwd=2, lty=2)
lines(xx, exp(A%*%x), col=4, lwd=2, lty=1)


plot(xx, log(y), t="p", pch=16)
lines(xx, etaT, col=2, lwd=2)
lines(xx, A%*%x0, col=3, lwd=2, lty=2)
lines(xx, A%*%x, col=4, lwd=2, lty=1)

par(mfrow=c(1,1))


## graphical view
fxFUN <- function(x){
  y.hat <- exp(A%*%x)
  fx <- sum(y*log(y.hat) - y.hat)
  fx <- as.numeric(fx)
  return(fx)
}
x1s <- seq(2, 4, length=1000)
l1 <- length(x1s)
x2s <- seq(1, 3, length=999)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- -fxFUN(x=c(x1s[i], x2s[j]))
  }
}


brek <- quantile(FX, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FX, breaks=brek, col=colo)
## eq constr
x2s.eq <- seq(1, 3, length=1000)
x1s.eq <- log(30 - exp(x2s.eq))
lines(x1s.eq, x2s.eq, col=4, pch=16, lwd=3)

## SQL fit
points(x[1], x[2], pch=16, col=1, cex=2)
## uncostrained fit
points(x0[1], x0[2], pch=17, col=3, cex=2)

cbind(x0,x)










## Poisson-likelihood example with two variables and offset
## Poisson-likelihood example with two variables and offset
## Poisson-likelihood example with two variables and offset
## Poisson-likelihood example with two variables and offset

## Similar to a Gompertz
rm(list = ls())
options(device="X11")
library(quadprog)

xx <- 30:100
m <- length(xx)
A <- cbind(1,xx)
xUC.T <- c(-10,0.1)
etaT <- A%*%xUC.T
muT <- exp(etaT)
E0 <- read.table("E.txt")
e <- E0[30:100,30]

yT <- e*muT
y <- rpois(m, yT)
lmx <- log(y/e)

par(mfrow=c(1,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)))
lines(xx, yT, col=2, lwd=2)
plot(xx, lmx)
lines(xx, etaT, col=2)
par(mfrow=c(1,1))


exp(xUC.T[1]) + exp(xUC.T[2])

# x0 <- lm(log((y+1)/(e+1))~xx)$coef
# x <- x0
# for(it in 1:20){
#   eta <- A%*%x
#   mu <- e*exp(eta)
#   z <- (y-mu)/mu + eta
#   W <- diag(c(mu))
#   fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
#   fx2 <- t(A)%*%W%*%A
#   d <- solve(fx2, -fx1)
#   x <- x+d
#   cat(x, "\n")
#   if(max(abs(d))<1e-6) break
# }
# cbind(x0, glm(y~xx, offset=log(e), family=poisson)$coef, x)

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx <- exp(x1) + exp(x2) - 1.1 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  hx1 <- rbind(exp(x1), exp(x2))
  return(hx1)
}

x0 <- as.vector(glm(y~xx, offset=log(e), family=poisson)$coef)
x <- x0
for(it in 1:20){
  eta <- A%*%x
  mu <- e*exp(eta)
  z <- (y-mu)/mu + eta
  W <- diag(c(mu))
  fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
  fx2 <- t(A)%*%W%*%A
  
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:2]
  x <- x+d
  cat(x, "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(x);hxFUN(x0);hxFUN(xUC.T)

etaC <- A%*%x


t(A)%*%W%*%A%*%x
t(A)%*%W%*%etaC


par(mfrow=c(1,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)))
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(A%*%x0), col=3, lwd=2, lty=2)
lines(xx, e*exp(A%*%x), col=4, lwd=2, lty=1)

plot(xx, lmx, t="p", pch=16)
lines(xx, etaT, col=2, lwd=2)
lines(xx, A%*%x0, col=3, lwd=2, lty=2)
lines(xx, A%*%x, col=4, lwd=2, lty=1)

par(mfrow=c(1,1))



## graphical view
fxFUN <- function(x){
  y.hat <- e*exp(A%*%x)
  fx <- sum(y*log(y.hat) - y.hat)
  fx <- -as.numeric(fx)
  return(fx)
}

x1s <- seq(-11, -9, length=100)
l1 <- length(x1s)
x2s <- seq(0.09, 0.11, length=99)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    FX[i,j] <- fxFUN(x=c(x1s[i], x2s[j]))
  }
}

brek <- quantile(FX, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FX, breaks=brek, col=colo)
## eq constr
x2s.eq <- seq(0.09, 0.11, length=1000)
x1s.eq <- log(1.1 - exp(x2s.eq))
lines(x1s.eq, x2s.eq, col=4, pch=16, cex=0.5)
## SQL fit
points(x[1], x[2], pch=16, col=1, cex=2)
## uncostrained fit
points(x0[1], x0[2], pch=17, col=3, cex=2)

cbind(x0,x)




















## Poisson-likelihood example with fertility law
## Poisson-likelihood example with fertility law
## Poisson-likelihood example with fertility law
## Poisson-likelihood example with fertility law

rm(list = ls())
options(device="X11")
library(quadprog)

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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
lmx <- log(y/e)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## observed TFR
TFRobs <- sum(y/e)


## penalty stuff
Im <- diag(m)
D <- diff(Im, diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

eta <- log((y+1)/(e+1))
for(it in 1:200){
  mu <- e*exp(eta)
  W <- diag(as.vector(mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- y - mu + mu * eta
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  d <- solve(fx2, -fx1)
  eta <- eta + d
  if(max(abs(d))<1e-4) break
  cat(it, max(abs(d)), "\n")
}
## only smooth
etaS <- eta
sum(exp(etaS));TFRobs

hxFUN <-  function(eta){
  eta1 <- eta[1]
  hx <- exp(eta1)
  for(i in 2:m){
    eta1 <- eta[i]
    hx <- hx + exp(eta1)
  }
  hx <- hx - 3## !!!!!
  return(hx)
}

hx1FUN <-  function(eta){
  eta1 <- eta[1]
  hx1 <- rbind(exp(eta1))
  for(i in 2:m){
    eta1 <- eta[i]
    hx1 <- rbind(hx1, exp(eta1))
  }
  return(hx1)
}

eta <- etaS
for(it in 1:20){
  mu <- e*exp(eta)
  W <- diag(as.vector(mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- y - mu + mu * eta
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  
  hx <- hxFUN(eta)
  hx1 <- hx1FUN(eta)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:m]
  eta <- eta+d
  cat(max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(eta);hxFUN(etaS)

TFRsmo <- sum(exp(etaS))
TFRcon <- sum(exp(eta))
cbind(TFRobs, TFRsmo, TFRcon)

par(mfrow=c(2,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(etaS), col=3, lwd=2, lty=2)
lines(xx, e*exp(eta), col=4, lwd=2, lty=1)

plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, yT/e, col=2, lwd=2)
lines(xx, exp(etaS), col=3, lwd=2, lty=2)
lines(xx, exp(eta), col=4, lwd=2, lty=1)

plot(xx, log(y/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))










## Poisson-likelihood example with fertility law
## with B-splines and TFR constraint
## Poisson-likelihood example with fertility law
## with B-splines and TFR constraint
## Poisson-likelihood example with fertility law
## with B-splines and TFR constraint
## Poisson-likelihood example with fertility law
## with B-splines and TFR constraint

rm(list = ls())
options(device="X11")
library(quadprog)
library(MortalitySmooth)
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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
lmx <- log(y/e)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## observed TFR
TFRobs <- sum(y/e)

B <- MortSmooth_bbase(x=xx, xl=min(xx), 
                      xr=max(xx),
                      ndx=floor(m/4), deg=3)
nb <- ncol(B)


## penalty stuff
D <- diff(diag(nb), diff=3)
tDD <- t(D)%*%D
lambda <- 10^2
P <- lambda*tDD

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
sum(exp(etaS));TFRobs

hxFUN <-  function(betas){
  eta <- B%*%betas
  exp.eta <- exp(eta)
  sum.exp.eta <- sum(exp.eta)
  hx <- sum.exp.eta - 4 # !!!
  return(hx)
}

betas <- betasS
hx1FUN <- function(betas){
  exp.eta <- exp(B%*%betas)
  hx1 <- sum(B[,1]*exp.eta)
  for(i in 2:nb){
    hx1 <- rbind(hx1, sum(B[,i]*exp.eta))
  }
  return(hx1)
}


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
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
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

etaC <- B%*%betas

par(mfrow=c(1,2))
plot(tBWBpP%*%betas)
plot(e*exp(etaC))


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







## Poisson-likelihood example with fertility law
## with B-splines and MAB and TFR constraints
## Poisson-likelihood example with fertility law
## with B-splines and MAB and TFR constraints
## Poisson-likelihood example with fertility law
## with B-splines and MAB and TFR constraints
## Poisson-likelihood example with fertility law
## with B-splines and MAB and TFR constraints

rm(list = ls())
options(device="X11")
library(quadprog)
library(MortalitySmooth)
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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
mx <- y/e
lmx <- log(mx)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## observed TFR
TFRobs <- sum(mx)
## observed MAB
MABobs <- sum(mx* (xx+0.5) ) / sum(mx)

## B-splines
B <- MortSmooth_bbase(x=xx, xl=min(xx), 
                      xr=max(xx),
                      ndx=floor(m/4), deg=3)
nb <- ncol(B)


## penalty stuff
D <- diff(diag(nb), diff=3)
tDD <- t(D)%*%D
lambda <- 10^2
P <- lambda*tDD

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
muS <- exp(etaS)

TFRsmo <- sum(muS)
cbind(TFRsmo, TFRobs)

MABsmo <- sum(muS* (xx+0.5) ) / sum(muS)
cbind(MABsmo, MABobs)


## constraints of TFR
TFRtar <- 4
hxFUN <-  function(betas){
  eta <- B%*%betas
  exp.eta <- exp(eta)
  sum.exp.eta <- sum(exp.eta)
  hx <- sum.exp.eta - TFRtar
  return(hx)
}
hx1FUN <- function(betas){
  exp.eta <- exp(B%*%betas)
  hx1 <- sum(B[,1]*exp.eta)
  for(i in 2:nb){
    hx1 <- rbind(hx1, sum(B[,i]*exp.eta))
  }
  return(hx1)
}


## constraints of MAB
MABtar <- 28
gxFUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  num <- sum(exp(eta)*xxx)
  den <- sum(exp(eta))
  gx <- num/den - MABtar ## !!!!!!!!!!!
  return(gx)
}



gx1FUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  mu <- exp(eta)
  mux <- mu*xxx
  sum.mux <- sum(mux)
  den <- sum(mu)
  den2 <- den^2
  
  gx1 <- matrix(0, nb, 1)
  for(i in 1:nb){
    Bmux.i <- B[,i]*mux
    part1 <- sum(Bmux.i)/den
    Bmu <- B[,i]*mu
    part2 <- ( sum.mux * sum(Bmu) ) / den2
    gx1[i,1] <- part1 - part2
  }
  return(gx1)
}



betas <- betasS
for(it in 1:100){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  gx <- gxFUN(betas)
  gx1 <- gx1FUN(betas)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(2)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
}
betasC <- betas
etaC <- eta
muC <- exp(etaC)

TFRcon <- sum(muC)
cbind(TFRsmo, TFRobs, TFRtar, TFRcon)

MABcon <- sum(muC* (xx+0.5) ) / sum(muC)
cbind(MABsmo, MABobs, MABtar, MABcon)


plot(xx, mx, t="h", lwd=2, 
     ylim=range(0, mx, muS, muC), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, muS, col=3, lwd=2, lty=2)
lines(xx, muC, col=4, lwd=2, lty=1)



par(mfrow=c(2,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(etaS), col=3, lwd=2, lty=2)
lines(xx, e*exp(eta), col=4, lwd=2, lty=1)

plot(xx, mx, t="h", lwd=2, 
     ylim=range(0, mx, muS, muC), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, muS, col=3, lwd=2, lty=2)
lines(xx, muC, col=4, lwd=2, lty=1)

plot(xx, log(y/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))





## Poisson-likelihood example with fertility law
## with B-splines and TFR, MAB and VAR constraints
## Poisson-likelihood example with fertility law
## with B-splines and TFR, MAB and VAR constraints
## Poisson-likelihood example with fertility law
## with B-splines and TFR, MAB and VAR constraints
## Poisson-likelihood example with fertility law
## with B-splines and TFR, MAB and VAR constraints


rm(list = ls())
options(device="X11")
library(quadprog)
library(MortalitySmooth)
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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
mx <- y/e
lmx <- log(mx)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## observed TFR
TFRobs <- sum(mx)
## observed MAB
MABobs <- sum(mx* (xx+0.5) ) / sum(mx)
## observed VAB
VABobs <- sum((((xx+0.5) - MABobs)^2)*mx)/sum(mx)

## B-splines
B <- MortSmooth_bbase(x=xx, xl=min(xx), 
                      xr=max(xx),
                      ndx=floor(m/4), deg=3)
nb <- ncol(B)


## penalty stuff
D <- diff(diag(nb), diff=3)
tDD <- t(D)%*%D
lambda <- 10^2
P <- lambda*tDD

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
muS <- exp(etaS)

TFRsmo <- sum(muS)
cbind(TFRsmo, TFRobs)

MABsmo <- sum(muS* (xx+0.5) ) / sum(muS)
cbind(MABsmo, MABobs)

VABsmo <- sum((((xx+0.5) - MABsmo)^2)*muS)/sum(muS)
cbind(VABsmo, VABobs)


## constraints of TFR
TFRtar <- 4
hxFUN <-  function(betas){
  eta <- B%*%betas
  exp.eta <- exp(eta)
  sum.exp.eta <- sum(exp.eta)
  hx <- sum.exp.eta - TFRtar
  return(hx)
}
hx1FUN <- function(betas){
  exp.eta <- exp(B%*%betas)
  hx1 <- sum(B[,1]*exp.eta)
  for(i in 2:nb){
    hx1 <- rbind(hx1, sum(B[,i]*exp.eta))
  }
  return(hx1)
}


## constraints of MAB
MABtar <- 28
gxFUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  num <- sum(exp(eta)*xxx)
  den <- sum(exp(eta))
  gx <- num/den - MABtar 
  return(gx)
}
gx1FUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  mu <- exp(eta)
  mux <- mu*xxx
  sum.mux <- sum(mux)
  den <- sum(mu)
  den2 <- den^2
  
  gx1 <- matrix(0, nb, 1)
  for(i in 1:nb){
    Bmux.i <- B[,i]*mux
    part1 <- sum(Bmux.i)/den
    Bmu <- B[,i]*mu
    part2 <- ( sum.mux * sum(Bmu) ) / den2
    gx1[i,1] <- part1 - part2
  }
  return(gx1)
}

## constraints of VAB
VABtar <- 40
uxFUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  mu <- exp(eta)
  mux <- mu*xxx
  tfr <- sum(mu)
  mab <- sum(mux)/tfr
  vab <- sum(((xxx - mab)^2)*mu) / tfr
  ux <- vab - VABtar
  return(ux)
}
ux1FUN <- function(betas){
  xxx <- xx+0.5
  eta <- B%*%betas
  mu <- exp(eta)
  mux <- mu*xxx
  tfr <- sum(mu)
  tfr2 <- tfr^2
  mab <- sum(mux)/tfr
  x.mab <- xxx-mab
  x.mab2 <- x.mab^2
  x.mab2.mu <- x.mab2*mu
  sum.x.mab2.mu <- sum(x.mab2.mu)
  ux1 <- matrix(0, nb, 1)
  for(i in 1:nb){
    part1 <- sum(B[,i]*x.mab2.mu)/tfr
    part2 <- (sum(B[,i]*mu) * sum.x.mab2.mu)/tfr2
    ux1[i,1] <- part1 - part2
  }
  return(ux1)
}
  




betas <- betasS
for(it in 1:100){
  eta <- B%*%betas
  mu <- e*exp(eta)
  z <- (y - mu)/mu + eta
  w <- c(mu)
  tBWB <- t(B) %*% (w * B)
  tBWBpP <- tBWB + P
  tBWz <- t(B) %*% (w * z)
  
  fx1 <- tBWBpP%*%betas - tBWz
  fx2 <- tBWBpP
  
  hx <- hxFUN(betas)
  hx1 <- hx1FUN(betas)
  gx <- gxFUN(betas)
  gx1 <- gx1FUN(betas)
  ux <- uxFUN(betas)
  ux1 <- ux1FUN(betas)
  
  Am <- cbind(hx1, gx1, ux1)
  bv <- -c(hx, gx, ux)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0*diag(3)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:nb]
  betas <- betas + d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-4) break
}
betasC <- betas
etaC <- eta
muC <- exp(etaC)

TFRcon <- sum(muC)
cbind(TFRsmo, TFRobs, TFRtar, TFRcon)

MABcon <- sum(muC* (xx+0.5) ) / sum(muC)
cbind(MABsmo, MABobs, MABtar, MABcon)

VABcon <- sum((((xx+0.5) - MABcon)^2)*muC)/sum(muC)
cbind(VABsmo, VABobs, VABtar, VABcon)

plot(xx, mx, t="h", lwd=2, 
     ylim=range(0, mx, muS, muC), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, muS, col=3, lwd=2, lty=2)
lines(xx, muC, col=4, lwd=2, lty=1)



par(mfrow=c(2,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(etaS), col=3, lwd=2, lty=2)
lines(xx, e*exp(eta), col=4, lwd=2, lty=1)

plot(xx, mx, t="h", lwd=2, 
     ylim=range(0, mx, muS, muC), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, muS, col=3, lwd=2, lty=2)
lines(xx, muC, col=4, lwd=2, lty=1)

plot(xx, log(y/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))









## Poisson-likelihood example with fertility law 
## and extrapolation
## Poisson-likelihood example with fertility law 
## and extrapolation
## Poisson-likelihood example with fertility law 
## and extrapolation
## Poisson-likelihood example with fertility law 
## and extrapolation

rm(list = ls())
options(device="X11")
library(quadprog)

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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
lmx <- log(y/e)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, y/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## assume uncompleted fertility pattern 
## with known/assumed total fertility rate  
## after a certain age we know nothing
## but the TFR

TFRtrue <- sum(exp(etaT))

wxx <- 35
whi <- xx>=wxx
y[whi] <- 999
e[whi] <- 9999
wei <- rep(1, m)
wei[whi] <- 0

## penalty stuff
Im <- diag(m)
D <- diff(Im, diff=2)
tDD <- t(D)%*%D
lambda <- 10^2
P <- lambda*tDD

eta <- log((y+1)/(e+1))
for(it in 1:200){
  mu <- e*exp(eta)
  W <- diag(as.vector(wei*mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- wei*(y - mu + mu * eta)
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  d <- solve(fx2, -fx1)
  eta <- eta + d
  if(max(abs(d))<1e-4) break
  cat(it, max(abs(d)), "\n")
}
## only smooth
etaS <- eta
TFRsmo <- sum(exp(etaS))
TFRsmo;TFRtrue


hxFUN <-  function(eta){
  eta1 <- eta[1]
  hx <- exp(eta1)
  for(i in 2:m){
    eta1 <- eta[i]
    hx <- hx + exp(eta1)
  }
  hx <- hx - TFRtrue## !!!!!
  return(hx)
}

hx1FUN <-  function(eta){
  eta1 <- eta[1]
  hx1 <- rbind(exp(eta1))
  for(i in 2:m){
    eta1 <- eta[i]
    hx1 <- rbind(hx1, exp(eta1))
  }
  return(hx1)
}

eta <- etaS
for(it in 1:20){
  mu <- e*exp(eta)
  W <- diag(as.vector(wei*mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- wei*(y - mu + mu * eta)
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  
  hx <- hxFUN(eta)
  hx1 <- hx1FUN(eta)
  Am <- hx1
  bv <- -hx
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), 0))
  RHS <- rbind(-fx1, bv)
  d <- solve(LHS, RHS)[1:m]
  eta <- eta+d
  cat(max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(eta);hxFUN(etaS)

TFRcon <- sum(exp(eta))
cbind(TFRsmo, TFRcon,TFRtrue)

par(mfrow=c(1,3))
yNA <- y
yNA[whi] <- NA

plot(xx, yNA/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, exp(etaS), col=3, lwd=2, lty=2)
lines(xx, exp(eta), col=4, lwd=2, lty=1)

plot(xx, log(yNA/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))








###########################################
## OLS + TWO nonlinear constraint
## OLS + TWO nonlinear constraint
## OLS + TWO nonlinear constraint
## OLS + TWO nonlinear constraint
## ########################################

rm(list = ls())
options(device="X11")
library(quadprog)
library(rgl)
xx1 <- 30:100
m1 <- length(xx1)
xx2 <- 30:99
m2 <- length(xx2)

A <- as.matrix(expand.grid(1, xx1, xx2))
xUC.T <- c(-10,-1,0.5)
yT <- A%*%xUC.T
y <- yT + rnorm(nrow(A),0,5)
YT <- matrix(yT, m1, m2)
Y <- matrix(y, m1, m2)

plot3d(A[,2], A[,3], y, type="p", col="red", 
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, YT, 
          back = 'line', front = 'line', col = 'orange', 
          lwd = 2.5, alpha = 0.5)


p <- persp(xx1, xx2, YT, theta=30, phi=30, 
             zlim=range(y,yT), 
             col="lightblue",expand = 0.8,shade = 0.2,
             xlab="x1", ylab="x2", zlab="y")
obs <- trans3d(A[,2], A[,3], y, p)
points(obs, col="red",pch=16, cex=0.5)


fxFUN <- function(x){
  y.hat <- A%*%x
  #fx <- sum( (y-y.hat)^2 )
  fx <- t(y-y.hat)%*%(y-y.hat)
  fx <- as.numeric(fx)
  return(fx)
}

fx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx1 <- 2 * (t(A)%*%A%*%x) - 2*(t(A)%*%y)
  return(fx1)
}

fx2FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  fx2 <- 2*(t(A)%*%A)
  return(fx2)
}

x0 <- solve(t(A)%*%A, t(A)%*%y)
Ylm <- matrix(A%*%x0, m1, m2)

# p <- persp3d(xx1, xx2, Ylm, theta=30, phi=30, 
#              zlim=range(y,yT), 
#              col="lightblue",expand = 0.5,shade = 0.2,
#              xlab="x1", ylab="x2", zlab="y")
# points3d(A[,2], A[,3], y, col=2, pch=16, cex=0.5)
# 
# 
# x <- c(-10.5, -1.2, 0.7)
# for(it in 1:20){
#   fx1 <- t(A)%*%A%*%x - t(A)%*%y
#   fx2 <- t(A)%*%A
#   d <- solve(fx2, -fx1)
#   x <- x + d
#   cat(x, "\n")
#   if(max(abs(d))<1e-6) break
# }
# x
# x0
# 



exp(xUC.T[1]) + exp(xUC.T[2]) + exp(xUC.T[3])

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  hx <- exp(x1) + exp(x2) + exp(x3) - 3 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  hx1 <- rbind(exp(x1), exp(x2), exp(x3))
  return(hx1)
}

# 
# 
# x <- x0
# for(it in 1:20){
#   fx <- fxFUN(x)
#   fx1 <- fx1FUN(x)
#   fx2 <- fx2FUN(x)
#   hx <- hxFUN(x)
#   hx1 <- hx1FUN(x)
#   
#   Am <- hx1
#   bv <- -hx
#   # opt <- solve.QP(Dmat=fx2, dvec=-fx1,
#   #                 Amat=Am, bvec=bv, meq=1)
#   # d <- opt$solution
#   LHS <- rbind(cbind(fx2, Am),
#                cbind(t(Am), 0))
#   RHS <- rbind(-fx1, bv)
#   d <- solve(LHS, RHS)[1:3]
#   
#   x <- x + d
#   cat(x, "\n")
#   if(max(abs(d))<1e-6) break
# }
# hxFUN(x)
# hxFUN(x0)
# 

gxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  num <- exp(x1)*1 + exp(x2)*2 + exp(x3)*3
  den <- exp(x1)+exp(x2)+exp(x3)
  gx <- num/den - 2.9 ## !!!!!!!!!!!
  return(gx)
}

numT <- exp(xUC.T[1])*1 + exp(xUC.T[2])*2 + exp(xUC.T[3])*3
denT <- exp(xUC.T[1]) + exp(xUC.T[2]) + exp(xUC.T[3])
numT/denT
gxFUN(xUC.T)
gxFUN(x0)

hxFUN(xUC.T)
hxFUN(x0)


gx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  num1 <- -(2*exp(x3) +   exp(x2)) * exp(x1)
  num2 <- -(  exp(x3) -   exp(x1)) * exp(x2)
  num3 <-  (  exp(x2) + 2*exp(x1)) * exp(x3)
  den <- (exp(x1)+exp(x2)+exp(x3))^2
  gx1.1 <- num1/den
  gx1.2 <- num2/den
  gx1.3 <- num3/den
  gx1 <- rbind(gx1.1, gx1.2, gx1.3)
  return(gx1)
}

x <- x0
for(it in 1:20){
  fx <- fxFUN(x)
  fx1 <- fx1FUN(x)
  fx2 <- fx2FUN(x)
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  gx <- gxFUN(x)
  gx1 <- gx1FUN(x)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  opt <- solve.QP(Dmat=fx2, dvec=-fx1,
                  Amat=Am, bvec=bv, meq=2)
  d <- opt$solution

  # LHS <- rbind(cbind(fx2, Am),
  #              cbind(t(Am), matrix(0,2,2)))
  # RHS <- c(-fx1, bv)
  # d <- solve(LHS, RHS)[1:3]
  # 
  x <- x + d
  cat(it, x, "\n")
  if(max(abs(d))<1e-6) break
}
gxFUN(x)
gxFUN(x0)
hxFUN(x)
hxFUN(x0)

exp(x0[1]) + exp(x0[2]) + exp(x0[3])
exp(x[1]) + exp(x[2]) + exp(x[3])

cbind(x,x0)

Yclm <- matrix(A%*%x, m1, m2)

par(mfrow=c(2,2))
whi <- floor(seq(1,m2,length=4))
for(i in 1:4){
  plot(xx1, Y[,whi[i]])
  title(main=paste("xx2=", xx2[whi[i]]))
  lines(xx1, Yclm[,whi[i]], col=4, lwd=2)
  lines(xx1, Ylm[,whi[i]], col=2, lty=2)
}
par(mfrow=c(2,2))
whi <- floor(seq(1,m1,length=4))
for(i in 1:4){
  plot(xx2, Y[whi[i],])
  title(main=paste("xx1=", xx1[whi[i]]))
  lines(xx2, Yclm[whi[i],], col=4, lwd=2)
  lines(xx2, Ylm[whi[i],], col=2, lty=2)
}
par(mfrow=c(1,1))

plot3d(A[,2], A[,3], y, type="p", col="red",
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, Ylm,
          back='line', front='line', col='orange',
          lwd=1.5, alpha=0.3)
surface3d(xx1, xx2, Yclm,
          back='line', front='line', col='cyan',
          lwd=1.5, alpha=0.5)

















## Poisson-likelihood example with THREE variables and two constraints
## Poisson-likelihood example with THREE variables and two constraints
## Poisson-likelihood example with THREE variables and two constraints
## Poisson-likelihood example with THREE variables and two constraints

rm(list = ls())
options(device="X11")
library(quadprog)
library(rgl)

xx1 <- seq(0,1,length=70)#sort(runif(70))
m1 <- length(xx1)
xx2 <- seq(0,1,length=69)#sort(runif(69))
m2 <- length(xx2)

A <- as.matrix(expand.grid(1, xx1, xx2))
xUC.T <- c(3,2,1)
etaT <- A%*%xUC.T
ETAT <- matrix(etaT, m1, m2)
muT <- exp(etaT)
range(muT)
y <- rpois(m1*m2, lambda=muT)
range(y)
MUT <- matrix(muT, m1, m2)
Y <- matrix(y, m1, m2)

## counts and expected values
plot3d(A[,2], A[,3], y, type="p", col="red", 
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, MUT, 
          back = 'line', front = 'line', col = 'orange', 
          lwd = 2.5, alpha = 0.5)

## log-counts and linear predictor
plot3d(A[,2], A[,3], log(y), type="p", col="red", 
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, ETAT, 
          back = 'line', front = 'line', col = 'orange', 
          lwd = 2.5, alpha = 0.5)

exp(xUC.T[1]) + exp(xUC.T[2]) + exp(xUC.T[3])

hxFUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  hx <- exp(x1) + exp(x2) + exp(x3) - 32 ## !!!!!
  return(hx)
}

hx1FUN <-  function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  hx1 <- rbind(exp(x1), exp(x2), exp(x3))
  return(hx1)
}

gxFUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  num <- exp(x1)*1 + exp(x2)*2 + exp(x3)*3
  den <- exp(x1)+exp(x2)+exp(x3)
  gx <- num/den - 1.7 ## !!!!!!!!!!!
  return(gx)
}

numT <- exp(xUC.T[1])*1 + exp(xUC.T[2])*2 + exp(xUC.T[3])*3
denT <- exp(xUC.T[1]) + exp(xUC.T[2]) + exp(xUC.T[3])
numT/denT
gxFUN(xUC.T)
gxFUN(x0)

hxFUN(xUC.T)
hxFUN(x0)


gx1FUN <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  num1 <- -(2*exp(x3) +   exp(x2)) * exp(x1)
  num2 <- -(  exp(x3) -   exp(x1)) * exp(x2)
  num3 <-  (  exp(x2) + 2*exp(x1)) * exp(x3)
  den <- (exp(x1)+exp(x2)+exp(x3))^2
  gx1.1 <- num1/den
  gx1.2 <- num2/den
  gx1.3 <- num3/den
  gx1 <- rbind(gx1.1, gx1.2, gx1.3)
  return(gx1)
}



x0 <- as.vector(glm(y~A[,2]+A[,3], family=poisson)$coef)
etaglm <- A%*%x0
muglm <- exp(etaglm)
ETAglm <- matrix(etaglm, m1, m2)
MUglm <- exp(ETAglm)

x <- x0
for(it in 1:20){
  eta <-  A%*%x
  mu <- exp(eta)
  z <- (y-mu)/mu + eta
  W <- diag(c(mu))
  fx1 <- t(A)%*%W%*%A%*%x - t(A)%*%W%*%z
  fx2 <- t(A)%*%W%*%A
  
  hx <- hxFUN(x)
  hx1 <- hx1FUN(x)
  gx <- gxFUN(x)
  gx1 <- gx1FUN(x)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), matrix(0,2,2)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:3]
  x <- x+d
  cat(it, x, "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(x);hxFUN(x0);hxFUN(xUC.T)
gxFUN(x);gxFUN(x0);gxFUN(xUC.T)

cbind(x,x0)

etacglm <- A%*%x
mucglm <- exp(etacglm)
ETAcglm <- matrix(etacglm, m1, m2)
MUcglm <- exp(ETAcglm)


## linear predictor
par(mfrow=c(2,2))
whi <- floor(seq(1,m2,length=4))
for(i in 1:4){
  plot(xx1, log(Y[,whi[i]]))
  title(main=paste("xx2=", xx2[whi[i]]))
  lines(xx1, ETAcglm[,whi[i]], col=4, lwd=2)
  lines(xx1, ETAglm[,whi[i]], col=2, lty=2)
}
par(mfrow=c(2,2))
whi <- floor(seq(1,m1,length=4))
for(i in 1:4){
  plot(xx2, log(Y[whi[i],]))
  title(main=paste("xx1=", xx1[whi[i]]))
  lines(xx2, ETAcglm[whi[i],], col=4, lwd=2)
  lines(xx2, ETAglm[whi[i],], col=2, lty=2)
}
par(mfrow=c(1,1))

## expected values
plot3d(A[,2], A[,3], y, type="p", col="red",
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, MUglm,
          back='line', front='line', col='orange',
          lwd=1.5, alpha=0.5)
surface3d(xx1, xx2, MUcglm,
          back='line', front='line', col='cyan',
          lwd=1.5, alpha=0.5)

## linear predictor
plot3d(A[,2], A[,3], log(y), type="p", col="red",
       xlab="xx1", ylab="xx2", zlab="y", site=5, lwd=15)
surface3d(xx1, xx2, ETAglm,
          back='line', front='line', col='orange',
          lwd=1.5, alpha=0.5)
surface3d(xx1, xx2, ETAcglm,
          back='line', front='line', col='blue',
          lwd=1.5, alpha=0.5)









## Poisson-likelihood example with fertility law
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law
## and constraints of TFR and mean age at childbearing (MAB)

rm(list = ls())
options(device="X11")
library(quadprog)

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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
mx <- y/e
lmx <- log(mx)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, mx, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## observed TFR
TFRobs <- sum(mx)
## observed MAB
MABobs <- sum(mx* (xx+0.5) ) / sum(mx)

## penalty stuff
Im <- diag(m)
D <- diff(Im, diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

eta <- log((y+1)/(e+1))
for(it in 1:200){
  mu <- e*exp(eta)
  W <- diag(as.vector(mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- y - mu + mu * eta
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  d <- solve(fx2, -fx1)
  eta <- eta + d
  if(max(abs(d))<1e-4) break
  cat(it, max(abs(d)), "\n")
}
## only smooth
etaS <- eta
TFRsmo <- sum(exp(etaS))
cbind(TFRobs, TFRsmo)

MABsmo <- sum(exp(etaS)* (xx+0.5) ) / sum(exp(etaS))
cbind(MABobs, MABsmo)

## constraints on TFR
hxFUN <-  function(eta){
  hx <- sum(exp(eta)) - 3 # !!!
  return(hx)
}

hx1FUN <-  function(eta){
  hx1 <- matrix(exp(eta), ncol=1)
  return(hx1)
}

## constraints of MAB
gxFUN <- function(eta){
  xxx <- xx+0.5
  num <- sum(exp(eta)*xxx)
  den <- sum(exp(eta))
  gx <- num/den - 29 ## !!!!!!!!!!!
  return(gx)
}

gx1FUN <- function(eta){
  xxx <- xx+0.5
  
  den <- sum(exp(eta))^2
  
  gx1 <- matrix(0, length(xxx), 1)
  for(i in 1:nrow(gx1)){
    num1 <- xxx-xxx[i]
    num2 <- sum(num1*exp(eta))
    num <- - num2*exp(eta[i])
    gx1[i,1] <- num/den
  }
  return(gx1)
}

eta <- etaS
for(it in 1:20){
  mu <- e*exp(eta)
  W <- diag(as.vector(mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- y - mu + mu * eta
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  
  hx <- hxFUN(eta)
  hx1 <- hx1FUN(eta)
  gx <- gxFUN(eta)
  gx1 <- gx1FUN(eta)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), matrix(0,2,2)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:m]
  eta <- eta+d
  cat(max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(eta);hxFUN(etaS)
gxFUN(eta);gxFUN(etaS)

TFRcon <- sum(exp(eta))
cbind(TFRobs, TFRsmo, TFRcon)
MABcon <- sum(exp(eta)* (xx+0.5) ) / sum(exp(eta))
cbind(MABobs, MABsmo, MABcon)


par(mfrow=c(2,2))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
lines(xx, e*exp(etaS), col=3, lwd=2, lty=2)
lines(xx, e*exp(eta), col=4, lwd=2, lty=1)

plot(xx, mx, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, exp(etaS), col=3, lwd=2, lty=2)
lines(xx, exp(eta), col=4, lwd=2, lty=1)

plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))






## Poisson-likelihood example with fertility law 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)

rm(list = ls())
options(device="X11")
library(quadprog)

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
aT <- 1.75
bT <- 3
cT <- 28

muT <- Hadwiger(aT, bT, cT, x=xx)
etaT <- log(muT)

e <- c(421660,415695,420043,415994,418983,413376,414893,414386,
       416835,413434,416274,415330,421366,429169,432178,427665,
       419672,361242,309296,308672,294653,274481,252437,290595,
       294925,298035,302802,305092,314129,318449,324701,332193,
       338598,329218,327132,328347,326924,327279,327114,321125,
       322818,329727,337702,285953)
e <- e/20

yT <- e*muT
y <- rpois(m, yT)
mx <- y/e
lmx <- log(mx)

par(mfrow=c(1,3))
plot(xx, y, t="h", lwd=2, ylim=c(0, max(y)), main="counts")
lines(xx, yT, col=2, lwd=2)
plot(xx, mx, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
plot(xx, lmx, main="log-rates")
lines(xx, etaT, col=2)
par(mfrow=c(1,1))

## assume uncompleted fertility pattern 
## with known/assumed total fertility rate  
## after a certain age we know nothing
## but the TFR and MAB

## observed TFR
TFRtrue <- sum(mx)
## observed MAB
MABtrue <- sum(mx* (xx+0.5) ) / sum(mx)

wxx <- 35
whi <- xx>=wxx
y[whi] <- 999
e[whi] <- 9999
wei <- rep(1, m)
wei[whi] <- 0

## penalty stuff
Im <- diag(m)
D <- diff(Im, diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

eta <- log((y+1)/(e+1))
for(it in 1:200){
  mu <- e*exp(eta)
  W <- diag(as.vector(wei*mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- wei*(y - mu + mu * eta)
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  d <- solve(fx2, -fx1)
  eta <- eta + d
  if(max(abs(d))<1e-4) break
  cat(it, max(abs(d)), "\n")
}
## only smooth
etaS <- eta
TFRsmo <- sum(exp(etaS))
cbind(TFRsmo,TFRtrue)

MABsmo <- sum(exp(etaS)* (xx+0.5) ) / sum(exp(etaS))
cbind(MABtrue, MABsmo)


## constraints on TFR
hxFUN <-  function(eta){
  hx <- sum(exp(eta)) - TFRtrue # !!!
  return(hx)
}

hx1FUN <-  function(eta){
  hx1 <- matrix(exp(eta), ncol=1)
  return(hx1)
}

## constraints of MAB
gxFUN <- function(eta){
  xxx <- xx+0.5
  num <- sum(exp(eta)*xxx)
  den <- sum(exp(eta))
  gx <- num/den - MABtrue ## !!!!!!!!!!!
  return(gx)
}

gx1FUN <- function(eta){
  xxx <- xx+0.5
  den <- sum(exp(eta))^2
  gx1 <- matrix(0, length(xxx), 1)
  for(i in 1:nrow(gx1)){
    num1 <- xxx-xxx[i]
    num2 <- sum(num1*exp(eta))
    num <- - num2*exp(eta[i])
    gx1[i,1] <- num/den
  }
  return(gx1)
}

eta <- etaS
for(it in 1:100){
  mu <- e*exp(eta)
  W <- diag(as.vector(wei*mu))
  tXWX <- W
  tXWXpP <- W + P
  tXWz <- wei*(y - mu + mu * eta)
  fx1 <- tXWXpP%*%eta - tXWz
  fx2 <- tXWXpP
  
  hx <- hxFUN(eta)
  hx1 <- hx1FUN(eta)
  gx <- gxFUN(eta)
  gx1 <- gx1FUN(eta)
  
  Am <- cbind(hx1, gx1)
  bv <- -c(hx, gx)
  
  LHS <- rbind(cbind(fx2, Am),
               cbind(t(Am), matrix(0,2,2)))
  RHS <- c(-fx1, bv)
  d <- solve(LHS, RHS)[1:m]
  eta <- eta+d
  cat(it, max(abs(d)), "\n")
  if(max(abs(d))<1e-6) break
}
hxFUN(eta);hxFUN(etaS)
gxFUN(eta);gxFUN(etaS)

TFRcon <- sum(exp(eta))
cbind(TFRsmo, TFRcon, TFRtrue)

MABcon <- sum(exp(eta)* (xx+0.5) ) / sum(exp(eta))
cbind(MABsmo, MABcon, MABtrue)


par(mfrow=c(1,3))
yNA <- y
yNA[whi] <- NA

plot(xx, yNA/e, t="h", lwd=2, ylim=c(0, max(y/e)), main="rates")
lines(xx, muT, col=2, lwd=2)
lines(xx, exp(etaS), col=3, lwd=2, lty=2)
lines(xx, exp(eta), col=4, lwd=2, lty=1)

plot(xx, log(yNA/e), main="log-rates")
lines(xx, etaT, col=2, lwd=2)
lines(xx, etaS, col=3, lwd=2, lty=2)
lines(xx, eta, col=4, lwd=2, lty=1)

plot(xx, eta-etaS, main="smooth - constraints")
par(mfrow=c(1,1))





## Poisson-likelihood example with fertility law 
## in 2D 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## in 2D 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## in 2D 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)
## Poisson-likelihood example with fertility law 
## in 2D 
## and extrapolation
## and constraints of TFR and mean age at childbearing (MAB)



# see ConstrainedFerility.R















































obj.FUN <- function(d, fx1, fx2, hx1, gx1){
  bla <- t(fx1)%*%d + 0.5 * t(d)%*%fx2%*%d
  return(bla)
}
eq.FUN <- function(d, fx1, fx2, hx1, gx1){
  A <- sum(hx1*d)
  return(A)
}
eq.B <- -hx
ineq.FUN <- function(d, fx1, fx2, hx1, gx1){
  A <- sum(gx1*d)
  return(A)
}
ineq.UB <- -gx

opt <- solnp(c(-1,0.5), obj.FUN, 
             eqfun=eq.FUN, eqB=eq.B,
             ineqfun=ineq.FUN, 
             ineqUB=ineq.UB, ineqLB=-100,
             fx1=fx1, fx2=fx2, hx1=hx1,
             gx1=gx1)
d <- opt$pars
eq.FUN(d=d, fx1=fx1, fx2=fx2, hx1=hx1, gx1=gx1)
eq.B
ineq.FUN(d=d, fx1=fx1, fx2=fx2, hx1=hx1, gx1=gx1)
ineq.UB


opt$pars
optim(c(0,0), blaFUN, fx1=fx1, fx2=fx2)


fx2%*%d

eigen(fx2FUN(x=c(1.07921, 1.4604)))$values

eigen(fx2)$values

A%*%xx

d[2] <- -d[1]/2

t(fx1)%*%d + 0.5 * t(d)%*%fx2%*%d

A%*%d


Dmat <- matrix(c(6.4595, -4.4044, -4.4044, 4.1579), 2, 2)
dvec <- c(1.78475, -2.17750)
A <- matrix(c(1.4604, 1.07921, -1, -1), 2, 2)
a <- c(0.42393, 1.53961)
solve.QP(Dmat=Dmat, dvec=-dvec, Amat=A, bvec=a, meq=1)



d <- c(-0.92079, 0.4604)
cbind(A%*%d, a)

1*d[1] + 2*d[2]

-1*d[1] + -1*d[2]

A%*%opt$solution


solve(fx2, fx1)

A%*%x
a


LHS0 <- cbind(fx2, t(A))
LHS1 <- cbind(A, matrix(0,2,2))
LHS <- rbind(LHS0, LHS1)
RHS <- c(-fx1, a)
coef <- solve(LHS, RHS)
xx <- coef[1:2]


  


  
  



## Example from W. Gander et al., 
## Scientific Computing - An Introduction using Maple and MATLAB,
## Springer International Publishing Switzerland 2014

rm(list = ls())
options(device="X11")


x <- c(0.288175,0.525650,1.000600,1.475550,1.713025,1.950500,2.187975,
       2.425450,2.662925,2.900400,3.137875,3.375350,3.612825,3.850300,
       4.087775,4.325250,4.562725,4.800200,5.037675,5.275150,5.512625,
       5.750100,5.845090,5.940080,6.035070,6.130060,6.225050,6.320040,
       6.415030,6.510020,6.605010,6.700000,6.794990,6.889980,6.984970,
       7.079960,7.174950,7.269940,7.364930,7.459920,7.554910,7.649900,
       7.744890,7.839880,7.934870,8.029860,8.219840,8.409820) 

y <- c(1.08181,1.08174,1.08190,1.08193,1.08191,1.08199,1.08199,
       1.08201,1.08198,1.08205,1.08200,1.08202,1.08197,1.08201,
       1.08198,1.08203,1.08204,1.08199,1.08195,1.08188,1.08193,
       1.08193,1.08197,1.08202,1.08203,1.08203,1.08207,1.08206,
       1.08214,1.08216,1.08229,1.08236,1.08240,1.08235,1.08245,
       1.08245,1.08245,1.08251,1.08263,1.08282,1.08316,1.08340,
       1.08355,1.08361,1.08373,1.08376,1.08390,1.08391)

plot(x, y)



## transform inequalities in equality constraints with slack variables
rm(list = ls())
options(device="X11")


S1 <- 100
S2 <- 80
L <- 1000
F1 <- 20
F2 <- 15
P1 <- 25
P2 <- 35







##
## Assume we want to minimize: -(0 5 0) %*% b + 1/2 b^T b
## under the constraints:      A^T b >= b0
## with b0 = (-8,2,0)^T
## and      (-4  2  0) 
##      A = (-3  1 -2)
##          ( 0  0  1)
## we can use solve.QP as follows:
##
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
qp <- solve.QP(Dmat,dvec,Amat,bvec=bvec)

bhat <- qp$solution





rm(list = ls())
options(device="X11")

Dmat <- matrix(c(6,2,2,2),2,2)
dvec <- c(-1, -6)

solve(t(Dmat)%*%Dmat, t(Dmat)%*%dvec)



M <- cbind(c(6,2,0,0),
           c(2,2,0,0),
           c(0,0,-2,-3),
           c(0,0,0,-1))
solve(t(M)%*%M, t(M)%*%c(1,6,4,0))

Amat <- matrix(c(2,1,0,3,0,1),3,2)
bvec <- c(4,0,0)
qp <- solve.QP(Dmat,dvec,t(Amat),bvec=bvec)
qp$solution

x1s <- seq(-2, 2, length=100)
l1 <- length(x1s)
x2s <- seq(-2, 2, length=99)
l2 <- length(x2s)

FX <- matrix(0, l1, l2)
rownames(FX) <- x1s
colnames(FX) <- x2s

for(i in 1:l1){
  for(j in 1:l2){
    y.hat <- Dmat%*%c(x1s[i], x2s[j])
    y <- dvec
    FX[i,j] <- as.numeric( t(y-y.hat)%*%(y-y.hat) )
  }
}

## delete below ineq constr
whiNA <- matrix(NA, l1, l2)
for(i in 1:l1){
  for(j in 1:l2){
    whiNA[i,j] <- ifelse(2*x1s[i]+3*x2s[j]>=4 & x1s[i] >=0 & x2s[j] >=0, 
                         1, NA)
  }
}


FXna <- FX*whiNA

pmin <- which(min(FXna, na.rm=TRUE)==FXna, arr.ind=TRUE)
xhat <- c(x1s[pmin[1]],x2s[pmin[2]])
2*qp$solution[1]+3*qp$solution[2]


3*xhat[1]^2 + xhat[2]^2 + 2*xhat[1]*xhat[2] + xhat[1] + 6*xhat[2] +2
xhat <- qp$solution
3*xhat[1]^2 + xhat[2]^2 + 2*xhat[1]*xhat[2] + xhat[1] + 6*xhat[2] +2

Dmat%*%xhat - dvec
Dmat%*%qp$solution- dvec


brek <- quantile(FXna, 
                 probs=c(0,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.8,1),
                 na.rm=TRUE)
colo <- rev(heat.colors(length(brek)-1))
image(x1s, x2s, FXna, breaks=brek, col=colo)
## ineq constr
x2s.eq <- seq(-7, 7, length=1000)
x1s.eq <- (4 - 3*x2s.eq)/2
lines(x1s.eq, x2s.eq, col=4, lwd=3)
abline(v=0, col=4, pch=16, lwd=3)
abline(h=0, col=4, pch=16, lwd=3)

points(qp$solution[1], qp$solution[2], pch=16, col=5, cex=1.8)
points(xhat[1], xhat[2], pch=16, col=5, cex=1.8)







LHS <- rbind(cbind(Dmat, t(Amat)),
             cbind(Amat, 0))
RHS <- rbind(-fx1, bv)
d <- solve(LHS, RHS)[1:2]




rm(list = ls())
options(device="X11")




library(quadprog)
# Simulating data from a linear model
set.seed(2019)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 0.3 * x1 - 0.3 * x2 + rnorm(n, sd=0.5)
X <- cbind(x1,x2)
# estimating a linear model with no constraints
r <- lm(y ~ -1 + x1 + x2)

## adding equality constraints

# e.g. the coefficients should sum up to b
A <- matrix(1,2,1) # constraint matrix
b <- 0.2
# using solve.QP
eqR <- solve.QP(t(X)%*%X, t(X)%*%y, A, b, meq=1)
# manual implementation 
LHS <- rbind(cbind(t(X)%*%X, A), cbind(t(A), 0))
RHS <- c(t(X)%*%y, b)
eqM <- solve(LHS, RHS)
# compare outcomes 
cbind(eqR$solution, eqM[1:2])
# check constraint
t(A)%*%eqR$solution; t(A)%*%eqM[1:2]

## adding inequality constraints

# e.g. the coefficients should be larger than b
A <- diag(2)
b <- matrix(c(0.452,-0.33),2,1)
ineqR <- solve.QP(t(X)%*%X, t(X)%*%y, A, b)
# check constraints
t(A)%*%ineqR$solution >= b
## manual 
# here I was not able to reproduce manually ineqR
# meaning solely using solve()






zeros <- matrix(0, 2,2)
mu <- 1e-6
R <- diag(2)
tXX <- t(X)%*%X
LHS1 <- cbind(tXX, zeros, t(A))
LHS2 <- cbind(tXX, mu*A, t(R))
LHS3 <- cbind(t(A), R, zeros)
LHS <- rbind(LHS1, LHS2, LHS3)
RHS <- c(t(X)%*%y, rep(0,2), b)

solve(LHS, RHS)


U <- rbind(cbind(X),
           diag(2))
LHS1 <- rbind(t(X)%*%X, U)
LHS2 <- cbind(U, matrix(0,102,2))
LHS <- rbind(LHS1,LHS2)

RHS <- c(t(X)%*%y, y, b)

solve(LHS)




# r: y - X%*%x
# min 1/2 ||r^2||
#   s.t.
# X%*%x + r = y
# A%*%x - w = b
# w >= 0

r1 <- c(2,3,2,1,0,0)
r2 <- c(1,1,2,0,1,0)
r3 <- c(-7,-8,-10,0,0,1)
M <- rbind(r1,r2,r3)
q <- c(1000,800,0)

solve(t(M)%*%M+1e-6*diag(6), t(M)%*%q)







## constrained least squares
## with (linear) constraints on fitted values
rm(list = ls())
options(device="X11")
library(quadprog)

x <- sort(runif(100))
m <- length(x)
X <- cbind(1,x,x^2)
bT <- c(1,2,-1)
yT <- X%*%bT
y <- yT + rnorm(m,0,0.1)
plot(x, y)
lines(x, yT, col=2)


tXX <- t(X)%*%X
tXy <- t(X)%*%y
bH <- solve(tXX, tXy)
yH <- X%*%bH
sum(yH)

# ## constraining the sum of the coefficients
# A <- matrix(1,1,2)
# a <- 3
# LHS <- rbind(cbind(tXX, t(A)), cbind(A,0))
# RHS <- rbind(tXy, a)
# bC <- solve(LHS, RHS)[1:2]
# A%*%bC

## constraining the sum of the fitted values
A <- matrix(1,1,m)
A <- A%*%X
a <- 100
LHS <- rbind(cbind(tXX, t(A)), cbind(A,0))
RHS <- rbind(tXy, a)
bC <- solve(LHS, RHS)[1:3]
yC <- X%*%bC
sum(yC)

cbind(bT, bH, bC)


plot(x, y)
lines(x, yT, col=2)
lines(x, yH, col=3, lwd=2, lty=2)
lines(x, yC, col=4, lwd=4)







## constrained least squares with B-splines
## with (linear) constraints on fitted values
rm(list = ls())
options(device="X11")
library(quadprog)

x <- sort(runif(100))
m <- length(x)
yT <- sin(x*5)
y <- yT + rnorm(m,0,0.1)
plot(x, y)
lines(x, yT, col=2)
sum(y)
## B-splines
B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                      ndx=10, deg=3)
nb <- ncol(B)


tBB <- t(B)%*%B
tBy <- t(B)%*%y
alphaH <- solve(tBB, tBy)
yH <- B%*%alphaH
sum(yH)

## constraining the sum of the fitted values
A <- matrix(1,1,m)
A <- A%*%B
a <- sum(yT)*2.5
LHS <- rbind(cbind(tBB, t(A)), cbind(A,0))
RHS <- rbind(tBy, a)
coef <- solve(LHS, RHS)
alphaC <- coef[1:nb]
yC <- B%*%alphaC
sum(yC)

cbind(alphaH, alphaC)

cbind(sum(yT), sum(yH), sum(yC))

plot(x, y)
lines(x, yT, col=2)
lines(x, yH, col=3, lwd=2, lty=2)
lines(x, yC, col=4, lwd=4)

quantile(yH-yC)
coef[nb+1]
