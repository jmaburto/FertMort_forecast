rm(list = ls())

library(MortalitySmooth)
library(HMDdata)


x <- 0:100
m <- length(x)
y <- 1838:2013
n <- length(y)

D <- selectHMDdata("Island",
                   "Death", "Males", x, y)
E <- selectHMDdata("Island",
                   "Exposure", "Males", x, y)
lMX <- log(D/E)


Wei <- matrix(1, m, n)
Wei[E==0] <- 0

fit <- Mort2Dsmooth(x=x, y=y,
                    Z=D, offset=log(E),
                    W=Wei),
                    method=3,
                    lambdas=c(10^-3.5, 10^2))
plot(fit)

lMX.hat <- fit$logmortality

ran <- c(-12, 2)
for(i in 1:n){
  plot(x, lMX[,i], ylim=ran,
       main=paste("year :", y[i]))
  lines(x, lMX.hat[,i], col=2, lwd=3)
  ## locator(1)
  Sys.sleep(0.1)
}

for(i in 1:m){
  ran <- range(lMX[i,], lMX.hat[i,],
               na.rm=TRUE, finite=TRUE)
  plot(y, lMX[i,], ylim=ran,
       main=paste("age :", x[i]))
  lines(y, lMX.hat[i,], col=2, lwd=3)
  ## locator(1)
  Sys.sleep(0.1)
}


































## END
