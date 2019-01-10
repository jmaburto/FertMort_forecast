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

## loading age-sepcific-fertility-rates from HFD (France)
ASFR0 <- read.table("FRATNPasfrRR.txt", header=TRUE, skip=2)
## dimensions
y <- unique(ASFR0$Year)
x <- 12:55
m <- length(x)
n <- length(y)
## asfr in matrix
ASFR <- matrix(ASFR0$ASFR, m, n)

## estimating Hadwiger curve for each year
## in a really simplistic way with OLS
objFUN <- function(par, asfr.obs, x){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  asfr.hat <- Hadwiger(a,b,c,x)
  RSS <- sum((asfr.hat - asfr.obs)^2)
  return(RSS)
}
PAR.hat <- matrix(0, n, 3)
colnames(PAR.hat) <- letters[1:3]
rownames(PAR.hat) <- y

st.val <- c(1.5, 3.8, 28)
for(i in 1:n){
  opt <- optim(st.val, objFUN, asfr.obs=ASFR[,i], x=x)
  st.val <- opt$par
  PAR.hat[i,] <- opt$par
  # plot(x, ASFR[,i])
  # lines(x, Hadwiger(opt$par[1], opt$par[2], opt$par[3], x=x), col=2)
  # locator(1)
}
par(mfrow=c(3,1))

## fitted parameters
plot(y, PAR.hat[,1])
plot(y, PAR.hat[,2])
plot(y, PAR.hat[,3])
par(mfrow=c(1,1))

## smoothing (again in a simplistic way) fitted parameters
## since we will use them for simulating data
PAR <- PAR.hat*0
PAR[,1] <- loess(PAR.hat[,1]~y, span=0.2)$fitted
PAR[,2] <- loess(PAR.hat[,2]~y, span=0.2)$fitted
PAR[,3] <- loess(PAR.hat[,3]~y, span=0.2)$fitted
  
## plotting estimated and smooth parameters
par(mfrow=c(3,1))
plot(y, PAR.hat[,1])
lines(y, PAR[,1], col=2)
plot(y, PAR.hat[,2])
lines(y, PAR[,2], col=2)
plot(y, PAR.hat[,3])
lines(y, PAR[,3], col=2)
par(mfrow=c(1,1))

## save smooth parameters
write.table(PAR, "PAR.Hadwiger.txt")


# ## how to read it
# bla <- read.table("PAR.Hadwiger.txt", header=TRUE)

