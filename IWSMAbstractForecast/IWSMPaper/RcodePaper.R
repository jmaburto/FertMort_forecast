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
load("/home/camarda_cg/WORK/FertMort/trunk/MortalityExamples/ITAfemMort.RData")

ranx <- range(years)
rany <- range(e01, e0.tar, e0.S)
pdf("CamardaE0.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Life expectancy over time"), 
      cex.main=1.8)
mtext(expression(paste(e^0, " (years)")), 
      2, cex=1.8, line=2.1)
mtext(expression(paste("year, ", t)), 
      1, cex=1.8, line=2.8, las=1)
abline(h=seq(-200, 100, 5), col="grey80", lty=3)
abline(v=seq(1900, 2100, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()

points(years1, e01, cex=1.3, pch=16, col="grey50")
points(yearsF, e0.tar, pch=3, cex=1, lwd=2)
lines(years, e0.C, col=4, lwd=4)
lines(years, e0.S, col=2, lwd=3, lty=2)
lege1 <- expression(paste("Actual"))
lege2 <- expression(paste("Future values from the UN"))
lege3 <- expression(paste("2D P-splines"))
lege4 <- expression(paste("SQP + 2D P-splines"))
abline(v=years[n1]+0.5, col="grey50")
legend("bottomright", inset=0.05, 
       c(lege1, lege2, lege3, lege4),
       lwd=c(1,2,3,4), lty=c(NA, NA, 2,1), col=c("grey50", 1,2,4),
       bg="grey95", cex=1.7,
       pch=c(16,3, NA, NA), pt.cex=c(1.3, 1,1,1))
dev.off()


ranx <- range(years)
rany <- range(ed1, ed.tar, ed.S)
pdf("CamardaEd.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Lifespan variation over time"), 
      cex.main=1.8)
mtext(expression(paste(e, " (years)")), 
      2, cex=1.8, line=2.1)
mtext(expression(paste("year, ", t)), 
      1, cex=1.8, line=2.8, las=1)
abline(h=seq(-200, 100, 1), col="grey80", lty=3)
abline(v=seq(1900, 2100, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()

points(years1, ed1, cex=1.3, pch=16, col="grey50")
points(yearsF, ed.tar, pch=3, cex=1, lwd=2)
lines(years, ed.C, col=4, lwd=4)
lines(years, ed.S, col=2, lwd=3, lty=2)
lege1 <- expression(paste("Actual"))
lege2 <- expression(paste("Future values from time-series"))
lege3 <- expression(paste("2D P-splines"))
lege4 <- expression(paste("SQP + 2D P-splines"))
abline(v=years[n1]+0.5, col="grey50")
legend("topright", inset=0.05, 
       c(lege1, lege2, lege3, lege4),
       lwd=c(1,2,3,4), lty=c(NA, NA, 2,1), col=c("grey50", 1,2,4),
       bg="grey95", cex=1.7,
       pch=c(16,3, NA, NA), pt.cex=c(1.3, 1,1,1))
dev.off()


i=n
ranx <- range(ages)
rany <- range(ETA.S[,i], ETA.C[,i], na.rm=TRUE, finite=TRUE)
pdf("CamardaMort2050.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("log-mortality, year", years[i]), 
      cex.main=1.8)
mtext(expression(eta), 2, cex=1.8, line=2.8, las=1)
mtext(expression(paste("age, ", a)), 1, cex=1.8, 
      line=2.8, las=1)
abline(h=seq(-20, 10, 2), col="grey80", lty=3)
abline(v=seq(0, 120, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()
lines(ages, ETA.S[,i], col=2, lwd=3, lty=2)
lines(ages, ETA.C[,i], col=4, lwd=4)
lege1 <- expression(paste("2D P-splines, ", e^0, "=92.2"))
lege2 <- expression(paste("SQP + 2D P-splines, ", e^0, "=89.47"))
legend("topleft", inset=0.05, c(lege1, lege2),
       lwd=c(3,4), lty=c(2,1), col=c(2,4),
       bg="grey95", cex=1.7)
dev.off()                    

ranx <- range(50:m)
rany <- c(0,0.055)
pdf("CamardaDens2050.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("age-at-death distribution, year", years[i]), 
      cex.main=1.8)
mtext(expression("density"), 2, cex=1.8, line=2.8, las=0)
mtext(expression(paste("age, ", a)), 1, cex=1.8, 
      line=2.8, las=1)
abline(h=seq(-20, 10, 0.01), col="grey80", lty=3)
abline(v=seq(0, 120, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()
lines(ages, F.C[,i], col=4, lwd=4)
lines(ages, F.S[,i], col=2, lwd=3, lty=2)
lege1 <- expression(paste("2D P-splines"))
lege2 <- expression(paste("SQP + 2D P-splines"))
legend("topleft", inset=0.05, c(lege1, lege2),
       lwd=c(3,4), lty=c(2,1), col=c(2,4),
       bg="grey95", cex=1.7)
dev.off()



rm(list = ls())
options(device="X11")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(quadprog)
library(rgl)
library(MortalitySmooth)
library(HMDdata)
library(vars)
library(forecast)
load("/home/camarda_cg/WORK/FertMort/trunk/FertilityExamples/ESPfert.RData")


ranx <- range(t)
rany <- range(TFR1, TFR.S, TFR.C, TFR.tar)
pdf("CamardaTFR.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Total Fertility Rate over time"), 
      cex.main=1.8)
mtext(expression(paste(TFR, " (children per woman)")), 
      2, cex=1.8, line=2.1)
mtext(expression(paste("year, ", t)), 
      1, cex=1.8, line=2.8, las=1)
abline(h=seq(-200, 100, 0.5), col="grey80", lty=3)
abline(v=seq(1900, 2100, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()

points(t1, TFR1, cex=1.3, pch=16, col="grey50")
points(tF, TFR.tar, pch=3, cex=1, lwd=2)
lines(t, TFR.C, col=4, lwd=4)
lines(t, TFR.S, col=2, lwd=3, lty=2)
lege1 <- expression(paste("Actual"))
lege2 <- expression(paste("Future values from time-series"))
lege3 <- expression(paste("2D P-splines"))
lege4 <- expression(paste("SQP + 2D P-splines"))
abline(v=t[n1]+0.5, col="grey50")
legend("topright", #inset=0.05,
       c(lege1, lege2, lege3, lege4),
       lwd=c(1,2,3,4), lty=c(NA, NA, 2,1), col=c("grey50", 1,2,4),
       bg="grey95", cex=1.7,
       pch=c(16,3, NA, NA), pt.cex=c(1.3, 1,1,1))
dev.off()

ranx <- range(t)
rany <- range(MAB1, MAB.S, MAB.C, MAB.tar)
pdf("CamardaMAB.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Mean age at childbearing over time"), 
      cex.main=1.8)
mtext(expression(paste(MAB, " (years)")), 
      2, cex=1.8, line=2.1)
mtext(expression(paste("year, ", t)), 
      1, cex=1.8, line=2.8, las=1)
abline(h=seq(-200, 100, 1), col="grey80", lty=3)
abline(v=seq(1900, 2100, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()

points(t1, MAB1, cex=1.3, pch=16, col="grey50")
points(tF, MAB.tar, pch=3, cex=1, lwd=2)
lines(t, MAB.C, col=4, lwd=4)
lines(t, MAB.S, col=2, lwd=3, lty=2)
lege1 <- expression(paste("Actual"))
lege2 <- expression(paste("Future values from time-series"))
lege3 <- expression(paste("2D P-splines"))
lege4 <- expression(paste("SQP + 2D P-splines"))
abline(v=t[n1]+0.5, col="grey50")
legend("topleft", #inset=0.05, 
       c(lege1, lege2, lege3, lege4),
       lwd=c(1,2,3,4), lty=c(NA, NA, 2,1), col=c("grey50", 1,2,4),
       bg="grey95", cex=1.7,
       pch=c(16,3, NA, NA), pt.cex=c(1.3, 1,1,1))
dev.off()


ranx <- range(t)
rany <- range(VAB1, VAB.S, VAB.C, VAB.tar)
pdf("CamardaVAB.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Variance in childbearing age over time"), 
      cex.main=1.8)
mtext(expression(paste(VAB, " (",years^2,")")), 
      2, cex=1.8, line=2.1)
mtext(expression(paste("year, ", t)), 
      1, cex=1.8, line=2.8, las=1)
abline(h=seq(-200, 100, 2), col="grey80", lty=3)
abline(v=seq(1900, 2100, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()

points(t1, VAB1, cex=1.3, pch=16, col="grey50")
points(tF, VAB.tar, pch=3, cex=1, lwd=2)
lines(t, VAB.C, col=4, lwd=4)
lines(t, VAB.S, col=2, lwd=3, lty=2)
lege1 <- expression(paste("Actual"))
lege2 <- expression(paste("Future values from time-series"))
lege3 <- expression(paste("2D P-splines"))
lege4 <- expression(paste("SQP + 2D P-splines"))
abline(v=t[n1]+0.5, col="grey50")
legend("bottomright", #inset=0.05, 
       c(lege1, lege2, lege3, lege4),
       lwd=c(1,2,3,4), lty=c(NA, NA, 2,1), col=c("grey50", 1,2,4),
       bg="grey95", cex=1.7,
       pch=c(16,3, NA, NA), pt.cex=c(1.3, 1,1,1))
dev.off()

i=n
ranx <- range(x)
rany <- range(MU.S[,i], MU.C[,i], na.rm=TRUE, finite=TRUE)
pdf("CamardaFert2050.pdf", width=8, height=8)
plot(1,1,t="n", xlim=ranx, ylim=rany, 
     xlab="", ylab="",
     axes=FALSE)
title(main=paste("Age-specific fertility rates, year", t[i]), 
      cex.main=1.8)
mtext(expression(mu), 2, cex=1.8, line=2.8, las=1)
mtext(expression(paste("age, ", a)), 1, cex=1.8, 
      line=2.8, las=1)
abline(h=seq(-20, 10, 0.05), col="grey80", lty=3)
abline(v=seq(0, 120, 20), col="grey80", lty=3)
axis(1, cex.axis=1.5)
axis(2, las=2, cex.axis=1.5);box()
lines(x, MU.S[,i], col=2, lwd=3, lty=2)
lines(x, MU.C[,i], col=4, lwd=4)
lege1 <- expression(paste("2D P-splines"))
lege2 <- expression(paste("SQP + 2D P-splines"))
legend("topleft", inset=0.05, c(lege1, lege2),
       lwd=c(3,4), lty=c(2,1), col=c(2,4),
       bg="grey95", cex=1.7)
dev.off()                    
