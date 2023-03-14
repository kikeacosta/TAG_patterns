## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)

path <- "/home/gccamarda/WORK/Applications/TAG_COVID19_UN/2021_04_09_TAG/ModelPattern/TAG_patterns/Slides/Figures/"


## individual
load("OUTCOMES/ASIindiv.RData")


pdf(paste(path, "AsiaIndiv.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(3,4))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(3,4,3,1))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=2)
  yy <- 10^seq(-10, 0)
  axis(2, at=log(yy), labels=yy, las=2)
  axis(1);box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
  }
  xxx <- c(x, rev(x))
  yyy <- c(ETA.LOW[,i], rev(ETA.UP[,i]))
  polygon(xxx, yyy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
  lines(x, ETA[,i], col=colocou[i], lwd=4)
  #locator(1)
}
dev.off()


## all together: only fit w/o CI
pdf(paste(path, "AsiaIndivFIT.pdf", sep=""), width = 14, height = 9)
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="ASIAN POPULATIONS", col.main=1, cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  lines(x, ETA[,i], col=colocou[i], lwd=4)
  #locator(1)
}
legend("topleft", inset=0.05,
       legend=REG, col=colocou, lwd=4, bg="grey80", ncol=2, cex=1.5)
dev.off()

## Region
load("OUTCOMES/ASIregion.RData")


pdf(paste(path, "AsiaRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(3,4))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(3,4,3,1))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=2)
  yy <- 10^seq(-10, 0)
  axis(2, at=log(yy), labels=yy, las=2)
  axis(1);box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
  }
  xxx <- c(x, rev(x))
  yyy <- c(ETA.LOW[,i], rev(ETA.UP[,i]))
  polygon(xxx, yyy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
  lines(x, ETA[,i], col=colocou[i], lwd=4)
  #locator(1)
}
dev.off()


## all together: only fit w/o CI
pdf(paste(path, "AsiaRegionFIT.pdf", sep=""), width = 14, height = 9)
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="ASIAN POPULATIONS", col.main=1, cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  lines(x, ETA[,i], col=colocou[i], lwd=4)
  #locator(1)
}
legend("topleft", inset=0.05,
       legend=REG, col=colocou, lwd=4, bg="grey80", ncol=2, cex=1.5)
dev.off()


load("OUTCOMES/ASIindiv.RData")
ETA.IND <- ETA
load("OUTCOMES/ASIregion.RData")
ETA.REG <- ETA

pdf(paste(path, "AsiaIndivRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(3,4))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(3,4,3,1))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=2)
  yy <- 10^seq(-10, 0)
  axis(2, at=log(yy), labels=yy, las=2)
  axis(1);box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=2)
  }
  lines(x, ETA.IND[,i], col="darkred", lwd=4)
  lines(x, ETA.REG[,i], col="darkblue", lwd=4)
  #locator(1)
}
dev.off()








load("OUTCOMES/EURindiv.RData")
ETA.IND <- ETA
load("OUTCOMES/EURregion.RData")
ETA.REG <- ETA

pdf(paste(path, "EuropeIndivRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(6,4))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(0.5,0.5,1.5,0.5))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=1.5)
  yy <- 10^seq(-10, 0)
  #axis(2, at=log(yy), labels=yy, las=2)
  #axis(1)
  box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=1)
  }
  lines(x, ETA.IND[,i], col="darkred", lwd=3)
  lines(x, ETA.REG[,i], col="darkblue", lwd=3)
  #locator(1)
}
dev.off()



load("OUTCOMES/AFRindiv.RData")
ETA.IND <- ETA
load("OUTCOMES/AFRregion.RData")
ETA.REG <- ETA

pdf(paste(path, "AfricaIndivRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(2,4))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(0.5,0.5,1.5,0.5))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=1.5)
  yy <- 10^seq(-10, 0)
  #axis(2, at=log(yy), labels=yy, las=2)
  #axis(1)
  box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=2)
  }
  lines(x, ETA.IND[,i], col="darkred", lwd=4)
  lines(x, ETA.REG[,i], col="darkblue", lwd=4)
  #locator(1)
}
dev.off()


load("OUTCOMES/LAMindiv.RData")
ETA.IND <- ETA
load("OUTCOMES/LAMregion.RData")
ETA.REG <- ETA

pdf(paste(path, "LatinAmericaIndivRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(3,5))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(0.5,0.5,1.5,0.5))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=1.5)
  yy <- 10^seq(-10, 0)
  #axis(2, at=log(yy), labels=yy, las=2)
  #axis(1)
  box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=2)
  }
  lines(x, ETA.IND[,i], col="darkred", lwd=4)
  lines(x, ETA.REG[,i], col="darkblue", lwd=4)
  #locator(1)
}
dev.off()


load("OUTCOMES/NAMindiv.RData")
ETA.IND <- ETA
load("OUTCOMES/NAMregion.RData")
ETA.REG <- ETA

pdf(paste(path, "NorthAmericaIndivRegion.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(1,3))
for(i in 1:nc){
  colocou <- rainbow_hcl(nc)
  rany <- c(-19, -3)
  ranx <- range(x)
  par(mar=c(0.5,0.5,1.5,0.5))
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=FALSE,xlab="", ylab="")
  #mtext("age", 1, cex=2, line=3.5)
  #mtext("mortality, log-scale", 2, cex=2, line=3.5)
  title(main=REG[i], col.main=colocou[i], cex.main=1.5)
  yy <- 10^seq(-10, 0)
  #axis(2, at=log(yy), labels=yy, las=2)
  #axis(1)
  box()
  abline(h=log(yy), col = "lightgray", lty = "dotted")
  abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
  xx <- c(x, rev(x))
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=2)
  }
  lines(x, ETA.IND[,i], col="darkred", lwd=4)
  lines(x, ETA.REG[,i], col="darkblue", lwd=4)
  #locator(1)
}
dev.off()





pdf(paste(path, "WorldIndiv.pdf", sep=""), width = 14, height = 9)
lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
coloreg <- brewer.pal(n = 6, name = "Dark2")
par(mfrow=c(2,3))

rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[1], col.main=coloreg[1], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
library(RColorBrewer)
coloreg <- brewer.pal(n = 6, name = "Dark2")#c("#89C5DA", "#DA5724", "#8A7C64", "#74D944", "#CE50CA", "#C0717C")
load("OUTCOMES/EURindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[1], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[2], col.main=coloreg[2], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/AFRindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[2], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[3], col.main=coloreg[3], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/ASIindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[3], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[4], col.main=coloreg[4], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/LAMindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[4], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[5], col.main=coloreg[5], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/NAMindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[5], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[6], col.main=coloreg[6], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/OCEindiv.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[6], lwd=3)
}
# lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
# legend("topleft", inset=0.05,
#        legend=lege, col=coloreg, lwd=3, cex=1.5)
dev.off()



pdf(paste(path, "WorldRegion.pdf", sep=""), width = 14, height = 9)
lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
coloreg <- brewer.pal(n = 6, name = "Dark2")
par(mfrow=c(2,3))

rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[1], col.main=coloreg[1], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
library(RColorBrewer)
coloreg <- brewer.pal(n = 6, name = "Dark2")#c("#89C5DA", "#DA5724", "#8A7C64", "#74D944", "#CE50CA", "#C0717C")
load("OUTCOMES/EURregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[1], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[2], col.main=coloreg[2], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/AFRregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[2], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[3], col.main=coloreg[3], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/ASIregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[3], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[4], col.main=coloreg[4], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/LAMregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[4], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[5], col.main=coloreg[5], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/NAMregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[5], lwd=3)
}
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[6], col.main=coloreg[6], cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
load("OUTCOMES/OCEregion.RData")
for(i in 1:nc){
  lines(x, ETA[,i], col=coloreg[6], lwd=3)
}
# lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
# legend("topleft", inset=0.05,
#        legend=lege, col=coloreg, lwd=3, cex=1.5)
dev.off()

## derivatives
source("EmpiricDeriv.R")


pdf(paste(path, "WorldROA.pdf", sep=""), width = 15, height = 8)
lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
coloreg <- brewer.pal(n = 6, name = "Dark2")

par(mfrow=c(1,2))
rany <- c(-19, -3)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Reference age-pattern", col.main=1, cex.main=1.5)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
library(RColorBrewer)
coloreg <- brewer.pal(n = 6, name = "Dark2")#c("#89C5DA", "#DA5724", "#8A7C64", "#74D944", "#CE50CA", "#C0717C")
load("OUTCOMES/EURregion.RData")
lines(x, ETA[,1], col=coloreg[1], lwd=4)
load("OUTCOMES/AFRregion.RData")
lines(x, ETA[,1], col=coloreg[2], lwd=4)
load("OUTCOMES/ASIregion.RData")
lines(x, ETA[,1], col=coloreg[3], lwd=4)
load("OUTCOMES/LAMregion.RData")
lines(x, ETA[,1], col=coloreg[4], lwd=4)
load("OUTCOMES/NAMregion.RData")
lines(x, ETA[,1], col=coloreg[5], lwd=4)
load("OUTCOMES/OCEregion.RData")
lines(x, ETA[,1], col=coloreg[6], lwd=4)
lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
legend("bottomright", inset=0.05,
      legend=lege, col=coloreg, lwd=3, cex=1.5, bg="grey80")



par(mfrow=c(3,2))
rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Rate of aging", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/EURregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[1], lwd=4)


load("OUTCOMES/AFRregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[2], lwd=4)
load("OUTCOMES/ASIregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[3], lwd=4)
load("OUTCOMES/LAMregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[4], lwd=4)
load("OUTCOMES/NAMregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[5], lwd=4)
load("OUTCOMES/OCEregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col=coloreg[6], lwd=4)

lege <- c("Europe", "Africa", "Asia", "Latin America", "North America", "Oceania")
legend("bottomright", inset=0.05,
       legend=lege, col=coloreg, lwd=3, cex=1.5, bg="grey80")
dev.off()



pdf(paste(path, "WorldROAdiff.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(2,3))
rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[1], col.main=coloreg[1], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/EURregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/EURindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)

rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[2], col.main=coloreg[2], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/AFRregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/AFRindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)

rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[3], col.main=coloreg[3], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/ASIregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/ASIindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)


rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[4], col.main=coloreg[4], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/LAMregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/LAMindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)


rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[5], col.main=coloreg[5], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/NAMregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/NAMindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)

rany <- c(-0.2, 0.2)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main=lege[6], col.main=coloreg[6], cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
load("OUTCOMES/OCEregion.RData")
eta1 <- emp.der(x=x, y=ETA[,1])
lines(x, eta1, col="darkblue", lwd=4)
load("OUTCOMES/OCEindiv.RData")
eta.mean <- apply(ETA, 1, mean)
eta1.mean <- emp.der(x=x, y=eta.mean)
lines(x, eta1.mean, col="darkred", lwd=4)
dev.off()



## END