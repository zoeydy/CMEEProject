setwd("/home/nesbit/Desktop/R/TheMulQuaBio-master/data")
install.packages("minpack.lm")
require("minpack.lm")

rm(list = ls())
graphics.off()

powMod <- function(x, a, b) {
  return(a * x^b)
}
MyData <- read.csv("../data/GenomeSize.csv")

head(MyData)
Data2Fit <- subset(MyData,Suborder == "Anisoptera")

Data2Fit <- Data2Fit[!is.na(Data2Fit$TotalLength), ]

plot(Data2Fit$TotalLength, Data2Fit$BodyWeight)
library("ggplot2")

ggplot(Data2Fit, aes(x = TotalLength, y = BodyWeight)) + 
  geom_point(size = (3),color="red") + theme_bw() + 
  labs(y="Body mass (mg)", x = "Wing length (mm)")
PowFit <- nlsLM(BodyWeight ~ powMod(TotalLength, a, b), data = Data2Fit, start = list(a = .1, b = .1))
summary(PowFit)
Lengths <- seq(min(Data2Fit$TotalLength),max(Data2Fit$TotalLength),len=200)
coef(PowFit)["a"]
coef(PowFit)["b"]
Predic2PlotPow <- powMod(Lengths,coef(PowFit)["a"],coef(PowFit)["b"])
plot(Data2Fit$TotalLength, Data2Fit$BodyWeight)
lines(Lengths, Predic2PlotPow, col = 'blue', lwd = 2.5)
confint(PowFit)
alb <- read.csv(file="../data/albatross_grow.csv")
alb <- subset(x=alb, !is.na(alb$wt))
plot(alb$age, alb$wt, xlab="age (days)", ylab="weight (g)", xlim=c(0, 100))
logistic1<-function(t, r, K, N0){
  N0*K*exp(r*t)/(K+N0*(exp(r*t)-1))
}

vonbert.w<-function(t, Winf, c, K){
  Winf*(1 - exp(-K*t) + c*exp(-K*t))^3
}
scale<-4000

alb.lin<-lm(wt/scale~age, data=alb)
alb.log<-nlsLM(wt/scale~logistic1(age, r, K, N0), start=list(K=1, r=0.1, N0=0.1), data=alb)

alb.vb<-nlsLM(wt/scale~vonbert.w(age, Winf, c, K), start=list(Winf=0.75, c=0.01, K=0.01), data=alb)
ages<-seq(0, 100, length=1000)

pred.lin<-predict(alb.lin, newdata = list(age=ages))*scale

pred.log<-predict(alb.log, newdata = list(age=ages))*scale

pred.vb<-predict(alb.vb, newdata = list(age=ages))*scale
plot(alb$age, alb$wt, xlab="age (days)", ylab="weight (g)", xlim=c(0,100))
lines(ages, pred.lin, col=2, lwd=2)
lines(ages, pred.log, col=3, lwd=2)
lines(ages, pred.vb, col=4, lwd=2)

legend("topleft", legend = c("linear", "logistic", "Von Bert"), lwd=2, lty=1, col=2:4)

par(mfrow=c(3,1), bty="n")
plot(alb$age, resid(alb.lin), main="LM resids", xlim=c(0,100))
plot(alb$age, resid(alb.log), main="Logisitic resids", xlim=c(0,100))
plot(alb$age, resid(alb.vb), main="VB resids", xlim=c(0,100))
n<-length(alb$wt)
list(lin=signif(sum(resid(alb.lin)^2)/(n-2*2), 3), 
     log= signif(sum(resid(alb.log)^2)/(n-2*3), 3), 
     vb= signif(sum(resid(alb.vb)^2)/(n-2*3), 3))
