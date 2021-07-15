setwd("/home/nesbit/Desktop/Data_RVK/Code")
library(tidyverse)
library(ggplot2)
require("minpack.lm")
rm(list=ls())
graphics.off()
alb <- read.csv(file="../Data/All_Figs/Bae_et_al_2014_processed.csv")
alb <- subset(x=alb, !is.na(alb$PopBio))
plot(alb$Time, alb$PopBio, xlab="age (H)", ylab="PopBio", xlim=c(0, 700))
    

    logistic1<-function(t, r, K, N0){
          N0*K*exp(r*t)/(K+N0*(exp(r*t)-1))
        }

    scale<-4000
    for(i in 1:22){
      df_subset <- subset(alb, alb$Species == unique(alb$Species)[i])    
      alb.log<-nlsLM(PopBio/scale~logistic1(Time, r, K, N0), start=list(K=1, r=0.1, N0=0.1), data=df_subset)
      ages<-seq(0, 650, length= q )
      q<-ncol(pred.log)
      pred.log<-predict(alb.log, newdata = list(age=ages))*scale
    plotted <- plot(df_subset$Time, df_subset$PopBio, xlab="age (H)", ylab="PopBio", xlim=c(0, 650))
      lines(ages, pred.log, col=3, lwd=2)
      print(plotted)
      }
    




length(pred.log)
