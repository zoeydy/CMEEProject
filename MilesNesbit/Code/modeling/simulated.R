require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)
library(dplyr)
library(growthcurver)
library(stringr)
library(taxize)
library(nls.multstart)
library(zoo)
library(broom)
library(growthcurver)
library(cowplot)
library(psych)
library(MASS)
library(fitdistrplus)


rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")

B_0 = 10
k = 1
E = 2.5
t = 1:100
B0_funct <- function(t,k,E,B_0){
  return(B_0*exp(-(E/(k*t))))
}

y_pred = B0_funct(t, k, E, B_0)
plot(t, y_pred)

b=100
x = (sample.int(101,size=100,replace=TRUE)-1)
m=-1

trade_funct <- function(m,x,b){
  return((m*x)+b)
}

y_pred = trade_funct(m,x,b)
y_pred = as.integer(y_pred)
x1= as.integer(x)

fakedata <- data.frame(y_pred, x1)

ggplot(fakedata, aes(x=x1,y=y_pred))+
  geom_point(position= 'jitter', size= 3)+ 
  ggtitle("Expected r vs K tradeoff")+
  xlab('predicted r')+
  ylab('predicted K')+
  geom_smooth(method='lm', formula= y~x)

