require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)
library(dplyr)
library(growthcurver)
library(stringr)
library(taxize)
library(nls.multstart)
library(zoo)

rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
files <- read.csv('../Results/Formated/Total_processed.csv')


#clean
drops <- c("X", "X.1")
files <- files[ , !(names(files) %in% drops)]


#select smith
df <- subset(files, Citation == 'Smith_et_al')

#divide time by 60
df$Time <- df$Time/60

#drop smith from files
files <- subset(files, Citation != 'Smith_et_al')
files <- rbind(files, df)

write.csv(files,"../Results/Formated/Total_processed.csv", row.names = FALSE)
