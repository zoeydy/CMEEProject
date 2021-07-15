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
files <- read.csv("../Results/Formated/3Total_processed_Log_RVK_3_models.csv")


med <-(unique(files$Medium))

write.csv(med, "../Results/Formated/mediums.csv")

mediums<- read.csv("../Results/Formated/mediumsnutrients.csv")

#clean
drops <- c("X", "Citation")
mediums <- mediums[ , !(names(mediums) %in% drops)]

#bind
total <- merge(files,mediums,by="Medium")

#save
write.csv(total, "../Results/Formated/4Total_processed_Log_RVK_3_models.csv")

