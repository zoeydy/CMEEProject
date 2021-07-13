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

#create log
files$PopBioLog <- log(files$PopBio)

log(0.283275711)

#reorder
col_order <- c("Time", "PopBio", "PopBioLog",
               "Temp", "Time_units", 'PopBio_unts', 'Species', 'Medium', 'Rep', 'Citation')
my_data2 <- files[, col_order]


#save
write.csv(my_data2,"../Results/Formated/Total_processedLog.csv", row.names = FALSE)
