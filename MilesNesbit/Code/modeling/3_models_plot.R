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
files <- read.csv("../Results/Formated/3Total_processed_Log_RVK_3_models.csv")

#clean the data
files$r <- NA
files$K <- NA
files$t_lag <- NA
files$N_0 <- NA
files$r.t.value <- NA
files$K.t.value <- NA
files$t_lag.t.value <- NA
files$N_0.t.value <- NA

for (i in unique(files$discrete_curve)){
  
  # Subset data by curve number
  ff_sub <- subset(files, discrete_curve == i)
  
  # Find estimated values
  r <- subset(ff_sub, parameter == "r_max")$Estimate[1]
  r_tval <- subset(ff_sub, parameter == "r_max")$t.value[1]
  Nmax <- subset(ff_sub, parameter == "N_max")$Estimate[1]
  Nmax_tval <- subset(ff_sub, parameter == "N_max")$t.value[1]
  N0 <- subset(ff_sub, parameter == "N_0")$Estimate[1]
  N0_tval <- subset(ff_sub, parameter == "N_0")$t.value[1]
  tlag <- subset(ff_sub, parameter == "t_lag")$Estimate[1]
  tlag_tval <- subset(ff_sub, parameter == "t_lag")$t.value[1]
  
  # Add data 
  index <- which(files$discrete_curve == i)
  for (j in index){
    files$r[j] <- r
    files$r.t.value[j] <- r_tval
    files$K[j] <- Nmax
    files$K.t.value[j] <- Nmax_tval
    files$N_0[j] <- N0
    files$N_0.t.value[j] <- N0_tval
    files$t_lag[j] <- tlag
    files$t_lag.t.value[j] <- tlag_tval
  }
  
}

#nest
ff <- files %>%
   nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts, -Estimate,
        -t.value, -discrete_curve, -model, -parameter, -r, -K, -t_lag, -N_0,-r.t.value,
        -K.t.value,-N_0.t.value,-t_lag.t.value, data = c("Time","PopBio"))

# Tidy ff
ff_new <- subset(ff, parameter == "r_max")


#drops
drops <- c("parameter", "t.value", "Estimate")
ff_new <- ff_new[ , !(names(ff_new) %in% drops)]


#try the plot
Absorbance <- ff_new %>%
  filter(str_detect(PopBio_unts, 'Abs'))

ABSD <- ff_new %>%
  filter(str_detect(PopBio_unts, 'OD'))

Absorbance <- rbind(Absorbance, ABSD)

drops2 <- c("data")
ff_new <- ff_new[ , !(names(ff_new) %in% drops2)]


#save ff_new
write.csv(ff_new,"../Results/Formated/4r_and_K.csv", row.names = FALSE)


#plot
ggplot(Absorbance, aes(x= log(r), y= K, color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Absorbance")
