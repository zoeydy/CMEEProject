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
files <- read.csv('../Results/Formated/Total_processedLog.csv')


#kill pop bio
files <- select(files,-c(PopBio))

#subsetting
bae_cb <-subset(subset(files, Temp == 5), Species == 'Chryseobacterium.balustinum')
drops <- c("Citation", "Temp", 'Species', 'Time_units', 'PopBio_unts', 'Medium', 'Rep')
bae_cb <- bae_cb[ , !(names(bae_cb) %in% drops)]


#specify the model functions:
gompertz_model <- function(t, r_max, N_max, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}


#mult start
fit_gompertz_multi <- nls_multstart(PopBioLog ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag),
                                    data = bae_cb,
                                    start_lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 0),
                                    start_upper = c(t_lag=20, r_max=10, N_0 = 6, N_max = 10),
                                    lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 1),
                                    iter = 500,
                                    supp_errors = "Y")
print(fit_gompertz_multi)


