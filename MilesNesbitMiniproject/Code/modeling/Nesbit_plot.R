# install.packages("devtools")
# require(devtools)
# install_version("cowplot", version = "0.9.0", repos = "http://cran.us.r-project.org")
library(tidyverse)
library(broom)
library(growthcurver)
library(cowplot)
library(psych)
library(MASS)
library(fitdistrplus)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")

final_df <- read.csv('../Results/Formated/Total_processed.csv')

drops <- c("X", "X.1")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#plotting series
final_df %>%
  ggplot(aes(x=Time,y=PopBio,colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~paste(Species,Medium, Rep), scales = "free")

