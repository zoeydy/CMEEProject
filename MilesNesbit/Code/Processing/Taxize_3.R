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

fullset <- read.csv("../Results/Formated/4Total_processed_Log_RVK_3_models.csv", stringsAsFactors = FALSE)
files2 <- read.csv('../Results/Formated/Total_processed.csv', stringsAsFactors = FALSE)

#clean
drops <- c("X", "X.1")
files <- files[ , !(names(files) %in% drops)]

#add empty column to df

df1 <- files2

files2 <- read.csv('../Results/Formated/Total_processed.csv')
files2 <- files2[ , !(names(files2) %in% drops)]



files <- read.csv('../Results/Formated/SpeciesList.csv', stringsAsFactors = FALSE)

df <- read.csv('../Results/Formated/FullListSpecies.csv', stringsAsFactors = FALSE)

df2 <- read.csv('../Results/Formated/FullListSpeciesStart.csv', stringsAsFactors = FALSE)

df3 <- read.csv('../Results/Formated/FullListSpeciesStart2.csv', stringsAsFactors = FALSE)

df4 <- read.csv('../Results/Formated/FullListSpeciesStart3.csv', stringsAsFactors = FALSE)


#clean df
df <- df %>% 
  rename(
    Name = Species
  )

df <- df[ , !(names(df) %in% drops)]

files <- files %>% 
  rename(
   Species = x
  )


files <- files %>% 
  rename(
    SpeciesNumber = X
  )

df4 <- df4 %>% 
  rename(
    Name = spp_species
  )

#merge
total2 <- merge(total,df4,by="Name")

df4final<- read.csv('../Results/Formated/FullListSpeciesFinal.csv', stringsAsFactors = FALSE)

total <- merge(fullset, df4final, by = 'Species')


write.csv(total, "../Results/Formated/5Total_processed_Log_RVK_3_models.csv")
