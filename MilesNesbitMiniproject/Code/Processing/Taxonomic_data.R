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
files <- read.csv('../Results/Formated/Total_processed.csv', stringsAsFactors = FALSE)

#clean
drops <- c("X", "X.1")
files <- files[ , !(names(files) %in% drops)]

#add empty column to df

df <- files

files <- read.csv('../Results/Formated/Total_processed.csv')
files <- files[ , !(names(files) %in% drops)]


#remove '.'
df$Species <- sub('\\.', ' ', df$Species)


#1
df$Species[1]

gnr_resolve("Chryseobacterium balustinum")

classification("Chryseobacterium balustinum", db = "ncbi")
df_subset <- subset(df , df$Species) [1]



#
a1 <- gnr_resolve(unique(df$Species)[1])
b1 <- classification(a1$matched_name[1], db = "ncbi")
c1 <- b1$`Chryseobacterium balustinum`$name
d1 <- b1$`Chryseobacterium balustinum`$rank
e1 <- rbind(d1, c1)
colnames(e1) <- as.character(d1)
e1 = e1[-1]


#loop it up
# for(i in 1:105) {  
#   tryCatch({
#      
#     a <- gnr_resolve(unique(df$Species)[i])
#     b <- classification(a$matched_name[1], db = "ncbi")
#     c <- b$`Chryseobacterium balustinum`$name
#     
#     
#     
#   }
# )
# }

test <- data.frame()
for(i in unique(df$Species)) {
  tryCatch({
    a <- gnr_resolve(i)
    b <- classification(a$matched_name[1], db = "ncbi")
    c <- b[[1]]
    d <- c[7:10, ]
    test <- rbind(test, d)
    print(paste("Value of d:", d))
  }
)
}
