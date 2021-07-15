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
files <- read.csv('../Results/Formated/SpeciesList.csv', stringsAsFactors = FALSE)

df <- read.csv('../Results/Formated/FullListSpecies.csv', stringsAsFactors = FALSE)


#clean df
df <- df %>% 
  rename(
   Name = Species
  )
   

df$spp_order <- NA
   df$spp_genus <- NA
   df$spp_family <- NA
   df$spp_species <- NA
   df$spp_class <- NA

   
drops <- c("X")
df <- df[ , !(names(df) %in% drops)]


#make second df
df2<- data_frame('NA')   
df2$spp_order <- NA
df2$spp_genus <- NA
df2$spp_family <- NA
df2$spp_species <- NA
df2$spp_class <- NA

df2 <- matrix(NA,80,5)

# # trial
# a1 <- gnr_resolve(unique(df$Species)[1])
# b1 <- classification(a1$matched_name[1], db = "ncbi")
# c1 <- b1$`Chryseobacterium balustinum`$name
# d1 <- b1$`Chryseobacterium balustinum`$rank
# e1 <- rbind(d1, c1)
# colnames(e1) <- as.character(d1)
# e1 = e1[-1]
# 
# 
# a <- 'd5c0045aff4bbe44fb5be7b6d8deaeb60109'
# 
# x <- data.frame()
# for(sppname in unique(df$Species)) {
#   tryCatch({
#     
#     # taxize::use_entrez('d5c0045aff4bbe44fb5be7b6d8deaeb60109')
#     c <- b[[1]]
#     d <- c[7:10, ]
#     test <- rbind(test, d)
#     print(paste("Value of d:", d))
#   }
#   )
# }

b <- classification(unique(df$Name), db = "ncbi")


# b <- merge(b,df, by ='Name')
# b %>% filter(!is.na(Value))


for (i in seq(1, length(b))){
  
# for (i in seq(1, 5)){
  working <- b[[i]]
  

if(!is.na(working))tryCatch({
  spp_species <- filter(working, rank=="species")$name
 }, error=function(e){spp_species <- get0("spp_species", ifnotfound = 'eh')})
if(!is.na(working)){tryCatch({
  spp_genus <- filter(working, rank=="genus")$name
 }, error=function(e){spp_genus <- get0("spp_genus", ifnotfound = 'eh')}) 
if(!is.na(working))tryCatch({
  spp_family <- filter(working, rank=="family")$name
 }, error=function(e){spp_family <- get0("spp_family", ifnotfound = 'eh')})
if(!is.na(working))tryCatch({
      spp_order <- filter(working, rank=="order")$name
    }, error=function(e){spp_order <- get0("spp_order", ifnotfound = 'eh')})
if(!is.na(working))tryCatch({
  spp_class <- filter(working, rank=="class")$name
 }, error=function(e){spp_class <- get0("spp_class", ifnotfound = 'eh')})



  
  spp_class <- get0("spp_class", ifnotfound = 'NA')
  spp_order <- get0("spp_order", ifnotfound = 'NA')
  spp_family <- get0("spp_family", ifnotfound = 'NA')
  spp_species <- get0("spp_species", ifnotfound = 'NA')
  
  z <- cbind(spp_class,spp_family,spp_genus,spp_order,spp_species)
  df2 <- rbind(df2, z)
#})
}
}


files2 <- read.csv("../Results/Formated/4Total_processed_Log_RVK_3_models.csv")
# 
# 
# 
# write.csv(df2, '../Results/Formated/FullListSpeciesStart2.csv')



df3 <- read.csv('../Results/Formated/FullListSpeciesStart3.csv', stringsAsFactors = FALSE)
