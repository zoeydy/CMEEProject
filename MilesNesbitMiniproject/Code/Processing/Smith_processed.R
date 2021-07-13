library(tidyverse)
library(ggplot2)
library(plyr)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")

df <- read.csv('../Data/Smith/combined_OD_data.csv',
               stringsAsFactors = F)

#change to correct set
drops <- c("ID")
final_df <- df[ , !(names(df) %in% drops)]
final_df <- rename(final_df, c("trait_name" = "PopBio_unts", "bacterial_genus" = "Species", "minute" = "Time", "Temperature" = "Temp",
                               "replicate" = "Rep", "trait_value" = "PopBio"))

head(final_df)
q <- final_df[,c('Species', 'Rep', 'Time', 'PopBio_unts', 'PopBio', 'Temp')]
w <- final_df[,c('Time','PopBio',"PopBio_unts", 'Temp','Species','Rep')]

w['Medium'] <- '10% LB broth'
w['Time_units'] <- 'Hours'
w['Citation'] <- 'Smith_et_al'

w <- w[,c('Time','PopBio','Temp','Time_units', 'PopBio_unts','Species','Medium', 'Rep', 'Citation')]


#plot
w %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#save
write.csv(w,"../Results/Formated/Smith_processed.csv")

