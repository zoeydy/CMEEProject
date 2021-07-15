library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)

#general start
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
df <- read.csv('../../../Documents/ggplot2.csv')


#change to long format
data_long <- gather(df, Image_Type, Median_NDVI, Coastal.Rocks:Hardwoods, factor_key=TRUE)
data_long
#plot
ggplot(data_long, aes(x= Month, y= Median_NDVI, group=Image_Type, colour=Image_Type)) + 
  geom_line(size=.75) + geom_point() +
  scale_x_discrete(limits=c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                            "Aug","Sep","Oct","Nov","Dec"))

