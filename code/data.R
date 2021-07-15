
rm(list = ls())

# require(dplyr)
library(ggplot2)
library(readr) # for reading several csv files in one time

# read the data
Data <- read.csv("../data/LogisticGrowthData.csv")
Meta <- read.csv("../data/LogisticGrowthMetaData.csv")
# taxo <- read.csv('../data/5total_processed.csv')

# delete duplicate rows
#data <- Data %>% dplyr::distinct(Time, PopBio, Temp, Time_units, PopBio_units, Species, Medium, Rep, Citation)

# insert ID column
Data$ID <- paste(Data$Species,"_",Data$Temp,"_",Data$Medium,"_",Data$Citation)
# delete the negative population
Data <- subset(Data, Data$PopBio > 0)
# log the population
Data$logN <- log(Data$PopBio)

# plot
# for (i in 1:length(unique(Data$ID))){
#   id <- unique(Data$ID)[i]
#   data <- subset(Data, Data$ID == id)

# plot by each ID
# FileName <- paste('../results/RawDataPlot/plot_ID_',i)
# pdf(file = FileName)
# print(
#   ggplot(data, aes(x = Time, y = logN)) +
#     geom_point(size = 3) +
#     labs(x = "Time (Hours)", y = "logarithum of Population size")
# )
# graphics.off()
#}


# save pop
write.csv(Data, "../data/pop.csv")

