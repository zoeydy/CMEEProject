rm(list = ls())
setwd("~/Documents/CMEEProject/code")
require(dplyr)
require(ggplot2)

# read the data
Data <- read.csv("../data/growth_rate_data.csv")

# insert ID column
Data$ID <- paste0(Data$Rep,"_",Data$Species,"_",Data$Temp,"_",Data$Medium,"_",Data$Citation)
length(unique(Data$ID))

# delete problematic data sets (checking species)
spe <- unique(Data$Species) # prob: c(105, 61, 60, 55, 51, 46, 41, 32)
for (i in c(105, 61, 60, 55, 51, 46, 41, 32)) {
  Data <- subset(Data, Data$Species != spe[i])
}
length(unique(Data$ID))

# delete duplicate rows
Data <- Data %>% dplyr::distinct(n,Time, PopBio, Temp, Time_units, PopBio_unts, Species, Medium, Rep, Citation,id,ID)
length(unique(Data$ID))

# delete the data set with negative population size
neg.df <- subset(Data, Data$PopBio < 0)
neg.id <- unique(neg.df$ID)
for (i in 1:length(neg.id)) {
  Data <- subset(Data, Data$ID != neg.id[i])
}
length(unique(Data$ID))

# log the population
Data <- Data[Data$PopBio > 0, ]
Data$logN <- log(Data$PopBio)
length(unique(Data$id))

# plot
for (i in 1:length(unique(Data$id))){
  idname <- unique(Data$id)[i]
  data <- subset(Data, Data$id == idname)
  Data[Data$id == idname, ]$n <- nrow(data)

  # # plot by each ID
  # FileName <- paste('../results/RawDataPlot/plot_ID_',i)
  # png(file = FileName)
  # print(
  #   ggplot(data, aes(x = Time, y = logN)) +
  #     geom_point(size = 3) +
  #     labs(x = "Time (Hours)", y = "logarithum of Population size")
  # )
  # graphics.off()
}

# delete data set has less than 6 points
Data <- Data[Data$n > 5,]
length(unique(Data$ID))

# add taxonomy information


# save pop
write.csv(Data, "../data/mini_pop.csv")
