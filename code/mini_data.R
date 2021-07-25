rm(list = ls())
setwd("~/Documents/Project/code")
require(dplyr)

# read the data
Data <- read.csv("../data/growth_rate_data.csv")

# insert ID column
Data$ID <- paste0(Data$Rep,"_",Data$Species,"_",Data$Temp,"_",Data$Medium,"_",Data$Citation)
length(unique(Data$ID))

# delete duplicate rows
Data <- Data %>% dplyr::distinct(n,Time, PopBio, Temp, Time_units, PopBio_unts, Species, Medium, Rep, Citation,id,ID)
length(unique(Data$ID))

# delete data set has less than 6 points
Data <- subset(Data, Data$n > 5)
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

# # rearrange the id
# Data$id <- rep(NA, nrow(Data))
# for (j in 1:length(unique(Data$ID))) {
#   idname <- unique(Data$ID)[j]
#   Data[Data$ID == idname,]$id <- j
# }


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
write.csv(Data, "../data/mini_pop.csv")
