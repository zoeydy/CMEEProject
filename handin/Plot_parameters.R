rm(list = ls())
graphics.off()

require(ggplot2)

# 1. read the data, starting value (parameters)
StartGom <- read.csv("../result/gompertz_Starting_Value.csv")

Data <- read.csv('../data/pop.csv')
Data <- Data[order(Data[,'ID'], Data[,'Time']),]

# 2. plot parameters for different species
for (i in 1:length(unique(Data$Species))) {
  #browser()
  data <- StartGom[StartGom$Species == unique(Data$Species)[i],]
  title1 <- paste("The r_max of", unique(Data$Species)[i])
  title2 <- paste("The K of", unique(Data$Species)[i])
  title3 <- paste("The t_lag of", unique(Data$Species)[i])
  
  FileName <- paste0("../result/hist_para/Parameters_", i,".png")
  png(file = FileName)
  
  par(mfcol=c(3,1)) 
  par(mfg = c(1,1)) 
  hist(data$r_max,main = title1, xlab="r_max value",ylab="Count")
  
  par(mfg = c(2,1)) 
  hist(data$K, main = title2, xlab="K value",ylab="Count")
  
  par(mfg = c(3,1))
  hist(data$t_lag, main = title3, xlab="t_lag value",ylab="Count")

  graphics.off()
}