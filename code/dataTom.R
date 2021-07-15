
rm(list = ls())
graphics.off()

Data <- read.csv('../data/growth_rate_data.csv')
Meta <- read_file("../data/growth_rate_meta_data")
# unique(Data$Citation)
# str(Data)

# add ID column
Data$ID <- paste0(Data$Species,"_",Data$Temp,"_",Data$Medium,"_",Data$Citation)

# plot the raw data
for (i in 1:length(unique(Data$ID))) {
  
  idname <- unique(Data$ID)[i]
  data <- subset(Data, Data$ID == idname)
  # plot
  fileName <- paste0('../results/RawDataPlot/',i)
  pdf(fileName)
  p <- ggplot(data , aes(x = Time, y = PopBio)) +
    geom_point(size = 2) +
    xlab('Time (Hours)') +
    ylab('Population Size')
  print(p)
  graphics.off()
  
}
