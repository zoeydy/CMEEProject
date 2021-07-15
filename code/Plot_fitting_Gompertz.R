

rm(list = ls())
graphics.off()

require(ggplot2)

# 1. read the data, starting value and compare models
StartGom <- read.csv("../data/gompertz_Starting_Value.csv")

PlotGom <- read.csv("../data/gompertz_plot_points.csv")

Data <- read.csv('../data/pop.csv')
Data <- Data[order(Data[,'ID'], Data[,'Time']),]

# 2. plot & get AIC of qubic model

for (i in 1:285){
  #browser()
  plot_df <- data.frame()
  id <- unique(Data$ID)[i]
  data <- Data[Data$ID == id,]
  # plot
  time2plot <- seq(min(data$Time), max(data$Time), length=1000)
  
  if (is.na(subset(StartGom, StartGom$ID==id)$AICc)){
    next                   
  }else{
    plot_gom <- data.frame(time = time2plot,logN = subset(PlotGom, PlotGom$ID==id)$pred_gom, model = rep("gompertz",1000))
    plot_df <- rbind(plot_df, plot_gom)
  }
  
  FileName <- paste0("../results/FitPlotGomp/plot_", i,".png")
  png(file = FileName)
  title <- paste("Fitting Gompertz Model plot, ID:",id)
  p <- ggplot(data, aes(x=Time, y=logN)) +
    geom_point() +
    labs(x = "Time (h)", y = "Logarithm of the population size (logN)") +
    ggtitle(title) +
    # annotate(geom = 'text', x = Inf, y = -Inf, hjust = 0, vjust = -0.2) +
    # ,label = comp_lable
    # stat_smooth(method = lm, level = 0.95, aes(colour="Cubic")) +
    geom_line(data = plot_df ,aes(x = time, y = logN, color = model), size=1)
    # scale_colour_manual(name="Model", values=c("darkblue", "darkred", "darkgreen"))
  print(p)
  graphics.off()
}

