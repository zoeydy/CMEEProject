
require(ggplot2)


StartGom <- read.csv("../data/mini_gomp_Start_value.csv")
StartBar <- read.csv("../data/mini_baranyi_Starting_Value.csv")

PlotGom <- read.csv("../data/mini_gomp_plot_points.csv")
PlotBar <- read.csv("../data/mini_baranyi_plot_points.csv")

Data <- read.csv('../data/mini_pop.csv')
Data <- Data[order(Data[,'ID'], Data[,'Time']),]

start_gom_na <- subset(StartGom, is.na(StartGom$AIC))
start_bar_na <- subset(StartBar, is.na(StartBar$AIC))
len.id.gom <- length(unique(start_gom_na$ID))
len.id.bar <- length(unique(start_bar_na$ID))

if (len.id.gom > len.id.bar) {
  len.id <- len.id.gom
  start_na <- start_gom_na
} else {
  len.id <- len.id.bar
  start_na <- start_bar_na
}

for (i in 1:len.id) {
  
  idname <- unique(start_na$ID)[i]
  dat <- Data[Data$ID == idname,]
  
  filepath <- paste0("../results/CheckData/",i,".png")
  png(filepath)
  
  ggplot(dat, aes(x = Time, y = logN)) +
    geom_point(size = 2) +
    labs(x = "Time (Hours)", y = "logarithum of Population size")
  
  ggsave(
        filename = filepath,
        width=10,
        height=10
      )
  
  graphics.off()
}
