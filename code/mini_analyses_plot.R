

rm(list = ls())
graphics.off()

setwd("~/Documents/Project/code")

require(ggplot2)

# 1. read the data, starting value and compare models
infos <- read.csv("../data/fit.id.info.csv")
plot.points <- read.csv("../data/fit.id.plot.csv")

Data <- read.csv('../data/mini_pop.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]

################
# plot fitting #
################
IDs <- unique(Data$id)
comp_df <- data.frame()
for (i in 1:length(IDs)){
  idname <- IDs[i]
  data <- Data[Data$id == idname,]
  plot.df <- plot.points[plot.points$id == idname, ]
  info.df <- infos[infos$id == idname, ]
  
  if (is.na(unique(plot.df$plot.point)[1])){
      print(paste("Fit was not successful for ID:",idname))
  }else{
    models <- unique(info.df$model)
    tit <- c()
    
    for (j in length(models)) {
      model <- models[j]
      # get plot notation
      tit.model <- paste0(model, ': AICc = ', info.df[info.df$model == model, ]$AICc,
                                            ', AIC = ', info.df[info.df$model == model, ]$AIC,
                                            ', BIC = ', info.df[info.df$model == model, ]$BIC,
                                            ', Rsquare = ', info.df[info.df$model == model, ]$rsq)
      tit <- c(tit, tit.model)
    }
    
    # get com_df
    comp.df.id <- data.frame(id = idname, 
                             bet.AICc = info.df[info.df$AICc == min(info.df$AICc),]$model,
                             bet.AIC = info.df[info.df$AIC == min(info.df$AIC),]$model,
                             bet.BIC = info.df[info.df$BIC == min(info.df$BIC),]$model,
                             bet.rsq = info.df[info.df$rsq == min(info.df$rsq),]$model)
    comp_df <- rbind(comp_df, comp.df.id)
    
    # plot fit plot
    title <- paste(tit[1],'\n', tit[2])
    FileName <- paste0("../results/fit/plot_", idname,".png")
    png(file = FileName)
    p <- ggplot(data, aes(x = Time, y = logN)) +
      geom_point(size = 1) +
      labs(x = "Time (h)", y = "Logarithm of the population size (logN)") +
      # ggtitle("Model comparison plot") +
      geom_line(data = plot.df, aes(x = time, y = plot.point, colour = model), size=1) +
      theme(legend.position = 'bottom') +
      annotate('text', label = title, x = min(data$Time), y = min(data$logN), hjust = -.2, vjust = 0) 
      # stat_smooth(method = lm, level = 0.95, aes(colour="Cubic")) +
      # scale_colour_manual(name="Model", values=c("darkblue", "darkred", "darkgreen"))
    print(p)
    graphics.off()
  }
}


####################
# r_max V.S. t_lag #
####################

ggplot(infos, aes(x = rmax, y = tlag, colour = model)) +
  geom_point(size = 1) 
  # stat_smooth(method = lm, level = 0.95) +
  # geom_smooth()


##################################
# plot comparing AIC,AICc,BIC,rsq#
##################################

# save the compareson data
write.csv(comp_df, "../data/compare_criteria.csv")

# visualize comparison
comp.item <- colnames(comp_df)[-1]
for (j in 1:length(comp.item)) {
  filename <- paste0("../results/comparasion_", comp.item[j], ".png")
  png(filename)
  p <- ggplot(comp_df, aes(x = comp_df[, j+1])) +
    geom_bar() +
    theme_bw() +
    labs(x = "Model", y = "Count") +
    # scale_x_discrete(labels=c("Gompertz Model", "Baranyi Model")) +
    geom_text(stat='count', aes(label=..count..), vjust=-1) 
  # ggtitle("Which one is the better model? (by comparing AIC)") +
  # theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))
  print(p)
  graphics.off()
}

  