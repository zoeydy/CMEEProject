rm(list = ls())
graphics.off()

setwd("~/Documents/Project/code")

############
# PACKAGES #
############
required_packages <- c("minpack.lm", "reshape2", "ggplot2", "gridExtra")
for (pkg in required_packages){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
#############
# LOAD DATA #
#############
Data <- read.csv("../data/mini_pop.csv")
Data <- Data[order(Data[,'id'], Data[,'Time']),]
############
# FUNCTION #
############
source("mini_model_fit_fun.R")


# Get vector of run IDs
IDs = unique(Data$id)

# fitting the model
fit.id.info <- data.frame()
fit.id.plot <- data.frame()
for (i in 1:length(IDs)) {
  
  if (i %% 20 == 0){
    print(paste0(i, ' / ', length(IDs),'fits completed'))
  }
  
  idname <- IDs[i]
  data <- Data[Data$id == idname,]
  # get inits
  samp.time = 1000
  samp.fact = 1
  starts <- sample_inits(data,samp.time,samp.fact)
  
  fit.df <- data.frame()
  plot.df <- data.frame()

  models <- fun_models
  fit.id.df <- data.frame()
  plot.id.df <- data.frame()
  for (j in 1:length(models)) {
    
    model.name <- models[j]
    
    fmla <- as.formula(paste0("logN ~", model.name, "(t=Time, r_max, N_0, K, t_lag)"))
    id.model.list <- fit_id_model(id = idname, model = model.name, fmla = fmla)
    df.id.model <- as.data.frame(id.model.list[[1]][1:8])
    
    if (nrow(df.id.model) != 0) {
      # is.list(id.model.list)
      names(df.id.model) <- c('N0', 'Nmax', 'tlag', 'rmax', 'AIC', 'AICc', 'BIC', 'rsq')
      df.id.model$id <- idname
      df.id.model$model <- model.name
      
      plot.id.model <- as.data.frame(id.model.list[[1]][9])
      names(plot.id.model) <- 'plot.point'
      plot.id.model$time <- seq(min(data$Time), max(data$Time), length.out = samp.time)
      plot.id.model$id <- idname
      plot.id.model$model <- model.name
    } else{
      df.id.model <- data.frame(N0 = NA, Nmax = NA, tlag = NA, rmax = NA, AIC = NA, AICc = NA, BIC = NA, rsq = NA, id = idname, model = model.name)
      plot.id.model <- data.frame(plot.point = NA, time = NA, id = idname, model = model.name)
    }
    
    fit.id.df <- rbind(fit.id.df, df.id.model)
    plot.id.df <- rbind(plot.id.df, plot.id.model)
    # fit.df <- rbind(fit.df, fit.id.df)
    # plot.df <- rbind(plot.df, fit.id.plot)
  }

  fit.id.info <- rbind(fit.id.info, fit.id.df)
  fit.id.plot <- rbind(fit.id.plot, plot.id.df)
}


write.csv(fit.id.info, "../data/fit.id.info.csv")
write.csv(fit.id.plot, "../data/fit.id.plot.csv")

