rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/sandbox/")

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
Data <- read.csv("../data/pop.csv")
Data <- Data[order(Data[,'id'], Data[,'Time']),]
############
# FUNCTION #
############
source("model_fit_fun.R")


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
  
  fmla <- as.formula(paste0("logN ~", " gompertz_model", "(t=Time, r_max, N_0, K, t_lag)"))
  id.model.list <- fit_id_model(id = idname, model = "gompertz_model", fmla = fmla)
  df.id.model <- as.data.frame(id.model.list[[1]][1:8])
  
  if (nrow(df.id.model) != 0) {
    # is.list(id.model.list)
    names(df.id.model) <- c('N0', 'Nmax', 'tlag', 'rmax', 'AIC', 'AICc', 'BIC', 'rsq')
    df.id.model$id <- idname
    df.id.model$model <- "gompertz_model"
    # get info about temperature, species and medium
    df.id.model$temp_c <- data$Temp[1]
    df.id.model$species <- data$Species[1]
    df.id.model$medium <- data$Medium[1]
    # add kelvin temperature
    df.id.model$temp_k <- data$Temp[1] + 273
    
    plot.id.model <- as.data.frame(id.model.list[[1]][9])
    names(plot.id.model) <- 'plot.point'
    plot.id.model$time <- seq(min(data$Time), max(data$Time), length.out = samp.time)
    plot.id.model$id <- idname
    plot.id.model$model <- "gompertz_model"
  } else{
    df.id.model <- data.frame(N0 = NA, Nmax = NA, tlag = NA, rmax = NA, AIC = NA, AICc = NA, BIC = NA, rsq = NA, 
                              id = idname, model = "gompertz_model", temp_c=data$Temp[1], species=data$Species[1], medium=data$Medium[1], temp_k=data$Temp[1] + 273)
    plot.id.model <- data.frame(plot.point = NA, time = NA, id = idname, model = "gompertz_model")
  }
  
  fit.id.info <- rbind(fit.id.info, df.id.model)
  fit.id.plot <- rbind(fit.id.plot, plot.id.model)
}


write.csv(fit.id.info, "../sandbox/gomp.info.csv")
write.csv(fit.id.plot, "../sandbox/gomp.plot.csv")

