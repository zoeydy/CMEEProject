## Script: MyBars.R
## Author: Sam Turner sat19@ic.ac.uk
## About: Fits all models and saves the model fit objects to ../data. Makes plots demonstrating model fit in log and linear space. 


# clear environment

rm(list=ls())

############
# PACKAGES #
############



required_packages <- c("minpack.lm", "reshape2", "ggplot2", "gridExtra")

for (pkg in required_packages){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

source("model_functions.R")
source("fitting_functions.R")
options(warn=-1);



#############
# LOAD DATA #
#############

# Load data
data <- read.table("../data/LogisticGrowthDataLogClean.csv", sep = "\t", header = T)
# set temperature as factor
data$Temp <- as.factor(data$Temp)

# Get vector of run IDs
IDs = unique(data$ID)

# make dataframe of amount of time removed for each id to make first data at t=0
removed_times <- data.frame(id = IDs, time_removed = rep(NA, length(IDs) ))
for (id in IDs){
  removed_times[id,] <- c(id, min(data[data$ID == id,]$Time))
}

# load initial value df
inits <- read.csv("../data/inits.csv")[,-1]


###################
## MODEL FITTING ##
###################

# fit models, and store model fit objects in list

fit_models <- function(IDs, models){
  fit_list <- list()
  for (id in IDs){
    if (id %% 20 == 0){
      print(paste0(as.character(id), ' / 287 fits completed'))
    }
    fit_list[[id]] <- list()
    fits <- fit_models_multi(id, models)
    for (model in models){
      if (is.list(fits[[model]])){
        fit_list[[id]][[model]] <- fits[[model]]
      }
      else{
        fit_list[[id]][[model]] <- 0
      }
    }
  }
  return(fit_list)
}

print("FITTING MODELS IN LOG SPACE...")
fit_list <- fit_models(IDs, all_models)

print("FITTING MODELS IN LINEAR SPACE...")
fit_list.linear <- fit_models(IDs, all_models.linear)

# save model fits
saveRDS(fit_list, file = "../data/fit_list.rds")
saveRDS(fit_list.linear, file = "../data/fit_list_linear.rds")

###################################
# demonstrations of model fitting #
###################################

print("Producing model fitting demonstration figures.")

# linear vs log residuals
g <- plot_fit_multi.compare(186, c("Logistic", "Gompertz"))
ggsave("../results/compare_log_lin_fit.pdf", plot = g, width = 10, height = 5)


# initial values
pdf("../results/initial_vals_check.pdf")
check_inits(182, non_linear_models)
graphics.off()



# Model fits in log space
g1 <- plot_fit_multi(140, all_models)  + theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal")
g2 <- plot_fit_multi(60, all_models)
g3 <- plot_fit_multi(132, all_models)
g4 <- plot_fit_multi(1, all_models)

mylegend <- make.legend(g1) 
total_plot <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                       g2 + theme(legend.position="none"),
                                       g3 + theme(legend.position="none"),
                                       g4 + theme(legend.position="none"),
                                       nrow=2),
                           mylegend, nrow=2,heights=c(10,1))

ggsave("../results/model_fits.pdf", device = "pdf",plot = total_plot, width =10, height = 9)


# Model fits in linear space
g1.lin <- plot_fit_multi(140, all_models.linear, T)  + theme(legend.direction = "horizontal", legend.position = "bottom",legend.box = "horizontal")
g2.lin <- plot_fit_multi(60, all_models, T)
g3.lin <- plot_fit_multi(132, all_models, T)
g4.lin <- plot_fit_multi(1, all_models, T)

mylegend <- make.legend(g1) 
total_plot <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                       g2 + theme(legend.position="none"),
                                       g3 + theme(legend.position="none"),
                                       g4 + theme(legend.position="none"),
                                       nrow=2),
                           mylegend, nrow=2,heights=c(10,1))

ggsave("../results/model_fits_linear.pdf", device = "pdf",plot = total_plot, width =10, height = 9)


