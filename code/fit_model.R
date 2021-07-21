

# clear environment
rm(list=ls())

############
# PACKAGES #
############
required_packages <- c("minpack.lm", "reshape2", "ggplot2", "gridExtra")
for (pkg in required_packages){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
source("model&fitting_functions.R")
options(warn=-1);

#############
# LOAD DATA #
#############
data <- read.csv("../data/processed_data.csv")
# set temperature as factor
data$Temp <- as.factor(data$Temp)

# Get vector of run IDs
IDs = unique(data$id)

# make dataframe of amount of time removed for each id to make first data at t=0
removed_times <- data.frame(id = IDs, time_removed = rep(NA, length(IDs) ))
for (id in IDs){
  removed_times[id,] <- c(id, min(data[data$id == id,]$Time))
}

# load starting value data frame
inits <- read.csv("../data/starting_value.csv")[,-1]


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

# ###################################
# # demonstrations of model fitting #
# ###################################
# 
# print("Producing model fitting demonstration figures.")
# 
# # linear vs log residuals
# g <- plot_fit_multi.compare(186, c("Logistic", "Gompertz"))
# ggsave("../results/compare_log_lin_fit.pdf", plot = g, width = 10, height = 5)
# 
# 
# # initial values
# pdf("../results/initial_vals_check.pdf")
# check_inits(182, non_linear_models)
# graphics.off()
# 
# 
# 
# # Model fits in log space
# g1 <- plot_fit_multi(140, all_models)  + theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal")
# g2 <- plot_fit_multi(60, all_models)
# g3 <- plot_fit_multi(132, all_models)
# g4 <- plot_fit_multi(1, all_models)
# 
# mylegend <- make.legend(g1) 
# total_plot <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
#                                        g2 + theme(legend.position="none"),
#                                        g3 + theme(legend.position="none"),
#                                        g4 + theme(legend.position="none"),
#                                        nrow=2),
#                            mylegend, nrow=2,heights=c(10,1))
# 
# ggsave("../results/model_fits.pdf", device = "pdf",plot = total_plot, width =10, height = 9)
# 
# 
# # Model fits in linear space
# g1.lin <- plot_fit_multi(140, all_models.linear, T)  + theme(legend.direction = "horizontal", legend.position = "bottom",legend.box = "horizontal")
# g2.lin <- plot_fit_multi(60, all_models, T)
# g3.lin <- plot_fit_multi(132, all_models, T)
# g4.lin <- plot_fit_multi(1, all_models, T)
# 
# mylegend <- make.legend(g1) 
# total_plot <- grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
#                                        g2 + theme(legend.position="none"),
#                                        g3 + theme(legend.position="none"),
#                                        g4 + theme(legend.position="none"),
#                                        nrow=2),
#                            mylegend, nrow=2,heights=c(10,1))
# 
# ggsave("../results/model_fits_linear.pdf", device = "pdf",plot = total_plot, width =10, height = 9)
# 



























# ########
# # mine #
# ########
# 
# #################
# # MODEL FITTING #
# #################
# # fit the model, generate the data frame for plotting and save the results
# starts <- start_value_df[start_value_df$id == id, ]
# # sample for getting best starting value
# SampTime <- 1000
# fact <- 2
# SampMin <- runif(SampTime,N0start-abs(N0start)*fact,N0start+abs(N0start)*fact)
# SampMax <- runif(SampTime,Kstart-abs(Kstart)*fact,Kstart+abs(Kstart)*fact)
# Samp_t_lag <- runif(SampTime, tlagStart-abs(tlagStart)*fact, tlagStart+abs(tlagStart)*fact)
# Samp_r_max <- runif(SampTime, 0 , rStart + abs(rStart)*fact)
# # fit Gompertz model
# gomp_plot_df <- data.frame()
# comp <- list()
# for (i in 1:length(unique(Data$ID))){
#   #browser()
#   # subset the Data by ID 
#   idname <- unique(Data$ID)[i]
#   dat <- subset(Data, Data$ID == idname)
#   
#   # for loop in each ID for finding starting value 
#   start_id_list <- list()
#   plot_id_list <- list()
#   for (j in 1:SampTime){
#     result <- GuessStart(dat,t = Time, r_maxS = Samp_r_max[j], N_0S = SampMin[j], KS = SampMax[j],t_lagS = Samp_t_lag[j])
#     start_id_list[[j]] <- result[1:8]
#     plot_id_list[[j]] <- result[9]
#   }
#   
#   # write the start_id_list into data frame
#   n.obs <- sapply(start_id_list, length)
#   seq.max <- seq_len(max(n.obs))
#   start_id_mat <- t(sapply(start_id_list, "[", i = seq.max))
#   start_id_df <- as.data.frame(start_id_mat)
#   start_id_df$index <- seq(1:1000)
#   colnames(start_id_df) <- c("AIC","AICc","BIC","RSq","r_max","N_0","K","t_lag","index")
#   index_AIC <- start_id_df$index[which.min(start_id_df$AIC)]
#   index_AICc <- start_id_df$index[which.min(start_id_df$AICc)]
#   index_BIC <- start_id_df$index[which.min(start_id_df$BIC)]
#   comp[[i]] <- all(sapply(list(index_AIC,index_BIC,index_AICc), function(x) x == index_AICc))
#   
#   # write the plot_id_list into data frame
#   plot_id_df <- as.data.frame(plot_id_list[index_AICc])
#   
#   if (nrow(start_id_df[index_AICc,]) == 0) {
#     start_id_df <- data.frame(ID=idname,Model="gompertz",AIC=NA,AICc=NA,BIC=NA,RSq=NA,r_max=NA,N_0=NA,K=NA,t_lag=NA,index=NA)
#     start_value_df <- rbind(start_value_df, start_id_df)
#     
#     plot_id_df <- data.frame(pred_gom=NA,ID=idname,model="gompertz")
#     gomp_plot_df <- rbind(gomp_plot_df, plot_id_df)
#   }else{
#     start_id_df <- data.frame(ID=idname,model="gompertz",start_id_df[index_AICc,])
#     colnames(start_id_df) <- c("ID","Model","AIC","AICc","BIC","RSq","r_max","N_0","K","t_lag","index")
#     start_value_df <- rbind(start_value_df, start_id_df)
#     
#     plot_id_df$ID <- idname
#     plot_id_df$model <- "gompertz"
#     colnames(plot_id_df) <- c("pred_gom", "ID", "model")
#     gomp_plot_df <- rbind(gomp_plot_df, plot_id_df)
#   }
# }
# # fit Baranyi model
# 
# 
# 
# write.csv(gomp_plot_df, "../results/plot_points.csv")