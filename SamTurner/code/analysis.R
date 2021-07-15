## Script: analysis.R
## Author: Sam Turner sat19@ic.ac.uk

## About: loads and analyses model fits. Produces figures for write up, and runs statistical models.


# clear environment

rm(list=ls())

############
# PACKAGES #
############

required_packages <- c("matrixStats",
                       "dplyr",
                       "minpack.lm",
                       "ggplot2",
                       "lme4",
                       "gridExtra",
                       "toOrdinal",
                       "reshape2"
                       )

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


# get all unique combinations of factors (species, temp, medium)
combos <- unique(cbind(as.vector(data$Species), as.vector(data$Temp),as.vector(data$Medium),as.vector(data$Rep)))
for (i in 1:dim(combos)[1]){
  species <- combos[i,1]
  temp <- combos[i,2]
  medium <- combos[i,3]
  repli <- combos[i,4]
}

# make dataframe of amount of time removed for each id to make first data at t=0
removed_times <- data.frame(id = IDs, time_removed = rep(NA, length(IDs) ))
for (id in IDs){
  removed_times[id,] <- c(id, min(data[data$ID == id,]$Time))
}


###################
# LOAD MODEL FITS #
###################

fit_list <- readRDS("../data/fit_list.rds")
fit_list.linear <- readRDS("../data/fit_list_linear.rds")

print("ANALYSING RESULTS...")
# visualise model fit success
get_fit_successes <- function(fit_list){

  fit_success_df <- array(NA, dim = c(length(IDs), 1+length(names(fit_list[[1]])) ) )
  fit_success_df[,1] <- IDs
  fit_success_df <- data.frame(fit_success_df)
  colnames(fit_success_df) <- c("IDs", names(fit_list[[1]]))
  
  for (id in IDs){
    for (model in names(fit_list[[1]])){
        if (is.list(fit_list[[id]][[model]])) fit_success_df[[id,model]] <- "FIT SUCCESS"
        else fit_success_df[[id,model]] <- "FIT FAIL"
    }
  }
  return(fit_success_df)
}

fit_success_df = get_fit_successes(fit_list)
fit_success_df.linear = get_fit_successes(fit_list.linear)

# successes in log space
long <- melt(fit_success_df, id.var = 'IDs')
ggplot(long, aes(x=variable, y=IDs)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_manual(values=c("red", "blue", "black")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# successes in linear space
long <- melt(fit_success_df.linear, id.var = 'IDs')
ggplot(long, aes(x=variable, y=IDs)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_manual(values=c("red", "blue", "black")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))



########################
# AICs and weights
########################

# calculate RSS for each model fit

get_RSS <- function(fit_list){
  RSS_df <- array(NA, dim = c(length(IDs), length(names(fit_list[[1]])) ) )
  RSS_df <- data.frame(RSS_df)
  colnames(RSS_df) <-  names(fit_list[[1]])
  
  for (id in IDs){
    for (model in names(fit_list[[1]])){
      if (is.list(fit_list[[id]][[model]])){
        RSS_df[[id,model]] <- sum(resid(fit_list[[id]][[model]])^2)
      }
      else{
        RSS_df[[id,model]] <- NA
      }
    }
  }
  return(RSS_df)
}

RSS_df <- get_RSS(fit_list)
RSS_df.linear <- get_RSS(fit_list.linear)





########################
# Remove bad data sets #
########################

# check which data sets are bad by visualising worst R2 values

# calculate R2 for each model fit
R2_df = data.frame("Logistic"=rep(NA,length(IDs)),"Gompertz"=rep(NA,length(IDs)),"Baranyi"=rep(NA,length(IDs)),"Buchanan"=rep(NA,length(IDs)),"Quadratic"=rep(NA,length(IDs)), "Cubic"=rep(NA,length(IDs)) )    

for (id in IDs){
  for (model in all_models){
    if (fit_success_df[id,model] == "FIT SUCCESS"){
      R2_df[id,model] <- 1 - RSS_df[id,model]  / sum((data[data$ID == id,]$logPopBio - mean(data[data$ID == id,]$logPopBio))^2)
    }
  }
}


worst_R2 <- order(rowMaxs(as.matrix(R2_df), na.rm = T))

# these are ids corresponding to data sets which do not show a meaningful times seris of population sizes, so hsould be removed from analysis.
bad_data <- c(243,21,17,231,227,5,248,20,222,14,109,110,232,105,104,111,113,106,34,275)

pdf("../results/removed_ids.pdf")
par(mfrow = c(8,8), mai = c(0.05,0.05,0,0))
for (i in 1:64){
  id <- worst_R2[[i]]
  tran <- c(min(data[data$ID == id,]$Time), max(data[data$ID == id,]$Time))
  popran <- c(min(data[data$ID == id,]$logPopBio), max(data[data$ID == id,]$logPopBio))
  plot( data[data$ID == id,]$Time, data[data$ID == id,]$logPopBio, xaxt='n', yaxt='n', pch = 20, cex = 0.5, xlim = tran + diff(tran)*c(-0.1, 0.1), ylim = popran + diff(popran)*c(-0.1, 0.1))
  if (id %in% bad_data){
    rect(xleft = tran[1]-diff(tran)*0.11, xright= tran[2]+diff(tran)*0.11, ybottom= popran[1]-diff(popran)*0.11, ytop=popran[2]+diff(popran)*0.11, border = "red", lw = 2)
  }
  text(x = max(data[data$ID == id,]$Time)-0.15*(max(data[data$ID == id,]$Time)-min(data[data$ID == id,]$Time)), y = min(data[data$ID == id,]$logPopBio)+0.1*(max(data[data$ID == id,]$logPopBio)-min(data[data$ID == id,]$logPopBio)),  round(max(as.matrix(R2_df)[id,], na.rm = T), 2), )
}
suppress <- dev.off()

IDs.filtered <- IDs[-bad_data]



#############
# AIC & BIC #
#############

# calculate AIC and BIC scores for each model fit

get_AIC <- function(RSS){
  AIC_df <- array(NA, dim = c(length(IDs), length(names(fit_list[[1]])) ) )
  AIC_df <- data.frame(AIC_df)
  colnames(AIC_df) <-  colnames(RSS)
  
  for (id in IDs){
    for (model in colnames(RSS)){
      if (!is.na(RSS[id,model]) ){
        n = dim(data[data$ID == id,])[1]
        k = model_params[[model]]
        AIC_df[[id,model]] <- 2 * k + n * log( RSS[id,model] / n )
      }
      
      else {
        AIC_df[[id,model]] <- NA
      }
    }
  }
  return(AIC_df)
}

get_BIC <- function(RSS){
  BIC_df <- array(NA, dim = c(length(IDs), length(names(fit_list[[1]])) ) )
  BIC_df <- data.frame(BIC_df)
  colnames(BIC_df) <-  colnames(RSS)
  
  for (id in IDs){
    for (model in colnames(RSS)){
      if (!is.na(RSS[id,model]) ){
        n = dim(data[data$ID == id,])[1]
        k = model_params[[model]]
        BIC_df[[id,model]] <- log(n) * k + n * log( RSS[id,model] / n )
      }
      
      else {
        BIC_df[[id,model]] <- NA
      }
    }
  }
  return(BIC_df)
}

AIC_df <- get_AIC(RSS_df)
AIC_df.linear <- get_AIC(RSS_df.linear)

BIC_df <- get_BIC(RSS_df)
BIC_df.linear <- get_BIC(RSS_df.linear)


# filter AICs, BICs, and IDs to remove bad data sets
IDs.filtered <- IDs[-bad_data]
AIC.filtered <- AIC_df[IDs.filtered,]
BIC.filtered <- BIC_df[IDs.filtered,]
AIC.filtered.linear <- AIC_df.linear[IDs.filtered,]
BIC.filtered.linear <- BIC_df.linear[IDs.filtered,]



##################
# Akaike weights #
##################

# calculate delta_AIC scores for each ID
get_deltas <- function(AICs){
  
  deltas_df <- array(NA, dim = c(length(IDs), 3 ) )
  deltas_df <- data.frame(deltas_df)
  colnames(deltas_df) <-  colnames(AICs)[2:4]
  
  for (model in colnames(deltas_df)){
  for (id in IDs){
    if (is.na(AICs[id,model])){deltas_df[[id,model]] <- NA}
    else{
      deltas_df[id,model] <- AICs[id,model] - min(AICs[id,])
    }
  }
  }
  return(deltas_df)
}

get_deltas.full <- function(AICs){
  
  deltas_df <- array(NA, dim = c(length(IDs), 6 ) )
  deltas_df <- data.frame(deltas_df)
  colnames(deltas_df) <-  colnames(AICs)
  
  for (model in colnames(deltas_df)){
    for (id in IDs){
      if (is.na(AICs[id,model])){deltas_df[[id,model]] <- NA}
      else{
        deltas_df[id,model] <- AICs[id,model] - min(AICs[id,])
      }
    }
  }
  return(deltas_df)
}

deltas.full <- get_deltas.full(AIC_df)
deltas.linear.full <- get_deltas.full(AIC_df.linear)

deltas <- get_deltas(AIC_df)
deltas.linear <- get_deltas(AIC_df.linear)

# calculate Akaike weights for each ID

get_weights <- function(deltas){
  weights_df <- array(NA, dim = c(length(IDs), dim(deltas)[2] ) )
  weights_df <- data.frame(weights_df)
  colnames(weights_df) <-  colnames(deltas)
  
  for (id in IDs){
    total <- sum(exp(-deltas[id,] / 2), na.rm = T)
    for (model in colnames(deltas)){
      if (is.null(deltas[id,model])){weights_df[[id,model]] <- 0}
      else{
        weights_df[id,model] <- exp(-deltas[id,model] / 2) / total
      }
    }
  }
  return(weights_df)
  
}

weights.full <- get_weights(deltas.full)
weights.linear.full <- get_weights(deltas.linear.full)

weights <- get_weights(deltas)
weights.linear <- get_weights(deltas.linear)

#######################
# PARAMETER ESTIMATES #
#######################

# Get data frames of parameter estimates for each model, in log space and in linear space:
# Nmax df #

Nmax_df = data.frame("Gompertz"=rep(NA,length(IDs)),"Baranyi"=rep(NA,length(IDs)),"Buchanan"=rep(NA,length(IDs)))    
Nmax_df.linear = data.frame("Gompertz.exp"=rep(NA,length(IDs)),"Baranyi.exp"=rep(NA,length(IDs)),"Buchanan.exp"=rep(NA,length(IDs)))    

for (id in IDs){
  for (model in non_linear_models){
    if (model != "Logistic" & model != "Logistic.exp"){
      if (fit_success_df[id,model] == "FIT SUCCESS"){
        Nmax_df[id,model] <- coef(fit_list[[id]][[model]])[[2]]
        Nmax_df.linear[id,paste0(model,".exp")] <- coef(fit_list.linear[[id]][[paste0(model,".exp")]])[[2]]
      }
      else {
        Nmax_df[[id,model]] <- NA
        Nmax_df.linear[[id,paste0(model,".exp")]] <- NA
      }
    }
  }
}


# Nmin df #

Nmin_df = data.frame("Gompertz"=rep(NA,length(IDs)),"Baranyi"=rep(NA,length(IDs)),"Buchanan"=rep(NA,length(IDs)))  
Nmin_df.linear = data.frame("Gompertz.exp"=rep(NA,length(IDs)),"Baranyi.exp"=rep(NA,length(IDs)),"Buchanan.exp"=rep(NA,length(IDs)))    

for (id in IDs){
  for (model in non_linear_models){
    if (model != "Logistic" & model != "Logistic.exp"){
      if (fit_success_df[id,model] == "FIT SUCCESS"){
        Nmin_df[id,model] <- coef(fit_list[[id]][[model]])[[1]]
        Nmin_df.linear[id,paste0(model,".exp")] <- coef(fit_list.linear[[id]][[paste0(model,".exp")]])[[1]]
      }
      else {
        Nmin_df[[id,model]] <- NA
        Nmin_df.linear[[id,paste0(model,".exp")]] <- NA
      }
    }
  }
}

# tlag df #

tlag_df = data.frame("Gompertz"=rep(NA,length(IDs)),"Baranyi"=rep(NA,length(IDs)),"Buchanan"=rep(NA,length(IDs)))  
tlag_df.linear = data.frame("Gompertz.exp"=rep(NA,length(IDs)),"Baranyi.exp"=rep(NA,length(IDs)),"Buchanan.exp"=rep(NA,length(IDs)))    

for (id in IDs){
  for (model in non_linear_models){
    if (model != "Logistic" & model != "Logistic.exp"){
      if (fit_success_df[id,model] == "FIT SUCCESS"){
        tlag_df[id,model] <- coef(fit_list[[id]][[model]])[[4]] + removed_times[id, "time_removed"]
        tlag_df.linear[id,paste0(model,".exp")] <- coef(fit_list.linear[[id]][[paste0(model,".exp")]])[[4]]  + removed_times[id, "time_removed"]
        
      }
      else {
        tlag_df[[id,model]] <- NA
        tlag_df.linear[[id,paste0(model,".exp")]] <- NA
        
      }
    }
  }
}


# rmax df #

rmax_df = data.frame("Gompertz"=rep(NA,length(IDs)),"Baranyi"=rep(NA,length(IDs)),"Buchanan"=rep(NA,length(IDs)))  
rmax_df.linear = data.frame("Gompertz.exp"=rep(NA,length(IDs)),"Baranyi.exp"=rep(NA,length(IDs)),"Buchanan.exp"=rep(NA,length(IDs)))    

for (id in IDs){
  for (model in non_linear_models){
    if (model != "Logistic" & model != "Logistic.exp"){
      if (fit_success_df[id,model] == "FIT SUCCESS"){
        rmax_df[id,model] <- coef(fit_list[[id]][[model]])[[3]]
        rmax_df.linear[id,paste0(model,".exp")] <- coef(fit_list.linear[[id]][[paste0(model,".exp")]])[[3]]
        
      }
      else {
        rmax_df[[id,model]] <- NA
        rmax_df.linear[[id,paste0(model,".exp")]] <- NA
        
      }
    }
  }
}


#####################################
#1 Akaike weighted parameter values #
#####################################

weighted_Nmax <- weighted_param_est(Nmax_df, weights)
weighted_Nmin <- weighted_param_est(Nmin_df, weights)
weighted_rmax <- weighted_param_est(rmax_df, weights)
weighted_tlag <- weighted_param_est(tlag_df, weights)

weighted_Nmax.linear <- weighted_param_est(Nmax_df.linear, weights.linear)
weighted_Nmin.linear <- weighted_param_est(Nmin_df.linear, weights.linear)
weighted_rmax.linear <- weighted_param_est(rmax_df.linear, weights.linear)
weighted_tlag.linear <- weighted_param_est(tlag_df.linear, weights.linear)

#########################
# Length of growth phase
#########################

# Find the nuber of points in the growth phase which can constrain tlag and rmax estimates.

# Number of points in growth phase
growth_phase_N <- c()

for (id in IDs.filtered){
  # time point early enough
  time_condition.early <- data[data$ID == id,]$Time  > (weighted_tlag[[id]]) *0.9
  time_condition.late <- data[data$ID == id,]$Time  < (weighted_tlag[[id]] + (weighted_Nmax[[id]] - weighted_Nmin[[id]])/weighted_rmax[[id]])*1.1
  
  # both conditions fulfilled
  growth <- sum( time_condition.early & time_condition.late )
  growth_phase_N <- c(growth_phase_N, growth)
}

growth.phase.size.frequencies <- table(growth_phase_N)


# IDs that have enough growth phase points
IDs.filtered.growth <- intersect(IDs.filtered, (1:287)[growth_phase_N > 1 ])
AIC.filtered.growth <- AIC_df[IDs.filtered.growth,]

################################
# Akaike weight distributuions #
################################

weights.full$IDs <- 1:287
weights.full.filtered <- weights.full[IDs.filtered,]

weights.linear.full$IDs <- 1:287
weights.linear.full.filtered <- weights.linear.full[IDs.filtered,]

# LOG SPACE FITS:

# AIC distributions and means plot.
weights_long <- melt(weights.full.filtered, na.rm = T, id.var = "IDs")
weights.means <- weights_long %>%
  group_by(variable) %>%
  summarize(z = mean(value))
g1 <- ggplot(weights_long, aes(x=value )) + geom_density() + theme_bw() + facet_grid(variable ~., scales = "free") + labs(x= "Akaike weight value", y= "Frequency density", title = "A")  +  geom_vline(aes(xintercept = z), weights.means, colour = "red", linetype = 2)


# LINEAR SPACE FITS:
colnames(weights.linear.full.filtered) <- colnames(weights.full.filtered)
weights_long.linear <- melt(weights.linear.full.filtered, na.rm = T, id.var = "IDs")
weights.means.linear <- weights_long.linear %>%
  group_by(variable) %>%
  summarize(z = mean(value))
g2 <- ggplot(weights_long.linear, aes(x=value )) + geom_density() + theme_bw() + facet_grid(variable ~., scales = "free") + labs(x= "Akaike weight value", y= "Frequency density", title = "B") +  geom_vline(aes(xintercept = z), weights.means.linear, colour = "red", linetype = 2)



g.total <- grid.arrange(g1,g2,nrow = 1)
ggsave("../results/Aw_distributions.pdf", plot=g.total, width = 8, height = 5)


# removed for project compilation

# print(weights.means)
# print(weights.means.linear)



######################
# Best fitting model #
######################

IDs.use <- IDs.filtered
AIC.analysis <- AIC_df[IDs.use,]
AIC.analysis.linear <- AIC_df.linear[IDs.use,]


# LOG SPACE

# get ranks of each mdoel for each id
ranks <- cbind(IDs.use,t(apply(AIC.analysis[,], 1, rank, ties.method='max')))

# consider only gompertz so we can compare 3 param vs 4 param
ranks.removed <- cbind(IDs.use,t(apply(AIC.analysis[,-c(2,4)], 1, rank, ties.method='max')))

# how often is each model best - removed for project compilation

# for (col in colnames(ranks)){
#   print(col)
#   print(sum(ranks[,col] == 1))
#   print(dim(ranks)[1])
# }

# best model bar chart
best.df <- data.frame(model = colnames(ranks[,-1]), best = colSums(ranks[,-1]==1))
g1 <- ggplot(data = best.df , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=best), position=position_dodge(width=0.9), vjust=-0.25)+ theme_bw()+ scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Number of best fit", title = "A") + ylim(0,1.1*max(best.df[,2]))

# second best bar
second.df <- data.frame(model = colnames(ranks[,-1]), best = colSums(ranks[,-1]==2))
g2 <- ggplot(data = second.df , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=best), position=position_dodge(width=0.9), vjust=-0.25)+ theme_bw()+ scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Number of second best fit", title = "B") + ylim(0,1.1*max(second.df[,2]))

# mean rank bar chart
best.df.means <- data.frame(model = colnames(ranks[,-1]), best = colMeans(ranks[,-1]))
g3 <- ggplot(data = best.df.means , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=round(best,2)), position=position_dodge(width=0.9), vjust=-0.25) + theme_bw() + scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Mean ranking", title = "C") + ylim(0,1.1*max(best.df.means[,2]))

# gompertz and buchanan removed ranks
best.df.removed <- data.frame(model = colnames(ranks.removed)[-1], best = colSums(ranks.removed==1)[-1])
g4 <- ggplot(data = best.df.removed , aes(x = model, y= best)) + geom_bar(stat = "identity")+ theme_bw()+ scale_x_discrete(limits = colnames(ranks.removed[,-1])) + geom_text(aes(label=round(best,2)), position=position_dodge(width=0.9), vjust=-0.25) + labs(x = "Model", y = "Number of best fit", title="D") + ylim(0,1.1*max(best.df.removed[,2]))

g <- arrangeGrob( g1,g2,g3,g4, nrow = 2)
ggsave("../results/best_fit_frequency_bar.pdf",plot = g, width = 10, height = 8, device = "pdf")



# calculate proportion of points in lag phase
tlag_ratios.filtered <- c()
for (id in IDs.use){
  tlag_ratios.filtered <- c(tlag_ratios.filtered, (sum(data[data$ID == id,]$Time < weighted_tlag[[id]]))  / length(data[data$ID == id,]$Time))
}

# ranking labels
labels <- c(toOrdinal(1),toOrdinal(2),toOrdinal(3),toOrdinal(4),toOrdinal(5),toOrdinal(6))

# Logistic fit ranking vs proportion of points in tlag phase
tlag_vs_logistic_fit_df <- data.frame(logistic_rank = ranks[,2], tlag = tlag_ratios.filtered)
mean_df <- data.frame(logistic_ranking = 1:6, mean <- tapply(tlag_ratios.filtered, ranks[,2], mean, na.rm = T))
g1 <- ggplot(data = tlag_vs_logistic_fit_df, aes(x=tlag, y = logistic_rank)) + geom_jitter(height = 0.1, width = 0, size = 0.3) + theme_bw() + labs(x = expression('Proportion of points in t'[lag]*' phase'), y = "Logistic AIC ranking", title = "A") + geom_point(data = mean_df, aes(x = mean, y = logistic_ranking), color = "red", size = 3) + scale_y_continuous(labels = labels, breaks = 1:6)

#  delta_AIC_logistic_vs_Gompertz against proportion of points in lag phase
AICdif.tlag <- data.frame(dif = AIC.analysis[,1]-AIC.analysis[,2], tlag = tlag_ratios.filtered)
g2 <- ggplot(data = AICdif.tlag, aes(x=tlag, y = dif)) + geom_point(size = .5) +geom_smooth(method = "lm") + theme_bw() + labs(x = expression('Proportion of points in t'[lag]*' phase'), y = expression(Delta~AIC[~Logistic~-~Gompertz]), title = "B")

g.total <- grid.arrange(g1,g2,nrow=1)
ggsave("../results/tlag_vs_logistic.pdf", plot=g.total, width = 8, height = 4)


#
#
# check if there is significant relationship
M.delta_AIC.tlag <- lm ( (AIC.analysis[,1]-AIC.analysis[,2])~ tlag_ratios.filtered)

# removed for project compilation
#summary(M.delta_AIC.tlag)
#plot(M)


#  delta_AIC_logistic_vs_gompertz against temp 

# get vector of temperatures in dataset
temp.vec <- c()
for (id in IDs.use){
  temp.vec <- c(temp.vec, as.numeric(as.character(data[data$ID == id,]$Temp)[1]))
}

AICdif.temp <- data.frame(dif = AIC.analysis[,1]-AIC.analysis[,2], temp = temp.vec)
pdf("../results/logistic_vs_gompertz_temp.pdf", width = 5, height = 5)
ggplot(data = AICdif.temp, aes(x=temp, y = dif)) + geom_jitter(width = .3, height = 0, size= 0.5) +geom_smooth(method = "lm") + theme_bw() + labs(x = "Temperature / °C", y = expression(Delta~AIC[~Logistic~-~Gompertz]))
graphics.off()

# check if there is significant relationship

M.delta_AIC.temp <- lm( (AIC.analysis[,1]-AIC.analysis[,2]) ~ temp.vec)
# removed for project compilation
#summary(M.delta_AIC.temp)
#plot(M.delta_AIC.temp)

# LINEAR SPACE

# get ranks of each mdoel for each id
ranks.linear <- cbind(IDs.use,t(apply(AIC.analysis.linear[,], 1, rank, ties.method='max')))

# consider only gompertz so we can compare 3 param vs 4 param
ranks.removed.linear <- cbind(IDs.use,t(apply(AIC.analysis.linear[,-c(2,4)], 1, rank, ties.method='max')))


# how often is each model best - removed for project compilation
# for (col in colnames(ranks.linear)){
#   print(col)
#   print(sum(ranks.linear[,col] == 1))
#   print(dim(ranks.linear)[1])
# }

# best model bar chart
best.df.linear <- data.frame(model = colnames(ranks[,-1]), best = colSums(ranks.linear[,-1]==1))
g1 <- ggplot(data = best.df.linear , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=best), position=position_dodge(width=0.9), vjust=-0.25)+ theme_bw()+ scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Number of best fit", title = "A") + ylim(0,1.1*max(best.df.linear[,2]))

# second best bar
second.df.linear <- data.frame(model = colnames(ranks[,-1]), best = colSums(ranks.linear[,-1]==2))
g2 <- ggplot(data = second.df.linear , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=best), position=position_dodge(width=0.9), vjust=-0.25)+ theme_bw()+ scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Number of second best fit", title = "B") + ylim(0,1.1*max(second.df.linear[,2]))

# mean rankings bar
best.df.means.linear <- data.frame(model = colnames(ranks[,-1]), best = colMeans(ranks.linear[,-1]))
g3 <- ggplot(data = best.df.means.linear , aes(x = model, y= best)) + geom_bar(stat = "identity", position = 'dodge') + geom_text(aes(label=round(best,2)), position=position_dodge(width=0.9), vjust=-0.25) + theme_bw() + scale_x_discrete(limits = colnames(ranks[,-1])) + labs(x = "Model", y = "Mean ranking", title = "C") + ylim(0,1.1*max(best.df.means.linear[,2]))

# gompertz and buchanan removed bar
best.df.removed.linear <- data.frame(model = colnames(ranks.removed)[-1], best = colSums(ranks.removed.linear==1)[-1])
g4 <- ggplot(data = best.df.removed.linear , aes(x = model, y= best)) + geom_bar(stat = "identity")+ theme_bw()+ scale_x_discrete(limits = colnames(ranks.removed[,-1])) + geom_text(aes(label=round(best,2)), position=position_dodge(width=0.9), vjust=-0.25) + labs(x = "Model", y = "Number of best fit", title="D") + ylim(0,1.1*max(best.df.removed.linear[,2]))

g <- arrangeGrob( g1,g2,g3,g4, nrow = 2)
ggsave("../results/best_fit_frequency_bar_linear.pdf",plot = g, width = 10, height = 8, device = "pdf")


# ranking vs proportion of points in tlag phase
tlag_vs_logistic_fit_df.linear <- data.frame(logistic_rank = ranks.linear[,2], tlag = tlag_ratios.filtered)
mean_df.linear <- data.frame(logistic_ranking = 1:6, mean <- tapply(tlag_ratios.filtered, ranks.linear[,2], mean, na.rm = T))
g1 <- ggplot(data = tlag_vs_logistic_fit_df.linear, aes(x=tlag, y = logistic_rank)) + geom_jitter(height = 0.1, width = 0, size = 0.3) + theme_bw() + labs(x = expression('Proportion of points in t'[lag]*' phase'), y = "Logistic AIC ranking", title = "A") + geom_point(data = mean_df.linear, aes(x = mean, y = logistic_ranking), color = "red", size = 3) + scale_y_continuous(labels = labels, breaks = 1:6)



#  delta_AIC_logistic_vs_Gompertz against proportion of points in lag phase
AICdif <- data.frame(dif = AIC.analysis.linear[,1]-AIC.analysis.linear[,2], tlag = tlag_ratios.filtered)
g2 <- ggplot(data = AICdif, aes(x=tlag, y = dif)) + geom_point(size = .5) +geom_smooth(method = "lm") + theme_bw() + labs(x = expression('Proportion of points in t'[lag]*' phase'), y = expression(Delta~AIC[~Logistic~-~Gompertz]), title = "B")
g.total <- grid.arrange(g1,g2,nrow=1)
ggsave("../results/tlag_vs_logistic_linear.pdf", plot=g.total, width = 8, height = 4)

# Check if there is a significant relationship

M.delta_AIC.tlag.linear <- lm ( (AIC.analysis.linear[,1]-AIC.analysis.linear[,2])~ tlag_ratios.filtered)
# removed for project compilation
#summary(M.delta_AIC.tlag.linear)

# delta_AIC_logistic_vs_Gompertz against temperature
AICdif <- data.frame(dif = AIC.analysis.linear[,1]-AIC.analysis.linear[,2], temp = temp.vec)
pdf("../results/logistic_vs_gompertz_temp_linear.pdf", width = 5, height = 5)
ggplot(data = AICdif, aes(x=temp, y = dif)) + geom_jitter(width = .3, height = 0, size= 0.5) +geom_smooth(method = "lm") + theme_bw() + labs(x = "Temperature / °C", y = expression(Delta~AIC[~Logistic~-~Gompertz]))
graphics.off()

M.delta_AIC.temp.linear <- lm( (AIC.analysis.linear[,1]-AIC.analysis.linear[,2]) ~ temp.vec)
# removed for project compilation
#summary(M.delta_AIC.temp.linear)
#plot(M.delta_AIC.temp.linear)

###############################
# PARAMETER ESTIMATE ANALYSIS #
###############################


# SPLIT IDs BASED ON UNIT

units_df <- data.frame(ID = (1:281),unit = rep(NA, 281))
for (id in IDs){
  units_df[units_df$ID == id, "unit"] <- as.character(data[data$ID == id,"PopBio_units"])[1]
}

# Make parameter dataframe
params.df <- data.frame("id"=NA, "rmax"=NA,"tlag"=NA,"Nmax"=NA,"rmax.linear"=NA,"tlag.linear"=NA,"Nmax.linear"=NA,  "Temp"=NA, "Medium"=NA, "Species"=NA, "Citation"=NA, "Unit" = NA)
l = 1
for (id in IDs.filtered.growth){
  
  rmax <- as.numeric(weighted_rmax[[id]])
  tlag <- as.numeric(weighted_tlag[[id]])
  Nmax <- as.numeric(weighted_Nmax[[id]])
  
  rmax.linear <- as.numeric(weighted_rmax.linear[[id]])
  tlag.linear <- as.numeric(weighted_tlag.linear[[id]])
  Nmax.linear <- as.numeric(weighted_Nmax.linear[[id]])
  
  temp <- as.numeric(as.character((data[data$ID == id,]$Temp)[1]))
  medium <- as.character((data[data$ID == id,]$Medium)[1])
  species <- as.character((data[data$ID == id,]$Species)[1])
  citation <- as.character((data[data$ID == id,]$Citation)[1])
  unit <- units_df[id,"unit"]
  params.df[l,] <-  c(id, rmax,tlag,Nmax,rmax.linear,tlag.linear,Nmax.linear, temp, medium, species, citation, unit)
  l=l+1
}

params.df$Temp <- as.factor(params.df$Temp)
params.df$Unit <- as.factor(params.df$Unit)
params.df$Medium <- as.factor(params.df$Medium)
params.df$Species <- as.factor(params.df$Species)

params.df$rmax <- as.numeric(params.df$rmax)
params.df$tlag <- as.numeric(params.df$tlag)
params.df$rmax.linear <- as.numeric(params.df$rmax.linear)
params.df$tlag.linear <- as.numeric(params.df$tlag.linear)

params.df <- params.df[params.df$rmax > 0,]



########
# LMMs #
########

# all observed temperatures
temps <- c(0,2,4,5,6,7,8,10,12,15,16,20,25,30,32,35,37)
temps.strings <- paste0("Temp", temps)

# LOG SPACE

# r_max

# run model
rmax.mixed <- lmer( log(rmax) ~ Temp + (1|Species) + (1|Medium), data = params.df)
result <- summary(rmax.mixed)$coef

# calculate estimates and s.e for each temperature - requires adding on intercept
estimates.rmax <- c(fixef(rmax.mixed)[1], fixef(rmax.mixed)[temps.strings[-1]])
estimates.rmax[2:length(estimates.rmax)] <- estimates.rmax[2:length(estimates.rmax)] + estimates.rmax[1]
s.e <- c(summary(rmax.mixed)$coef[1,2],summary(rmax.mixed)$coef[,2][temps.strings[-1]])

# plot points
g1 <- ggplot(data = data.frame(Temperature = temps, rmax = estimates.rmax), aes(x= Temperature, y = rmax))+geom_point() + theme_bw() + labs(y = expression(log(mu[max])), x = "Temperature / °C", title = "A")

# plot confidence intervals
for (i in 1:17){
  df <- data.frame(x =rep(temps[i],2), y = c(estimates.rmax[i]+1.96*s.e[i], estimates.rmax[i]-1.96**s.e[i]))
  g1 <- g1 + geom_line(data=df, aes(x=x, y=y), size = .2)
}

# tlag

# run model
tlag.mixed <- lmer( tlag ~ Temp + (1|Species) + (1|Medium), data = params.df)
result <- summary(tlag.mixed)$coef

# calculate estimates and s.e for each temperature - requires adding on intercept
estimates.tlag <- c(fixef(tlag.mixed)[1], fixef(tlag.mixed)[temps.strings[-1]])
estimates.tlag[2:length(estimates.tlag)] <- estimates.tlag[2:length(estimates.tlag)] + estimates.tlag[1]
s.e <- c(summary(tlag.mixed)$coef[1,2],summary(tlag.mixed)$coef[,2][temps.strings[-1]])

# plot points
g2 <- ggplot(data = data.frame(Temperature = temps, tlag = estimates.tlag), aes(x= Temperature, y = tlag))+geom_point() + theme_bw() + labs(y = expression(t[lag]), x = "Temperature / °C", title = "B")

# plot confidence intervals
for (i in 1:17){
  df <- data.frame(x =rep(temps[i],2), y = c(estimates.tlag[i]+1.96*s.e[i], estimates.tlag[i]-1.96*s.e[i]))
  g2 <- g2 + geom_line(data=df, aes(x=x, y=y), size = .2)
}

# make final plot
g <- grid.arrange(g1,g2, nrow = 1)
ggsave("../results/tlags.pdf", g, width = 10, height = 5)


# LINEAR SPACE

# rmax

# run model
rmax.mixed.linear <- lmer( log(rmax.linear) ~ Temp + (1|Species) + (1|Medium), data = params.df)
result <- summary(rmax.mixed.linear)$coef

# calculate estimates and s.e for each temperature - requires adding on intercept
estimates.rmax.linear <- c(fixef(rmax.mixed.linear)[1], fixef(rmax.mixed.linear)[temps.strings[-1]])
estimates.rmax.linear[2:length(estimates.rmax.linear)] <- estimates.rmax.linear[2:length(estimates.rmax.linear)] + estimates.rmax.linear[1]
s.e <- c(summary(rmax.mixed)$coef[1,2],summary(rmax.mixed)$coef[,2][temps.strings[-1]])

# plot points
g1 <- ggplot(data = data.frame(Temperature = temps, rmax = estimates.rmax.linear), aes(x= Temperature, y = rmax))+geom_point() + theme_bw() + labs(y = expression(log(mu[max])), x = "Temperature / °C", title = "A")

# plot confidence intervals
for (i in 1:17){
  df <- data.frame(x =rep(temps[i],2), y = c(estimates.rmax.linear[i]+1.96*s.e[i], estimates.rmax.linear[i]-1.96**s.e[i]))
  g1 <- g1 + geom_line(data=df, aes(x=x, y=y), size = .2)
}

# tlag

# run model
 # medium random effext approx = 0, so singular fit:
tlag.mixed.linear <- lmer( tlag.linear ~ Temp + (1|Species) + (1|Medium), data = params.df, control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
result <- summary(tlag.mixed)$coef

# calculate estimates and s.e for each temperature - requires adding on intercept
estimates.tlag.linear <- c(fixef(tlag.mixed.linear)[1], fixef(tlag.mixed.linear)[temps.strings[-1]])
estimates.tlag.linear[2:length(estimates.tlag.linear)] <- estimates.tlag.linear[2:length(estimates.tlag.linear)] + estimates.tlag.linear[1]
s.e <- c(summary(tlag.mixed)$coef[1,2],summary(tlag.mixed)$coef[,2][temps.strings[-1]])

# plot points
g2 <- ggplot(data = data.frame(Temperature = temps, tlag = estimates.tlag.linear), aes(x= Temperature, y = tlag))+geom_point() + theme_bw() + labs(y = expression(t[lag]), x = "Temperature / °C", title = "B")

# plot confidence intervals
for (i in 1:17){

  df <- data.frame(x =rep(temps[i],2), y = c(estimates.tlag.linear[i]+1.96*s.e[i], estimates.tlag.linear[i]-1.96*s.e[i]))
  g2 <- g2 + geom_line(data=df, aes(x=x, y=y), size = .2)
}

# make final plot
g <- grid.arrange(g1,g2, nrow = 1)
ggsave("../results/tlags_linear.pdf", g, width = 10, height = 5)


####################
# PARAMETER BIASES #
####################

# check that differences are normally distributed
par(mfrow = c(2,2))
for (i in 1:3){
  for (j in 1:3){
    if (i>j){
      n1 <- names(tlag_df)[i]
      n2 <- names(tlag_df)[j]
      hist(tlag_df[IDs.use,i] - tlag_df[IDs.use,j], xlab = paste(n1, "-", n2), main = paste(n1, "-", n2) , breaks = 100)
      
    }

  }
}

G.vs.Ba.t <- t.test(tlag_df[IDs.use,1], tlag_df[IDs.use,2], paired = T)
G.vs.Bu.t <- t.test(tlag_df[IDs.use,1], tlag_df[IDs.use,3], paired = T)
Ba.vs.Bu.t <- t.test(tlag_df[IDs.use,2], tlag_df[IDs.use,3], paired = T)

G.vs.Ba.r <- t.test(rmax_df[IDs.use,1], rmax_df[IDs.use,2], paired = T)
G.vs.Bu.r <- t.test(rmax_df[IDs.use,1], rmax_df[IDs.use,3], paired = T)
Ba.vs.Bu.r <-t.test(rmax_df[IDs.use,2], rmax_df[IDs.use,3], paired = T)


# not printed for project compilation

# G.vs.Ba.t$estimate
# G.vs.Bu.t$estimate
# Ba.vs.Bu.t$estimate
# t.mean <- mean(unlist(weighted_tlag), na.rm = T)
# G.vs.Ba.t$estimate / t.mean
# G.vs.Bu.t$estimate / t.mean
# Ba.vs.Bu.t$estimate / t.mean
# 
# 
# 
# G.vs.Ba.r$estimate
# G.vs.Bu.r$estimate
# Ba.vs.Bu.r$estimate
# r.mean <- mean(unlist(weighted_rmax), na.rm = T)
# G.vs.Ba.r$estimate / r.mean
# G.vs.Bu.r$estimate / r.mean
# Ba.vs.Bu.r$estimate / r.mean

# fit in linear space vs logistic space

Mod <- t.test(unlist(weighted_tlag),unlist(weighted_tlag.linear), paired=T)
