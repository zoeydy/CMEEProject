rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")


###################
# require package #
###################
require(ggplot2)
library(gridExtra)
library(grid)


#################
# read the data #
#################
# Data <- read.csv('../data/mini_pop.csv')
# Data <- Data[order(Data[,'id'], Data[,'Time']),]
# IDs <- unique(Data$id)
# plot data
plot.df <- read.csv("../data/gomp.plot.csv")
# info
infos <- read.csv("../data/gomp.info.csv")
# delete the non-positive parameters
info0 <- subset(infos, infos$rmax > 0 & infos$tlag > 0)
# set digits for temp.c
info0[, "temp_c"] <- round(info0[, "temp_c"], digits = 2)
# change colnames
info0$temp.c <- info0$temp_c
info0$temp.k <- info0$temp_k
# set boltzmann constant 8.617*10^(-5) eV K^(-1)
K = 8.617*10^(-5)
# add 1/KT
info0$mi.one.over.KT <- (-1)/(K*(info0$temp_k))
info0$log.r <- log(info0$rmax)
info0$log.one.tlag <- log(1/info0$tlag)


# ###############################################################
# standardize species name and get validate species data frame #
# ###############################################################
# (more than 3 data points in each data sets subset by species name and has more than one temperature) #

# fit Arrhenius model, get lnA and E calculated from rmax and tlag with SE (within species), using data has T>=30°C
spes <- unique(info0$species)

# same length as spes
stand.spe <- c("Pantoea.agglomerans", "Clavibacter.michiganensis","Stenotrophomonas.maltophilia","Klebsiella.pneumonia","Dickeya.zeae",
               "Pectobacterium.carotovorum","Pantoea.agglomerans","Acinetobacter.clacoaceticus","Stenotrophomonas.maltophilia","Klebsiella.pneumonia",
               "Clavibacter.michiganensis","Chryseobacterium.balustinum","Enterobacter","Bacillus.pumilus","Pseudomonas.fluorescens",
               "Pseudomonas.fluorescens","Acinetobacter.clacoaceticus","Stenotrophomonas.maltophilia","Dickeya.zeae","Bacillus.pumilus",
               "Pantoea.agglomerans","Acinetobacter.clacoaceticus","Tetraselmis.tetrahele","Soil.Microbial.Community","Staphylococcus",
               "Pseudomonas","Aerobic.Psychotropic","Aerobic.Mesophilic" ,"Spoilage","Escherichia.coli",
               "Salmonella.Typhimurium","Curtobacterium.psychrophilum","Cytophaga.antarctica","Cytophaga.xantha","Spirillum.pleomorphum",
               "Micrococcus.cryophilus","Pseudomonas.flourescens","pseudomonas","yeasts.and.moulds","enterobacteriaceae",
               "pseudomonas.sp","enterobacteriaceae","pseudomonas","yeasts.moulds","pseudomonas",
               "yeasts.moulds","enterobacteriaceae","lactic.acid.bacteria","lactic.acid.bacteria","enterobacteriaceae",
               "yeasts.and.moulds","pseudonomads","brochothrix.thermosphacta","aerobic.bacteria","Serratia.marcescens",
               "Arthrobacter","Arthrobacter","Arthrobacter","Arthrobacter.aurescens","Arthrobacter.aurescens",
               "Arthrobacter.globiformis","Arthrobacter.simplex","Lactobacillus.plantarum","Weissella.viridescens","Lactobacillus.sakei",
               "Oscillatoria.agardhii","Oscillatoria.agardhii","Pseudomonas","Pseudomonas.flourescens","Lactobaciulus.plantarum",
               "Flavobacterium","Pseudomonas","Rahnella","Raoultella","Sphingobacterium",
               "Rhizobium","Arthrobacter","Citrobacter","Pseudomonas",
               "Curtobacterium","Microbacterium","IsoMix","Rhizobium","Roseomonas","Novosphingobium",
               "Novosphingobium","Raoultella","Serratia","Chryseobacterium","Chryseobacterium",
               "Paenibacillus","Paenibacillus","Bacillus")
# standardize the species name
info.sta <- data.frame()
for (i in 1:length(spes)) {
  sta1 <- info0[info0$species == spes[i],]
  sta1$sta.spe <- stand.spe[i]
  info.sta <- rbind(info.sta, sta1)
}



####################
# define functions #
####################
# function to calculate the confidence interval
ci_calc <- function(x){
  MEAN <- mean(x)
  SD <- sd(x)
  n <- length(x)
  return(qnorm(0.975)*SD/sqrt(n))
}
# function to save plot
save_plot <- function(filepath, plo){
  pdf(paste0(filepath,".pdf"))
  print(plo)
  graphics.off()
}
save_png <- function(filepath, plo){
  png(paste0(filepath,"png"))
  print(plo)
  graphics.off()
}

# function to group data by temperature in different groups
group_data <- function(dat, n){
  rang <- c(seq(min(dat$temp.c), max(dat$temp.c) + 1, length.out =n+1), max(dat$temp.c))
  temp.info <- data.frame()
  for (i in 1:n) {
    df <- dat[dat$temp.c >= rang[i] & dat$temp.c < rang[i+1], ]
    df$temp.c.group <- paste0(rang[i], "~", rang[i], " °C") 
    df$mean.temp.c <- mean(df$temp.c)
    df$CI.temp.c <- ci_calc(df$temp.c)
    temp.info <- rbind(temp.info, df)
  }
  return(temp.info)
}
# function to calculate the mean log(rmax) and lot(1/tlag) with CI under one specific temperature
mean_rate <- function(dat){
  temps <- unique(dat$temp.c)
  df <- data.frame()
  for (i in 1:length(temps)) {
    tempdf <- dat[dat$temp.c == temps[i], ]
    meandf <- data.frame(temp.c = tempdf$temp.c[1], mean.log.r = mean(tempdf$log.r), CI.r = ci_calc(tempdf$log.r), 
                         mean.log.t = mean(tempdf$log.one.t), CI.t = ci_calc(tempdf$log.one.t))
    df <- rbind(df, meandf)
  }
  return(df)
}
# function to calculate the mean rate with CI using data grouped by temperature
mean_rate_group <- function(dat) {
  return.df <- data.frame()
  temp.group <- unique(dat$temp.c.group)
  for (i in 1:length(temp.group)) {
    df <- dat[dat$temp.c.group == temp.group[i],]
    meandf <- data.frame(temp.c.group = df$temp.c.group[1], 
                         mean.temp.c = df$mean.temp.c[1], CI.temp.c = df$CI.temp.c[1],
                         mean.log.r = mean(df$log.r), CI.r = ci_calc(df$log.r), 
                         mean.log.t = mean(df$log.one.t), CI.t = ci_calc(df$log.one.t),
                         mean.mi.KT = mean(df$mi.one.over.KT),
                         CI.mi.KT = ci_calc(df$mi.one.over.KT)
                         )
    return.df <- rbind(return.df, meandf)
  }
  return(return.df)
}
# function to plot grouped data frame with CI
group_plot <- function(dat){
  p1 <- ggplot(data = dat, aes(x = mean.mi.KT, y = mean.log.r)) +
    geom_point() +
    stat_smooth(data = dat[dat$mean.mi.KT <= -1/(K*(273+30)), ],formula = y~x, method = lm, se = TRUE)+
    geom_errorbar(aes(ymin = mean.log.r - CI.t, ymax = mean.log.r + CI.t, width = .1)) +
    geom_errorbar(aes(xmin = mean.mi.KT - CI.mi.KT, xmax = mean.mi.KT + CI.mi.KT, width = .1)) +
    labs(y="log(tlag)(1/h) +/- CI ", x = "Mean of -1/KT(eV) +/- CI ") +
    geom_ribbon(data = dat[dat$mean.mi.KT >= max(dat[dat$mean.mi.KT <= -1/(K*(273+30)), ]$mean.mi.KT),], aes(ymin = -Inf, ymax = Inf, alpha = .3)) +
    theme(legend.position = 'none') 
    #scale_x_continuous(name = "Temperature Range", labels = dat$GROUP, breaks=dat$mean.mi.KT) +
    # labs(tag = "A")
  p2 <- ggplot(data = dat, aes(x = mean.mi.KT, y = mean.log.t)) +
    geom_point() +
    stat_smooth(data = dat[dat$mean.mi.KT <= -1/(K*(273+30)), ],formula = y~x, method = lm, se = TRUE)+
    geom_errorbar(aes(ymin = mean.log.t - CI.t, ymax = mean.log.t + CI.t, width = .1)) +
    geom_errorbar(aes(xmin = mean.mi.KT - CI.mi.KT, xmax = mean.mi.KT + CI.mi.KT, width = .1)) +
    labs(y="log(tlag)(1/h) +/- CI ", x = "Mean of -1/KT(eV) +/- CI ") +
    geom_ribbon(data = dat[dat$mean.mi.KT >= max(dat[dat$mean.mi.KT <= -1/(K*(273+30)), ]$mean.mi.KT),], aes(ymin = -Inf, ymax = Inf, alpha = .3)) +
    theme(legend.position = 'none') 
    #scale_x_continuous(name = "Temperature Range", labels = dat$GROUP, breaks=dat$mean.mi.KT) +
    # labs(tag = "B")
  grid.arrange(p1,p2,nrow = 1)
}


#################################
# group plot rate VS temp_group #
#################################
info4 <- group_data(info.sta, 4)
info5 <- group_data(info.sta, 5)
info6 <- group_data(info.sta, 6)
info7 <- group_data(info.sta, 7)
info8 <- group_data(info.sta, 8)
mean.info4 <- mean_rate_group(info4)
mean.info5 <- mean_rate_group(info5)
mean.info6 <- mean_rate_group(info6)
mean.info7 <- mean_rate_group(info7)
mean.info8 <- mean_rate_group(info8)
g1 <- group_plot(mean.info4)
g2 <- group_plot(mean.info5)
g3 <- group_plot(mean.info6)
g4 <- group_plot(mean.info7)
g5 <- group_plot(mean.info8)
# 1. plot
pdf("../results/arrhenius/mean_rate_temp_group.pdf")
p <- grid.arrange(
  g2,g3,g4,g5, nrow = 4,
  top = "Across species mean fitness value VS Temperature(Grouped)",
  bottom = textGrob(
    "CI is calculated from normal distribution",
    gp = gpar(fontface = 6, fontsize = 9),
    hjust = 1,
    x = 1
  )
)
print(p)
graphics.off()
# # 2. calculation: fit arrhenius on mean value in each temperature group (use mean.info5)
# double_mean <- function(dat){
#   dat$mi.one.over.KT <- -1/(K*(dat$mean.temp.c+273))
#   fit_double_mean_r <- lm(dat[mean.info5$dat <= 20,], formula = mean.log.r ~ mi.one.over.KT)
#   sum_double_mean_r <- summary(fit_double_mean_r)
#   return(data.frame(summary(fit_double_mean_r)$coefficients)[2,c(1,2,4)])
# }
# e_double_mean5 <- double_mean(mean.info5)

##########################################################################################
# log of rmax and 1/tlag in mean value with CI VS temperature(°C) less than 20°C  #
########################################################################################## 
mean_df <- data.frame()
temp.c.s <- unique(info.sta$temp.c)
for (i in 1:length(temp.c.s)) {
  df1 <- subset(info.sta, info.sta$temp.c == temp.c.s[i])
  mean1 <- data.frame(temp.c = df1$temp.c[1],
                      mean.rmax = mean(log(df1$rmax)),
                      mean.1tlag = mean(log(1/df1$tlag)),
                      ci.rmax = ci_calc(log(df1$rmax)),
                      ci.1tlag = ci_calc(log(1/df1$tlag)))
  mean_df <- rbind(mean_df, mean1)
}
# 1.plot
# rmax
r.temp <- mean_df[mean_df$mean.rmax == max(mean_df$mean.rmax),]$temp.c
r.temp <- 30
p <- ggplot(data = mean_df, aes(x = temp.c, y =mean.rmax )) +
  geom_point() +
  stat_smooth(data = mean_df[mean_df$temp.c <= r.temp,], 
              formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.rmax - ci.rmax, 
                    ymax=mean.rmax + ci.rmax,
                    width = 0.1)) +
  geom_ribbon(data = mean_df[mean_df$temp.c >= r.temp,], aes(ymin=-Inf, ymax=Inf, alpha = .5 )) +
  theme(legend.position = 'none') +
  labs(title = "Mean of log(rmax) VS temperature(°C) less than T_opt", x = "Temperature(°C)", y = "Mean of log(rmax) +/- CI")
save_plot("../results/arrhenius/mean_log_r_temp",p)

# tlag
# r.temp <- mean_df[mean_df$mean.rmax == max(mean_df$mean.rmax),]$temp.c
t.temp <- 30
p <- ggplot(data = mean_df, aes(x = temp.c, y =mean.1tlag )) +
  geom_point() +
  stat_smooth(data = mean_df[mean_df$temp.c <= t.temp,],formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.1tlag - ci.1tlag, 
                    ymax=mean.1tlag + ci.1tlag,
                    width = 0.1)) +
  geom_ribbon(data = mean_df[mean_df$temp.c >= t.temp, ], aes(ymin = -Inf, ymax = Inf, alpha = 0.1)) +
  theme(legend.position = 'none') +
  labs(title = "Mean of log(1/tlag) VS temperature(°C) less than T_opt", x = "Temperature(°C)",y = "Mean of log(1/tlag) +/- CI")
save_plot("../results/arrhenius/mean_log_tlag_temp",p)


########################################################
# grouping by species get E with CI #
########################################################
arrhe.df <- data.frame()
# vali.spe <- data.frame()
vali.info <- data.frame()
spes <- unique(info.sta$sta.spe)
for (i in 1:length(spes)) {
  # subset
  spe.df <- info.sta[info.sta$sta.spe == spes[i],]
  
  if (nrow(spe.df) < 3 | length(unique(spe.df$temp.c)) < 2) {
    arrhe.df <- arrhe.df
  }else{
    # vali.spe <- rbind(vali.spe, spe.df)
    vali.info <- rbind(vali.info, spe.df)
    
    fit.r <- lm(data = spe.df, formula = log.r~mi.one.over.KT)
    fit.r.df <- data.frame(summary(fit.r)$coefficients)
    fit.r.df$id <- rep(spe.df$id[1],2)
    fit.r.df$from = rep("rmax",2)
    fit.r.df$param = c("lnA","E")
    fit.r.df$sta.spe = spe.df$sta.spe[1]
    
    fit.t <- lm(data = spe.df, formula = log.one.tlag~mi.one.over.KT)
    fit.t.df <- data.frame(summary(fit.t)$coefficients)
    fit.t.df$id = rep(spe.df$id[1],2)
    fit.t.df$from = rep("tlag",2)
    fit.t.df$param = c("lnA","E")
    fit.t.df$sta.spe = spe.df$sta.spe[1]
    
    arrhe.df <- rbind(arrhe.df, fit.r.df, fit.t.df)
  }
}

###### get fit_validated_data_fram
fit.vali <- subset(arrhe.df, arrhe.df$Pr...t.. < 0.05 & arrhe.df$Estimate > 0)
spes <- unique(fit.vali$sta.spe)
fit.vali.info <- data.frame()
for (i in 1:length(spes)) {
  infodf <- info.sta[info.sta$sta.spe == spes[i],]
  fit.vali.info <- rbind(fit.vali.info, infodf)
  # #########################################
  # # plot arrhenius model (within species) #
  # #########################################
  # if (i == 25|26) {
  #   tpk.r <- spe.df[spe.df$log.r == max(spe.df$log.r),]$temp.c
  #   p <- ggplot(data = spe.df, aes(x = mi.one.over.KT, y = log.r)) +
  #     geom_point() +
  #     stat_smooth(data = spe.df[spe.df$temp.c <= tpk.r,], method = 'lm', formula = y~x) +
  #     geom_ribbon(data = spe.df[spe.df$temp.c >= tpk.r,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
  #   save.png(paste0("../results/arrhenius/arrhe_fit/r_spe",i), p)
  # 
  #   tpk.t <- spe.df[spe.df$log.one.t == max(spe.df$log.one.t),]$temp.c
  #   p <- ggplot(data = spe.df, aes(x = mi.one.over.KT, y = log.one.t)) +
  #     geom_point() +
  #     stat_smooth(data = spe.df[spe.df$temp.c <= tpk.t,], method = 'lm', formula = y~x) +
  #     geom_ribbon(data = spe.df[spe.df$temp.c >= tpk.t,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
  #   save.png(paste0("../results/arrhenius/arrhe_fit/t_spe",i), p)
  # } else{
    p <- ggplot(data = infodf, aes(x = mi.one.over.KT, y = log.r)) +
      geom_point() +
      stat_smooth(data = infodf[infodf$temp.c <= 30,], method = 'lm', formula = y~x)
  if (nrow(infodf[infodf$temp.c >= 30,] >= 1)) {
    p <- p + geom_ribbon(data = infodf[infodf$temp.c >= 30,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
    save_png(paste0("../results/arrhenius/arrhe_fit/r_spe",i), p)
  }else{
    save_png(paste0("../results/arrhenius/arrhe_fit/r_spe",i), p)
  }
      

    p <- ggplot(data = infodf, aes(x = mi.one.over.KT, y = log.one.tlag)) +
      geom_point() +
      stat_smooth(data = infodf[infodf$temp.c <= 30,], method = 'lm', formula = y~x) +
      labs(title="log(1/tlag) VS 1/KT",x="1/KT", y = "log(1/tlag)")
    if (nrow(infodf[infodf$temp.c >= 30,])>=1) {
      p <- p + geom_ribbon(data = infodf[infodf$temp.c >= 30,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3) 
      save_png(paste0("../results/arrhenius/arrhe_fit/t_spe",i), p)
    }else {
      save_png(paste0("../results/arrhenius/arrhe_fit/t_spe",i), p)
    }
      
      
    
  # }
}

# # get the E_L and corresponding CI using whole data set
# fit.r.spe <- lm(info0, formula = log.r ~ mi.one.over.KT)
# fit.t.spe <- lm(info0, formula = log.one.t~mi.one.over.KT)
# # summary(fit.r)$coefficients
# # summary(fit.t)$coefficients
# # CI = mean +/- 2SE
# CI.r.spe <- confint(fit.r, level = 0.95)
# CI.t.spe <- confint(fit.t, level = 0.95)



##########################################################
# plot E_s and lnA_s with E_L  and lnA_L(acrose species) #
##########################################################
# plot the E and lnA estimated from rmax and tlag, using data with temperature not larger than 30°C 
# and grouped by species validated by checking data points has more than one validate temperature and 3 points
# add dashed line representing across species E and lnA value with 95% level CI

# get the E_L and corresponding CI using validate species data set
fit.r.spe <- lm(fit.vali.info, formula = log.r ~ mi.one.over.KT)
fit.t.spe <- lm(fit.vali.info, formula = log.one.tlag~mi.one.over.KT)
CI.r.spe <- confint(fit.r.spe, level = 0.95)
CI.t.spe <- confint(fit.t.spe, level = 0.95)

edf <- fit.vali[fit.vali$param == "E",]
# save E value tables to latex
E_r <- subset(edf, edf$from == "rmax")
E_t <- subset(edf, edf$from == "tlag")

E_rmax <- data.frame(Species = E_r$sta.spe, "E value" = E_r$Estimate, 
                     "Confidence Interval" = 2*E_r$Std..Error, "P value" = E_r$Pr...t..)
E_rmax_ave <- data.frame(Species = "Mean Value", "E value" = mean(E_r$Estimate), 
                         "Confidence Interval" = ci_calc(E_r$Estimate), "P value" = mean(E_r$Pr...t..))
E_rmax <- rbind(E_rmax, E_rmax_ave)
E_tlag <- data.frame(Species = E_t$sta.spe, "E value" = E_t$Estimate, 
                     "Confidence Interval" = 2*E_t$Std..Error, "P value" = E_t$Pr...t..)
E_tlag_ave <- data.frame(Species = "Mean Value", "E value" = mean(E_t$Estimate), 
                         "Confidence Interval" = ci_calc(E_t$Estimate), "P value" = mean(E_t$Pr...t..))
E_tlag <- rbind(E_tlag, E_tlag_ave)
library(xtable)
print(xtable(E_rmax, type = "latex", digits = -10), file = "../results/arrhenius/E_rmax_table.tex")
print(xtable(E_tlag, type = "latex", digits = -10), file = "../results/arrhenius/E_tlag_table.tex")

mean.e.r <- mean(E_r$Estimate)
mean.e.t <- mean(E_t$Estimate)
mean.e.r.ci <- ci_calc(E_r$Estimate)
mean.e.t.ci <- ci_calc(E_t$Estimate)
p <- ggplot(data = edf, aes(x = Estimate, y = sta.spe, colour = from)) +
  geom_point()+
  geom_errorbar(data = edf, aes(xmin = Estimate - 2*Std..Error, xmax = Estimate + 2*Std..Error)) +
  # geom_rect(aes(xmin=CI.r.spe[2,1], xmax=CI.r.spe[2,2], ymin=-Inf, ymax=Inf), color = 'red') +
  # geom_rect(aes(xmin=CI.t.spe[2,1], xmax=CI.t.spe[2,2], ymin=-Inf, ymax=Inf), color = 'blue') +
  # add CI of E_L
  annotate("rect", xmin=CI.r.spe[2,1], xmax=CI.r.spe[2,2], ymin=-Inf, ymax=Inf, alpha = .2) +
  annotate("rect", xmin=CI.t.spe[2,1], xmax=CI.t.spe[2,2], ymin=-Inf, ymax=Inf, alpha = .2) +
  # add E_L slope from fitting the whole validate data set
  geom_vline(xintercept=fit.r.spe$coefficients[2], linetype="dashed", color = "red") +
  geom_vline(xintercept=fit.t.spe$coefficients[2], linetype="dashed", color = "blue") +
  annotate("text", label = "E_L estimated from rmax",
             size = 2, x = fit.r.spe$coefficients[2]+3, y = "Pseudomonas") +
  geom_segment(aes(x = fit.r.spe$coefficients[2], y = "Pseudomonas", xend = fit.r.spe$coefficients[2]+2, yend = "Pseudomonas"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  annotate("text", label = "E_L estimated from tlag",size = 2, x = fit.t.spe$coefficients[2]+2, y = "pseudomonas.sp") +
  geom_segment(aes(x = fit.t.spe$coefficients[2], y = "pseudomonas.sp", xend = fit.t.spe$coefficients[2]+1, yend = "pseudomonas.sp"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  # add E_s(mean of E_species)
  geom_vline(xintercept=mean.e.r, linetype="dashed", color = "darkred") +
  geom_vline(xintercept=mean.e.t, linetype="dashed", color = "darkblue") +
  annotate("text", label = "E_S estimated from rmax",size = 2, x = mean.e.r+2, y = "Staphylococcus") +
  geom_segment(aes(x = mean.e.r, y = "Staphylococcus", xend = mean.e.r+1, yend = "Staphylococcus"), 
               arrow = arrow(length = unit(0.1, "cm")), colour = 'black') +
  annotate("text", label = "E_S estimated from tlag", size = 2, x = mean.e.t+2, y = "Lactobaciulus.plantarum") +
  geom_segment(aes(x = mean.e.t, y = "Lactobaciulus.plantarum", xend = mean.e.t+1, yend = "Lactobaciulus.plantarum"), 
               arrow = arrow(length = unit(0.1, "cm")), colour = 'black') +
  # add CI of E_s
  annotate("rect", xmin=mean.e.r-mean.e.r.ci, xmax=mean.e.r+mean.e.r.ci,ymin=-Inf, ymax=Inf, alpha = .2) +
  annotate("rect", xmin=mean.e.t-mean.e.t.ci, xmax=mean.e.t+mean.e.t.ci,ymin=-Inf, ymax=Inf, alpha = .2) 
  # scale_x_continuous(name = "Thermal Sensitivity", limits = c(-1.3*10^(-18), 1*10^(-18)))
# save_plot("../results/arrhenius/E_spe", p)
save_plot("../results/arrhenius/E_spe", p)

# plot histogram

p <- ggplot(edf, aes(x=Estimate, fill=from)) +
  geom_histogram() +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #theme_ipsum() +
  # geom_bar(stat="count")
  labs(fill="")

################################################
# add mean and medium data
################################################
save_plot("../results/arrhenius/E_hist", p)

######################################
# fitter plot ES, ES_average, and EL #
######################################
E_long <- data.frame(
  SL = "Long-term",
  from = c('rmax','tlag'),
  E = c(fit.r.spe$coefficients[2],fit.t.spe$coefficients[2]),
  E_min = c(CI.r.spe[2,1], CI.t.spe[2,1]),
  E_max = c(CI.r.spe[2,2], CI.t.spe[2,2]))
E_short <- data.frame(
  SL = "Short-term",
  from = c('rmax','tlag'),
  E = c(E_rmax[11,2],E_tlag[16,2]),
  E_min = c(E_rmax[11,2]-E_rmax[11,3], E_rmax[11,2]-E_rmax[11,3]),
  E_max = c(E_rmax[11,2]+E_rmax[11,3], E_rmax[11,2]+E_rmax[11,3]),
  medi = c(median(E_rmax$E.value[-length(E_rmax$E.value)]), median(E_tlag$E.value[-length(E_tlag$E.value)])))

hist_df <- data.frame(value = c(E_short$E, E_short$medi), 
                      label = c("Mean E estimated from rmax", "Mean E estimated from tlag",
                      "Median E estimated from rmax","Median E estimated from tlag"))

E_short <- E_short[,-ncol(E_short)]
E_plot <- rbind(E_short,E_long)

# histogram of the activation energy 
hist.linetype <- c("twodash", "twodash","dotted", "dotted")
hist.color <- c("#D55E00","royalblue3","#D55E00","royalblue3",'#D55E00',"royalblue3")
p <- ggplot(edf, aes(x=Estimate, color=from, fill = from, alpha = .5)) +
  geom_density() +
  geom_vline(data = hist_df, aes(xintercept=value, color=label)
             , linetype = hist.linetype, colour = c("#D55E00","royalblue3","#D55E00","royalblue3")
             ) +
  scale_color_manual(values=hist.color)
# twodash linetype represents mean value, dotted linetype represents median value
save_plot("../results/arrhenius/E_hist", p)
  
# activation energy acorss species short- and long-term plot
Species_Overlap_Plot  <- ggplot(E_long, aes(x=from)) +
  geom_point(data=edf, aes(y=Estimate, colour='3', shape = '3'), alpha = 0.7, position = position_jitter(width = 0.4, height = 0.0)) +
  geom_errorbar(data = E_short,aes(ymax=E_max, ymin=E_min), width=0.15, size=1, position="dodge") +
  geom_errorbar(data = E_long,aes(ymax=E_max, ymin=E_min), width=0.15, size=1, position="dodge") +
  geom_point(data = E_plot,aes(y=E, colour=SL, shape = SL), size = 4, position="dodge") +  
  xlab('Parameter') +
  ylab('Activation Energy (E)') +
  theme(legend.title=element_blank(),
        legend.justification=c(0,1), legend.position=c(0,1),
        plot.margin = unit(c(1,0,1,1), "cm"))+
  scale_colour_manual(name = "Colour", values =c('#E69F00','#D55E00','royalblue3'),
                      labels = c(expression(italic('E'['S'])),
                                 expression(italic('E'['L'])),
                                 quote(italic('\U0112'[S])))
  ) +
  # scale_shape(guide = FALSE)+
  #geom_hline(aes(yintercept=0.65), linetype='dotted') +
  scale_shape_manual(name = "Shape",values =c(19,15,17),
                     labels = c(expression(italic('E'['S'])),
                                expression(italic('E'['L'])),
                                quote(italic('\U0112'[S])))
  ) +
  coord_cartesian(ylim = c(0,5))
  # scale_shape_manual(guide = 'legend',
  #                    values =c('3'=19,
  #                              'ES'=15,
  #                              'EL'=17
  #                              ),
  #                    labels = c(expression(italic('E'['S'])),
  #                               expression(italic('E'['L'])),
  #                               quote(italic('\U0112'[S])))
  # ) +
  # main_theme
save_plot("../results/arrhenius/E_comperasion", Species_Overlap_Plot)


# # lnA
# adf <- fit.vali[fit.vali$param == "lnA",]
# mean.a.r <- mean(adf[adf$from == "rmax",]$Estimate)
# mean.a.t <- mean(adf[adf$from == "tlag",]$Estimate)
# mean.a.r.ci <- ci_calc(adf[adf$from == "rmax",]$Estimate)
# mean.a.t.ci <- ci_calc(adf[adf$from == "tlag",]$Estimate)
# p <- ggplot(data = adf, aes(y = species, x = Estimate, colour = from)) +
#   geom_point()+
#   geom_errorbar(data = adf, aes(xmin = Estimate - 2*Std..Error, xmax = Estimate + 2*Std..Error)) +
#   # add CI of lnA_L
#   annotate("rect", xmin=CI.r.spe[1,1], xmax=CI.r.spe[1,2], ymin=-Inf, ymax=Inf, alpha = .2) +
#   annotate("rect", xmin=CI.t.spe[1,1], xmax=CI.t.spe[1,2], ymin=-Inf, ymax=Inf, alpha = .2) +
#   # add lnA_L slope from fitting the whole validate data set
#   geom_vline(xintercept=fit.r.spe$coefficients[1], linetype="dashed", color = "red") +
#   geom_vline(xintercept=fit.t.spe$coefficients[1], linetype="dashed", color = "blue") +
#   annotate("text", label = "lnA_L estimated from rmax",size = 2, x = fit.r.spe$coefficients[1]-130, y = "Stenotrophomonas.maltophilia") +
#   geom_segment(aes(x = fit.r.spe$coefficients[1], y = "Stenotrophomonas.maltophilia", xend = fit.r.spe$coefficients[1]-50, yend = "Stenotrophomonas.maltophilia"), 
#                arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
#   annotate("text", label = "lnA_L estimated from tlag", size = 2,x = fit.r.spe$coefficients[1]+150, y = "Tetraselmis.tetrahele") +
#   geom_segment(aes(x = fit.t.spe$coefficients[1], y = "Tetraselmis.tetrahele", xend = fit.r.spe$coefficients[1]+70, yend = "Tetraselmis.tetrahele"), 
#                arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
#   # add lnA_s(mean of E_species)
#   geom_vline(xintercept=mean.a.r, linetype="dashed", color = "darkred") +
#   geom_vline(xintercept=mean.a.t, linetype="dashed", color = "darkblue") +
#   annotate("text", label = "lnA_S estimated from rmax",size = 2, x = mean.a.r-130, y = "Serratia.marcescens") +
#   geom_segment(aes(x = mean.a.r, y = "Serratia.marcescens", xend = mean.a.r-50, yend = "Serratia.marcescens"), 
#                arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
#   annotate("text", label = "lnA_S estimated from tlag",size = 2, x = mean.a.t+150, y = "Oscillatoria.agardhii") +
#   geom_segment(aes(x = mean.a.t, y = "Oscillatoria.agardhii", xend = mean.a.t+50, yend = "Oscillatoria.agardhii"), 
#                arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
#   # add CI of lnA_s
#   annotate("rect", xmin=mean.a.r-mean.a.r.ci, xmax=mean.a.r+mean.a.r.ci,ymin=-Inf, ymax=Inf, alpha = .2) +
#   annotate("rect", xmin=mean.a.t-mean.a.t.ci, xmax=mean.a.t+mean.a.t.ci,ymin=-Inf, ymax=Inf, alpha = .2) 
#   
# # save_plot("../results/arrhenius/lnA_spe", p)
# save_plot("../results/arrhenius/lnA_spe_nostaspe", p)

################################################
# Sharpe-Schoolfield fitting #
################################################
