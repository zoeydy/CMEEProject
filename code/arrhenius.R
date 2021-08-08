rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")


#################
# read the data #
#################
# Data <- read.csv('../data/mini_pop.csv')
# Data <- Data[order(Data[,'id'], Data[,'Time']),]
infos <- read.csv("../data/gomp.info.csv")
# delete the non-positive parameters
info0 <- subset(infos, infos$rmax > 0 & infos$tlag > 0)
# set boltzmann constant 8.617*10^(-5) eV K^(-1)
K = 8.617*10^(-5)
# add 1/KT
info0$one.over.KT <- 1/(K*(info0$temp_k))
info0$log.r <- log(info0$rmax)
info0$log.one.t <- log(1/info0$tlag)


####################
# define functions #
####################
# function to calculate the confidence interval
ci.calc <- function(x){
  MEAN <- mean(x)
  SD <- sd(x)
  n <- length(x)
  return(qnorm(0.975)*SD/sqrt(n))
}
# function to calculate the mean log(rmax) and lot(1/tlag) with CI under one specific temperature
MEAN.rate <- function(dat){
  temps <- unique(dat$temp_c)
  df <- data.frame()
  for (i in 1:length(temps)) {
    tempdf <- dat[dat$temp_c == temps[i], ]
    meandf <- data.frame(temp_c = tempdf$temp_c[1], mean.log.r = mean(tempdf$log.r), CI.r = ci.calc(tempdf$log.r), 
                         mean.log.t = mean(tempdf$log.one.t), CI.t = ci.calc(tempdf$log.one.t))
    df <- rbind(df, meandf)
  }
  return(df)
}
# function to save plot
save.plot <- function(filepath, plo){
  pdf(paste0(filepath,".pdf"))
  print(plo)
  graphics.off()
}
save.png <- function(filepath, plo){
  png(paste0(filepath,"png"))
  print(plo)
  graphics.off()
}

################################################################
# standardize species name and get validate species data frame #
################################################################
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

########################################################
# fit arrhenius model by species get E and lnA with SE #
########################################################
arrhe.df <- data.frame()
vali.spe <- data.frame()
spes <- unique(info.sta$sta.spe)
for (i in 1:length(spes)) {
  # subset
  spe.df <- info.sta[info.sta$species == spes[i],]
  
  if (nrow(spe.df) < 3 | length(unique(spe.df$temp_c)) < 2) {
    arrhe.df <- arrhe.df
  }else{
    vali.spe <- rbind(vali.spe, spe.df)
    
    fit.r <- lm(data = spe.df, formula = log.r~one.over.KT)
    fit.r.df <- data.frame(summary(fit.r)$coefficients)
    fit.r.df$id <- rep(spe.df$id[1],2)
    fit.r.df$from = rep("rmax",2)
    fit.r.df$param = c("lnA","E")
    fit.r.df$sta.spe = spe.df$sta.spe[1]
    
    fit.t <- lm(data = spe.df, formula = log.one.t~one.over.KT)
    fit.t.df <- data.frame(summary(fit.t)$coefficients)
    fit.t.df$id = rep(spe.df$id[1],2)
    fit.t.df$from = rep("tlag",2)
    fit.t.df$param = c("lnA","E")
    fit.t.df$sta.spe = spe.df$sta.spe[1]
    
    arrhe.df <- rbind(arrhe.df, fit.r.df, fit.t.df)
    
    #########################################
    # plot arrhenius model (within species) #
    #########################################
    if (i == 25|26) {
      tpk.r <- spe.df[spe.df$log.r == max(spe.df$log.r),]$temp_c
      p <- ggplot(data = spe.df, aes(x = one.over.KT, y = log.r)) +
        geom_point() +
        stat_smooth(data = spe.df[spe.df$temp_c <= tpk.r,], method = 'lm', formula = y~x) +
        geom_ribbon(data = spe.df[spe.df$temp_c >= tpk.r,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
      save.png(paste0("../results/arrhenius/arrhe_fit/r_spe",i), p)

      tpk.t <- spe.df[spe.df$log.one.t == max(spe.df$log.one.t),]$temp_c
      p <- ggplot(data = spe.df, aes(x = one.over.KT, y = log.one.t)) +
        geom_point() +
        stat_smooth(data = spe.df[spe.df$temp_c <= tpk.t,], method = 'lm', formula = y~x) +
        geom_ribbon(data = spe.df[spe.df$temp_c >= tpk.t,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
      save.png(paste0("../results/arrhenius/arrhe_fit/t_spe",i), p)
    } else{
      p <- ggplot(data = spe.df, aes(x = one.over.KT, y = log.r)) +
        geom_point() +
        stat_smooth(data = spe.df[spe.df$temp_c <= 30,], method = 'lm', formula = y~x) +
        geom_ribbon(data = spe.df[spe.df$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
      save.png(paste0("../results/arrhenius/arrhe_fit/r_spe",i), p)

      p <- ggplot(data = spe.df, aes(x = one.over.KT, y = log.one.t)) +
        geom_point() +
        stat_smooth(data = spe.df[spe.df$temp_c <= 30,], method = 'lm', formula = y~x) +
        geom_ribbon(data = spe.df[spe.df$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), alpha = 0.3)
        # labs(title="log(1/tlag) VS 1/KT (<30°C)",x="1/KT", y = "log(1/tlag)")
      save.png(paste0("../results/arrhenius/arrhe_fit/t_spe",i), p)
    }
  }
}


# # get the E_L and corresponding CI using whole data set
# fit.r.spe <- lm(info0, formula = log.r ~ one.over.KT)
# fit.t.spe <- lm(info0, formula = log.one.t~one.over.KT)
# # summary(fit.r)$coefficients
# # summary(fit.t)$coefficients
# # CI = mean +/- 2SE
# CI.r.spe <- confint(fit.r, level = 0.95)
# CI.t.spe <- confint(fit.t, level = 0.95)

# get the E_L and corresponding CI using validate species data set
fit.r.spe <- lm(vali.spe, formula = log.r ~ one.over.KT)
fit.t.spe <- lm(vali.spe, formula = log.one.t~one.over.KT)
CI.r.spe <- confint(fit.r.spe, level = 0.95)
CI.t.spe <- confint(fit.t.spe, level = 0.95)

##########################################################
# plot E_s and lnA_s with E_L  and lnA_L(acrose species) #
##########################################################
# plot the E and lnA estimated from rmax and tlag, using data with temperature not larger than 30°C 
# and grouped by species validated by checking data points has more than one validate temperature and 3 points
# add dashed line representing across species E and lnA value with 95% level CI

edf <- arrhe.df[arrhe.df$param == "E",]
# save E value tables to latex
E_r <- subset(edf, edf$from == "rmax")
E_t <- subset(edf, edf$from == "tlag")
E_rmax <- data.frame(Species = E_r$sta.spe, "E value" = E_r$Estimate, 
                     "Confidence Interval" = 2*E_r$Std..Error, "P value" = E_r$Pr...t..)
E_tlag <- data.frame(Species = E_t$sta.spe, "E value" = E_t$Estimate, 
                     "Confidence Interval" = 2*E_t$Std..Error, "P value" = E_t$Pr...t..)
library(xtable)
print(xtable(E_rmax, type = "latex", digits = -10), file = "../results/arrhenius/E_rmax_table.tex")
print(xtable(E_tlag, type = "latex", digits = -10), file = "../results/arrhenius/E_tlag_table.tex")

mean.e.r <- mean(edf[edf$from == "rmax",]$Estimate)
mean.e.t <- mean(edf[edf$from == "tlag",]$Estimate)
mean.e.r.ci <- ci.calc(edf[edf$from == "rmax",]$Estimate)
mean.e.t.ci <- ci.calc(edf[edf$from == "tlag",]$Estimate)
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
  annotate("text", label = "E_L estimated from rmax",size = 2, x = 7*10^(-19), y = "Stenotrophomonas.maltophilia") +
  geom_segment(aes(x = fit.r.spe$coefficients[2], y = "Stenotrophomonas.maltophilia", xend = 3.5*10^(-19), yend = "Stenotrophomonas.maltophilia"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  annotate("text", label = "E_L estimated from tlag",size = 2, x = -7*10^(-19), y = "Tetraselmis.tetrahele") +
  geom_segment(aes(x = fit.t.spe$coefficients[2], y = "Tetraselmis.tetrahele", xend = -3.5*10^(-19), yend = "Tetraselmis.tetrahele"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  # add E_s(mean of E_species)
  geom_vline(xintercept=mean.e.r, linetype="dashed", color = "darkred") +
  geom_vline(xintercept=mean.e.t, linetype="dashed", color = "darkblue") +
  annotate("text", label = "E_S estimated from rmax",size = 2, x = 7*10^(-19), y = "Serratia.marcescens") +
  geom_segment(aes(x = mean.e.r, y = "Serratia.marcescens", xend = 3.5*10^(-19), yend = "Serratia.marcescens"), 
               arrow = arrow(length = unit(0.1, "cm")), colour = 'black') +
  annotate("text", label = "E_S estimated from tlag",size = 2, x = -7*10^(-19), y = "Oscillatoria.agardhii") +
  geom_segment(aes(x = mean.e.t, y = "Oscillatoria.agardhii", xend = -3.5*10^(-19), yend = "Oscillatoria.agardhii"), 
               arrow = arrow(length = unit(0.1, "cm")), colour = 'black') +
  # add CI of E_s
  annotate("rect", xmin=mean.e.r-mean.e.r.ci, xmax=mean.e.r+mean.e.r.ci,ymin=-Inf, ymax=Inf, alpha = .2) +
  annotate("rect", xmin=mean.e.t-mean.e.t.ci, xmax=mean.e.t+mean.e.t.ci,ymin=-Inf, ymax=Inf, alpha = .2) 
  # scale_x_continuous(name = "Thermal Sensitivity", limits = c(-1.3*10^(-18), 1*10^(-18)))

save.plot("../results/arrhenius/E_spe", p)


adf <- arrhe.df[arrhe.df$param == "lnA",]
mean.a.r <- mean(adf[adf$from == "rmax",]$Estimate)
mean.a.t <- mean(adf[adf$from == "tlag",]$Estimate)
mean.a.r.ci <- ci.calc(adf[adf$from == "rmax",]$Estimate)
mean.a.t.ci <- ci.calc(adf[adf$from == "tlag",]$Estimate)
p <- ggplot(data = adf, aes(y = sta.spe, x = Estimate, colour = from)) +
  geom_point()+
  geom_errorbar(data = adf, aes(xmin = Estimate - 2*Std..Error, xmax = Estimate + 2*Std..Error)) +
  # add CI of lnA_L
  annotate("rect", xmin=CI.r.spe[1,1], xmax=CI.r.spe[1,2], ymin=-Inf, ymax=Inf, alpha = .2) +
  annotate("rect", xmin=CI.t.spe[1,1], xmax=CI.t.spe[1,2], ymin=-Inf, ymax=Inf, alpha = .2) +
  # add lnA_L slope from fitting the whole validate data set
  geom_vline(xintercept=fit.r.spe$coefficients[1], linetype="dashed", color = "red") +
  geom_vline(xintercept=fit.t.spe$coefficients[1], linetype="dashed", color = "blue") +
  
  annotate("text", label = "lnA_L estimated from rmax",size = 2, x = fit.r.spe$coefficients[1]-130, y = "Stenotrophomonas.maltophilia") +
  geom_segment(aes(x = fit.r.spe$coefficients[1], y = "Stenotrophomonas.maltophilia", xend = fit.r.spe$coefficients[1]-50, yend = "Stenotrophomonas.maltophilia"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  annotate("text", label = "lnA_L estimated from tlag", size = 2,x = fit.r.spe$coefficients[1]+150, y = "Tetraselmis.tetrahele") +
  geom_segment(aes(x = fit.t.spe$coefficients[1], y = "Tetraselmis.tetrahele", xend = fit.r.spe$coefficients[1]+70, yend = "Tetraselmis.tetrahele"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  # add lnA_s(mean of E_species)
  geom_vline(xintercept=mean.a.r, linetype="dashed", color = "darkred") +
  geom_vline(xintercept=mean.a.t, linetype="dashed", color = "darkblue") +
  annotate("text", label = "lnA_S estimated from rmax",size = 2, x = mean.a.r-130, y = "Serratia.marcescens") +
  geom_segment(aes(x = mean.a.r, y = "Serratia.marcescens", xend = mean.a.r-50, yend = "Serratia.marcescens"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  annotate("text", label = "lnA_S estimated from tlag",size = 2, x = mean.a.t+150, y = "Oscillatoria.agardhii") +
  geom_segment(aes(x = mean.a.t, y = "Oscillatoria.agardhii", xend = mean.a.t+50, yend = "Oscillatoria.agardhii"), 
               arrow = arrow(length = unit(0.15, "cm")), colour = 'black') +
  # add CI of lnA_s
  annotate("rect", xmin=mean.a.r-mean.a.r.ci, xmax=mean.a.r+mean.a.r.ci,ymin=-Inf, ymax=Inf, alpha = .2) +
  annotate("rect", xmin=mean.a.t-mean.a.t.ci, xmax=mean.a.t+mean.a.t.ci,ymin=-Inf, ymax=Inf, alpha = .2) 
  
save.plot("../results/arrhenius/lnA_spe", p)




#### fit lm() get different parameter fitting in ggplot
# fittest <- lm(info0[info0$temp_c <= 30, ], formula =  log.r ~ temp_c)
# summary(fittest)
# ggplot(data = info0, aes(x = temp_c, y =log.r )) +
#   geom_point() +
#   stat_smooth(data = info0[info0$temp_c <= 30,], formula = y~x, method = lm) +
#   geom_ribbon(data = info0[info0$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), fill = 'grey', alpha = 0.5)+
#   stat_regline_equation(label.x = 10, label.y = -10)