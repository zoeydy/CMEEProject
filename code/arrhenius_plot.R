


rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")

require(ggplot2)
library(ggpubr)

#############
# read data #
#############
infos <- read.csv("../data/gomp.info.csv")
# delete the non-positive parameters
info0 <- subset(infos, infos$rmax > 0 & infos$tlag > 0)
# set boltzmann constant 8.617*10^(-5) eV K^(-1)
K = 8.617*10^(-5)
# add 1/KT
info0$one.over.KT <- 1/(K*(info0$temp_k))

plot.df <- read.csv("../data/gomp.plot.csv")

Data <- read.csv('../data/mini_pop.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]

IDs <- unique(Data$id)


###################
# define function #
###################
# function to calculate the confidence interval
ci.calc <- function(x){
  MEAN <- mean(x)
  SD <- sd(x)
  n <- length(x)
  return(qnorm(0.975)*SD/sqrt(n))
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


# delete the non-positive parameters
info0 <- subset(infos, infos$rmax > 0 & infos$tlag > 0)


########
# plot #
########
##################################################
# log of 1/tlag and rmax VS temperature(celsius) #
##################################################
p <- ggplot(data = info0, aes(x = temp_c, y =log(rmax) )) +
  geom_point() +
  stat_smooth(data = info0[info0$temp_c <= 30,], formula = y~x, method = lm) +
  geom_ribbon(data = info0[info0$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), fill = 'grey', alpha = 0.5)+
  stat_regline_equation(label.x = 10, label.y = -10)

save.plot("../results/arrhenius/log_r_T", p)
p <- ggplot(data = info0, aes(x = temp_c, y =log(1/tlag) )) +
  geom_point() +
  stat_smooth(data = info0[info0$temp_c <= 30,], formula = y~x, method = lm) +
  geom_ribbon(data = info0[info0$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), fill = 'grey', alpha = 0.5)
save.plot("../results/arrhenius/log_1tlag_T", p)
# plot within species #
# calculate the info.spe data frame for plot by species with point more than 2 #
spe_info <- unique(infos$species)
info.spe <- data.frame()
for (i in 1:length(spe_info)) {
  dat <- infos[infos$species == spe_info[i],]
  if(nrow(dat) < 3){
    info.spe <- info.spe
  }else{
    info.spe <- rbind(info.spe, dat)
  }
}
# plot

p <- ggplot(info.spe, aes(x = temp_c, y = log(rmax), color = species, fill = species)) +
  geom_point() +
  stat_smooth(data = info.spe[info.spe$temp_c <= 30,] , formula = y~x, method = lm) +
  # geom_ribbon(data = info.spe[info.spe$temp_c >= 30,], aes(ymin = min(log(rmax)), ymax = max(log(rmax))), fill = 'grey', alpha = 0.5) +
  labs(title="log(rmax) VS temperature(°C) less than 30°C, group by species",x="Temperature(°C)", y = "log(rmax)") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~species)
save.plot("../results/arrhenius/spe_log_r_temp",p)


p <- ggplot(info.spe, aes(x = temp_c, y = log(1/tlag), color = species, fill = species)) +
  geom_point() +
  stat_smooth(data = info.spe[info.spe$temp_c <= 30,] , formula = y~x, method = lm) +
  # geom_ribbon(data = info.spe[info.spe$temp_c >= 30,], aes(ymin = min(log(1/tlag)), ymax = max(log(1/tlag))), fill = 'grey', alpha = 0.5) +
  labs(title="log(1/tlag) VS temperature(°C) less than 30°C, group by species",x="Temperature(°C)", y = "log(1/tlag)") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~species)
save.plot("../results/arrhenius/spe_log_tlag_temp",p)

##########################################################################################
# log of 1/tlag VS temperature(°C) less than 30°C in mean value with confidence interval #
########################################################################################## 
mean_df <- data.frame()
temp.c.s <- unique(info0$temp_c)
for (i in 1:length(temp.c.s)) {
  df1 <- subset(info0, info0$temp_c == temp.c.s[i])
  mean1 <- data.frame(temp.c = df1$temp_c[1],
                      mean.rmax = mean(log(df1$rmax)),
                      mean.1tlag = mean(log(1/df1$tlag)),
                      ci.rmax = ci.calc(log(df1$rmax)),
                      ci.1tlag = ci.calc(log(1/df1$tlag)))
  mean_df <- rbind(mean_df, mean1)
}
# rmax
r.temp <- mean_df[mean_df$mean.rmax == max(mean_df$mean.rmax),]$temp.c

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
save.plot("../results/arrhenius/mean_log_r_temp",p)

# tlag
t.temp <- 15
p <- ggplot(data = mean_df, aes(x = temp.c, y =mean.1tlag )) +
  geom_point() +
  stat_smooth(data = mean_df[mean_df$temp.c <= t.temp,],formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.1tlag - ci.1tlag, 
                    ymax=mean.1tlag + ci.1tlag,
                    width = 0.1)) +
  geom_ribbon(data = mean_df[mean_df$temp.c >= t.temp, ], aes(ymin = -Inf, ymax = Inf, alpha = 0.1)) +
  theme(legend.position = 'none') +
  labs(title = "Mean of log(1/tlag) VS temperature(°C) less than T_opt", x = "Temperature(°C)",y = "Mean of log(1/tlag) +/- CI")
save.plot("../results/arrhenius/mean_log_tlag_temp",p)



#############################################################################################
# log of rmax and 1/tlag VS temperature(°C) less than T_opt in different temperature groups #
#############################################################################################
# devide the temperature into 4 groups ((-5:10, 10:20, 20:30, 30:40))
df1 <- subset(info0, info0$temp_c >= -5 & info0$temp_c < 10)
df2 <- subset(info0, info0$temp_c >= 10 & info0$temp_c < 20)
df3 <- subset(info0, info0$temp_c >= 20 & info0$temp_c < 30)
df4 <- subset(info0, info0$temp_c >= 30 & info0$temp_c < 40)
df1$temp_c_group <- "-5~10 °C"
df2$temp_c_group <- "10~20 °C"
df3$temp_c_group <- "20~30 °C"
df4$temp_c_group <- "30~40 °C"
info_4 <- rbind(df1,df2,df3,df4)
# devide the temperature into 5 groups ((-5:5, 5:10, 10:20, 20:30, 30:40))
dm1 <- subset(info0, info0$temp_c >= -5 & info0$temp_c < 5)
dm2 <- subset(info0, info0$temp_c >= 5 & info0$temp_c < 10)
dm3 <- subset(info0, info0$temp_c >= 10 & info0$temp_c < 20)
dm4 <- subset(info0, info0$temp_c >= 20 & info0$temp_c < 30)
dm5 <- subset(info0, info0$temp_c >= 30 & info0$temp_c < 40)
dm1$temp_c_group <- "-5~5 °C"
dm2$temp_c_group <- "5~10 °C"
dm3$temp_c_group <- "10~20 °C"
dm4$temp_c_group <- "20~30 °C"
dm5$temp_c_group <- "30~40 °C"
info_5 <- rbind(dm1,dm2,dm3,dm4,dm5)
info_5$temp_c_group <- factor(info_5$temp_c_group, levels = unique(info_5$temp_c_group))
# devide the temperature into 5 groups ((-5:0, 0:8, 8:16, 16:32, 32:40))
dmm1 <- subset(info0, info0$temp_c >= -5 & info0$temp_c < 0)
dmm2 <- subset(info0, info0$temp_c >= 0 & info0$temp_c < 8)
dmm3 <- subset(info0, info0$temp_c >= 8 & info0$temp_c < 16)
dmm4 <- subset(info0, info0$temp_c >= 16 & info0$temp_c < 24)
dmm5 <- subset(info0, info0$temp_c >= 24 & info0$temp_c < 32)
dmm6 <- subset(info0, info0$temp_c >= 32 & info0$temp_c < 40)
dmm1$temp_c_group <- "-5~0 °C"
dmm2$temp_c_group <- "0~8 °C"
dmm3$temp_c_group <- "8~16 °C"
dmm4$temp_c_group <- "16~24 °C"
dmm5$temp_c_group <- "24~32 °C"
dmm6$temp_c_group <- "32~40 °C"
info_6 <- rbind(dmm1,dmm2,dmm3,dmm4,dmm5)
info_6$temp_c_group <- factor(info_6$temp_c_group, levels = unique(info_6$temp_c_group))


# function to get dataframe by grouping info0 into specific groups by temp_c_group 
group_temp <- function(dat, param){
  temp <- unique(dat$temp_c_group)
  temp.df <- data.frame()
  for (i in 1:length(temp)) {
    info.temp <- subset(dat, dat$temp_c_group == temp[i])
    if (param == "rmax") {
      temp.per <- data.frame(MEAN = mean(log(info.temp[,param])), 
                             CI = ci.calc(log(info.temp[,param])), 
                             MEAN.TEMP = mean(info.temp$temp_c),
                             GROUP = info.temp$temp_c_group[1])
    }
    else{
      temp.per <- data.frame(MEAN = mean(log(1/info.temp[,param])),
                             CI = ci.calc(log(1/info.temp[,param])),
                             MEAN.TEMP = mean(info.temp$temp_c),
                             GROUP = info.temp$temp_c_group[1])
    }
    temp.df <- rbind(temp.df, temp.per)
  }
  return(temp.df)
}

# function to plot grouped data frame with CI
group_plot <- function(dat, param, lab){
  temp.max <- dat[dat$MEAN == max(dat$MEAN), ]$MEAN.TEMP
  ggplot(data = dat, aes(x = MEAN.TEMP, y = MEAN)) +
    geom_point() +
    stat_smooth(data = dat[dat$MEAN.TEMP <= temp.max, ],formula = y~x, method = lm, se = TRUE)+
    geom_errorbar(aes(ymin = MEAN - CI, ymax = MEAN + CI, width = .1)) +
    # labs(x="Temperature Range", y="MEAN (1/h) +/- CI ") +
    labs(y="MEAN (1/h) +/- CI ") +
    geom_ribbon(data = dat[dat$MEAN.TEMP >= temp.max,],
                aes(ymin = -Inf, ymax = Inf, alpha = .3)) +
    theme(legend.position = 'none') +
    scale_x_continuous(name = "Temperature Range", labels = dat$GROUP, breaks=dat$MEAN.TEMP) +
    labs(tag = lab)
}

# group plot rate VS temp_group
g1 <- group_plot(group_temp(info_4,"tlag"), "tlag", "A")
g2 <- group_plot(group_temp(info_4,"rmax"), "rmax", "B")
g3 <- group_plot(group_temp(info_5,"tlag"), "tlag", "C")
g4 <- group_plot(group_temp(info_5,"rmax"), "rmax", "D")
g5 <- group_plot(group_temp(info_6,"tlag"), "tlag", "E")
g6 <- group_plot(group_temp(info_6,"rmax"), "rmax", "F")
# left is tlag, right is rmax
pdf("../results/arrhenius/mean_rate_temp_group.pdf")
p <- grid.arrange(
  g1,g2,g3,g4,g5,g6,
  #labels = LETTERS[1:6],
  nrow = 3,
  top = "Across species fitness value VS Temperature(Grouped)",
  bottom = textGrob(
    "CI is calculated from normal distribution",
    gp = gpar(fontface = 6, fontsize = 9),
    hjust = 1,
    x = 1
  )
)
print(p)
graphics.off()


######################################
# Arrhenius plot of everything above #
######################################

##################################################
# log of 1/tlag and rmax VS 1/KT #
##################################################
# function to plot y parameters VS 1/temp
p <- ggplot(data = info0, aes(x = one.over.KT, y =log(rmax) )) +
  geom_point() +
  stat_smooth(data = info0[info0$temp_c <= 30,], formula = y~x, method = lm) +
  geom_ribbon(data = info0[info0$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), fill = 'grey', alpha = 0.5) +
  theme(legend.position = 'none')
save.plot("../results/arrhenius/log_r_KT",p)

p <- ggplot(data = info0, aes(x = one.over.KT, y =log(1/tlag) )) +
  geom_point() +
  stat_smooth(data = info0[info0$temp_c <= 30,], formula = y~x, method = lm) +
  geom_ribbon(data = info0[info0$temp_c >= 30,], aes(ymin = -Inf, ymax = Inf), fill = 'grey', alpha = 0.5) +
  theme(legend.position = 'none')
save.plot("../results/arrhenius/log_1tlag_KT",p)

##############################################################
# log of 1/tlag VS 1/KT less than T_pk in mean value with CI #
##############################################################
mean_df <- data.frame()
temp.c.s <- unique(info0$temp_c)
for (i in 1:length(temp.c.s)) {
  df1 <- subset(info0, info0$temp_c == temp.c.s[i])
  mean1 <- data.frame(temp.c = df1$temp_c[1],
                      mean.rmax = mean(log(df1$rmax)),
                      mean.1tlag = mean(log(1/df1$tlag)),
                      ci.rmax = ci.calc(log(df1$rmax)),
                      ci.1tlag = ci.calc(log(1/df1$tlag)))
  mean_df <- rbind(mean_df, mean1)
}
# rmax
r.temp <- mean_df[mean_df$mean.rmax == max(mean_df$mean.rmax),]$temp.c
p <- ggplot(data = mean_df, aes(x = one.over.KT, y =mean.rmax )) +
  geom_point() +
  stat_smooth(data = mean_df[mean_df$temp.c <= r.temp,],
              formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.rmax - ci.rmax,
                    ymax=mean.rmax + ci.rmax,
                    width = 0.1)) +
  geom_ribbon(data = mean_df[mean_df$temp.c >= r.temp,], aes(ymin=-Inf, ymax=Inf, alpha = .5 )) +
  theme(legend.position = 'none') +
  labs(title = "Mean of log(rmax) VS temperature(°C) less than T_opt", x = "Temperature(°C)", y = "Mean of log(rmax) +/- CI")
save.plot("../results/arrhenius/mean_log_r_KT", p)

# tlag
t.temp <- mean_df[mean_df$mean.1tlag == max(mean_df$mean.1tlag),]$temp.c
p <- ggplot(data = mean_df, aes(x = one.over.KT, y =mean.1tlag )) +
  geom_point() +
  stat_smooth(data = mean_df[mean_df$temp.c <= t.temp,],formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.1tlag - ci.1tlag,
                    ymax=mean.1tlag + ci.1tlag,
                    width = 0.1)) +
  geom_ribbon(data = mean_df[mean_df$temp.c >= t.temp, ], aes(ymin = -Inf, ymax = Inf, alpha = 0.1)) +
  theme(legend.position = 'none') +
  labs(title = "Mean of log(1/tlag) VS 1/KT", x = "1/KT",y = "Mean of log(1/tlag) +/- CI")
save.plot("../results/arrhenius/mean_log_tlag_KT",p)


