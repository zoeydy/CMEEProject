


rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")

require(ggplot2)
library(ggpubr)

# set boltzmann constant 1.38064852 × 10-23 m2 kg s-2 K-1
K = 1.38064852 * 10^(-23)


# read the data, starting value and compare models
infos <- read.csv("../data/gomp.info.csv")
plot.df <- read.csv("../data/gomp.plot.csv")

Data <- read.csv('../data/mini_pop.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]

IDs <- unique(Data$id)
spes <- unique(Data$Species)
temps <- unique(info$temp_group)

###################
# define function #
###################
# function to calculate standard error
se <- function(x){
  sqrt(var(x))/length(x)
}
# function to get dataframe by grouping info_0 into specific groups by temp_c_group 
group_temp <- function(dat, param){
  temp <- unique(dat$temp_c_group)
  temp.df <- data.frame()
  for (i in 1:length(temp)) {
    info.temp <- subset(dat, dat$temp_c_group == temp[i])
    if (param == "rmax") {
      temp.per <- data.frame(MEAN = mean(log(info.temp[,param])), 
                             SE = se(log(info.temp[,param])), 
                             GROUP = info.temp$temp_c_group[1])
    }
    else{
      temp.per <- data.frame(MEAN = mean(log(1/info.temp[,param])),
                             SE = se(log(1/info.temp[,param])),
                             GROUP = info.temp$temp_c_group[1])
    }
    temp.df <- rbind(temp.df, temp.per)
  }
  return(temp.df)
}


# function to plot grouped data frame with SE
group_plot <- function(dat, param){
  ggplot(data = dat, aes(x = GROUP, y = MEAN)) +
    geom_point() +
    stat_smooth(formula = y~x, method = lm, se = TRUE)+
    geom_errorbar(aes(ymin = MEAN - SE, ymax = MEAN + SE, width = .1)) +
    labs(title=paste0("Rate(measured by )", param, " VS Temperature"), 
         x="Temperature Group", 
         y="MEAN rate +/- SE (1/h)")
}
# group_plot(group_temp(info_4, "tlag"), "tlag")

# function to plot y parameters VS 1/temp
arrhe_one <- function(dat, param){
  if(param == "tlag"){
    ggplot(data = dat, aes(x = 1/temp_k, y = log(1/dat[,param]))) +
      geom_point() +
      stat_smooth(method = lm, formula = y~x)
  }else{
    ggplot(data = dat, aes(x = 1/temp_k, y = log(dat[,param]))) +
      geom_point() +
      stat_smooth(method = lm, formula = y~x) 
  }
}


# function to save plot
save_plot <- function(path, plotfun){
  pdf(path)
  print(plotfun)
  graphics.off()
}


############################
# get complete information #
############################
info <- data.frame()

for (i in 1:length(IDs)){
  # subset info and data by idname
  idname <- IDs[i]
  data <- Data[Data$id == idname,]
  # plot.id <- plot.df[plot.df$id == idname, ]
  info.id <- infos[infos$id == idname, ]
  # get info about temperature, species and medium
  info.id$temp_c <- data$Temp[1]
  info.id$species <- data$Species[1]
  info.id$medium <- data$Medium[1]
  # add kelvin temperature
  info.id$temp_k <- info.id$temp_c + 273
  # rbind the data frame
  info <- rbind(info, info.id)
}

# delete the non-positive parameters
info_0 <- subset(info, info$rmax > 0 & info$tlag > 0)
# devide the temperature into 4 groups ((-5:10, 10:20, 20:30, 30:40))
df1 <- subset(info_0, info_0$temp_c >= -5 & info_0$temp_c < 10)
df2 <- subset(info_0, info_0$temp_c >= 10 & info_0$temp_c < 20)
df3 <- subset(info_0, info_0$temp_c >= 20 & info_0$temp_c < 30)
df4 <- subset(info_0, info_0$temp_c >= 30 & info_0$temp_c < 40)
df1$temp_c_group <- "-5~10 °C"
df2$temp_c_group <- "10~20 °C"
df3$temp_c_group <- "20~30 °C"
df4$temp_c_group <- "30~40 °C"
info_4 <- rbind(df1,df2,df3,df4)
# devide the temperature into 5 groups ((-5:5, 5:10, 10:20, 20:30, 30:40))
dm1 <- subset(info_0, info_0$temp_c >= -5 & info_0$temp_c < 5)
dm2 <- subset(info_0, info_0$temp_c >= 5 & info_0$temp_c < 10)
dm3 <- subset(info_0, info_0$temp_c >= 10 & info_0$temp_c < 20)
dm4 <- subset(info_0, info_0$temp_c >= 20 & info_0$temp_c < 30)
dm5 <- subset(info_0, info_0$temp_c >= 30 & info_0$temp_c < 40)
dm1$temp_c_group <- "-5~5 °C"
dm2$temp_c_group <- "5~10 °C"
dm3$temp_c_group <- "10~20 °C"
dm4$temp_c_group <- "20~30 °C"
dm5$temp_c_group <- "30~40 °C"
info_5 <- rbind(dm1,dm2,dm3,dm4,dm5)
info_5$temp_c_group <- factor(info_5$temp_c_group, levels = unique(info_5$temp_c_group))


########
# plot #
########

# log of 1/tlag and rmax VS temperature(celsius)
ggplot(data = info_0, aes(x = temp_c, y =log(rmax) )) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm)
ggplot(data = info_0, aes(x = temp_c, y =log(1/tlag) )) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm) 

# log of 1/tlag and rmax VS temperature(celsius) in mean value with error bar
mean_df <- data.frame()
temp.c.s <- unique(info_0$temp_c)
for (i in 1:length(temp.c.s)) {
  df1 <- subset(info_0, info_0$temp_c == temp.c.s[i])
  mean1 <- data.frame(temp.c = df1$temp_c[1],
                      mean.rmax = mean(log(df1$rmax)),
                      mean.1tlag = mean(log(1/df1$tlag)),
                      se.rmax = se(log(df1$rmax)),
                      se.1tlag = se(log(1/df1$tlag)))
  mean_df <- rbind(mean_df, mean1)
}
ggplot(data = mean_df, aes(x = temp.c, y =mean.rmax )) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.rmax - se.rmax, 
                    ymax=mean.rmax + se.rmax,
                    width = 0.1))
ggplot(data = mean_df, aes(x = temp.c, y =mean.1tlag )) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm) +
  geom_errorbar(aes(ymin=mean.1tlag - se.1tlag, 
                    ymax=mean.1tlag + se.1tlag,
                    width = 0.1))

# group plot rate VS temp_group
group_plot(group_temp(info_4,"tlag"), "tlag")
group_plot(group_temp(info_4,"rmax"), "rmax")
group_plot(group_temp(info_5,"tlag"), "tlag")
group_plot(group_temp(info_5,"rmax"), "rmax")

# Arrhenius plot of 1 parameter: rate VS 1/temp
arrhe_one(info_0, "tlag")
arrhe_one(info_0, "rmax")

# Arrhenius plot of 2 parameters: ramx*tlag VS 1/temp
ggplot(data = info_0, aes(x = 1/temp_k, y = log(rmax*tlag))) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm) +
  geom_text(x = 0.003, y = -10, label = lm_eqn(info_0, x = 1/info_0$temp_k, y = log(info_0$rmax*info_0$tlag)), parse = TRUE)
lm_eqn <- function(df, x, y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

