setwd("/home/nesbit/Desktop/Data_RVK/Code")
rm(list=ls())
graphics.off()
yoot<-read.csv("../Data/All_Figs/LogisticGrowthData.csv")
yeet<-print(yoot)
head(yeet)
plot(yeet)
library(tidyverse)
library(ggplot2)
require("minpack.lm")
library(plot3D)
#install.packages("plot3D")

yiit<-read.csv("../Data/All_Figs/Bae_et_al_2014_processed.csv")

    RVK<-function(t, r, K, N0){
      N0*K*exp(r*t)/(K+N0*(exp(r*t)-1))
    }
    
head(yiit)
plot(yiit)
boop <- plot(yiit$Species, yiit$PopBio)
Data2Fit <- subset(yiit , Temp == "5")
PopBioLim <- c(0, 0.6)
scatter3D(Data2Fit$Time, Data2Fit$PopBio, Data2Fit$X, main= "Bae et al 3D representation of time, species, and popbio",
          ylab= "PopBio",
          xlab= "Time",
          zlab= "Species #")
# assembly %>%
# group_by(Time, Popbio)%>%
# ggplot(aes(Data2Fit, fill=Temp, color=Temp)) +
# geom_histogram(bins=5)+
# facet_wrap(~Temp + Run, nrow=11)

# PopFit<- nlsLM(PopBio ~ RVK(Time, R, K), data = df_subset, start = list(R = 0, K = 0))
# Lengths <- seq(min(df_subset$Time),max(df_subset$Time),len=1000)
# Predic2PlotPow <- RVK(Lengths,coef(PowFit)["a"],coef(PowFit)["b"])
f <-ggplot(Data2Fit, aes(x= Time, y= PopBio)) + geom_point(size = (1),color="red") + theme_bw() + 
  labs(y="PopBio", x = "Time (Hrs)") +
  facet_wrap(~Species)
# lines(Lengths, Predic2PlotPow, col = 'blue', lwd = 2.5)
plot(f)
head(f)
head(Data2Fit)
VTime <- seq(0, 650, length = 1000)
scale <- 4000

for(i in 1:22){
  #subset data
    df_subset <- subset(Data2Fit,Data2Fit$Species == unique(Data2Fit$Species)[i])
    Q = nrow(Data2Fit)

  #fitting curve
  PopFit<- nlsLM(PopBio ~ RVK(Time, r, K, N0), data = df_subset, 
                 start = list(N0= min(df_subset$PopBio), r = min(df_subset$PopBio), K = max(df_subset$PopBio)))

  Lengths <- seq(min(df_subset$Time),max(df_subset$Time),len=Q)
  Predic2PlotPow <- RVK(Lengths, coef(PopFit)["r"],coef(PopFit)["K"], coef(PopFit)["N0"])

    G<-ggplot(Data2Fit, aes(x= Time, y= PopBio)) + geom_point(size = (1),color="red") + theme_bw() + 
    labs(y="PopBio", x = "Time (Hrs)") +
    geom_line()+
    # geom_line(aes(nlsLM(PopBio ~ RVK(Time, r, K, N0), data = df_subset, 
    #                     start = list(N0= min(df_subset$PopBio), r = min(df_subset$PopBio), K = max(df_subset$PopBio)))))+
    facet_wrap(~Species)
   
  #Saving output
  print((G))
}
# df_subset <- subset(Data2Fit,Data2Fit$Species == unique(Data2Fit$Species)[i])
# Q = nrow(Data2Fit)
# PopFit<- nlsLM(PopBio ~ RVK(Time, r, K, N0), data = df_subset, 
#                start = list(N0= min(df_subset$PopBio), r = min(df_subset$PopBio), K = max(df_subset$PopBio)))
# 
# Lengths <- seq(min(df_subset$Time),max(df_subset$Time),len=Q)
# Predic2PlotPow <- RVK(Lengths, coef(PopFit)["r"],coef(PopFit)["K"], coef(PopFit)["N0"])
# 
# pred.log<-predict(PopFit, newdata = list(VTime=VTime))*scale
# G<-ggplot(Data2Fit, aes(x= Time, y= PopBio)) + geom_point(size = (1),color="red") + theme_bw() + 
#   labs(y="PopBio", x = "Time (Hrs)") +
#   geom_line(aes(y= pred.log))+
#   facet_wrap(~Species)
# print(G)
# df <- Data2Fit %>%
#   nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts) %>%
#   mutate(fits = map(data, ~ lm(.$Time, .$PopBio)))
# 
# 
#   mutate(r = unlist(map(fits,~ .$vals$r)),
#          r_p = unlist(map(fits,~ .$vals$r_p)),
#          K = unlist(map(fits,~ .$vals$k)),
#          K_p = unlist(map(fits,~ .$vals$k_p))) %>%
#   filter(r_p < .1 , K_p < .1) %>%
#   mutate(r = ifelse(r <= 0,1,r),
#          K = ifelse(K <= 0,1,K),
#          a = r/K)

# mistake_not<- print('Mistake Not My Current State Of Joshing Gentle Peevishness
#                     For The Awesome And Terrible Majesty Of The Towering Seas Of
#                     Ire That Are Themselves The Mere Milquetoast Shallows Fringing
#                     My Vast Oceans Of Wrath')
# 
# alb <- read.csv(file="../Data/All_Figs/Bae_et_al_2014_processed.csv")
# alb <- subset(x=alb, !is.na(alb$PopBio))
# plot(alb$Time, alb$PopBio, xlab="age (H)", ylab="PopBio", xlim=c(0, 700))
#     logistic1<-function(t, r, K, N0){
#       N0*K*exp(r*t)/(K+N0*(exp(r*t)-1))
#     }


