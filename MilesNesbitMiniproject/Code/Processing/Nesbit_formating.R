library(tidyverse)
library(ggplot2)
library(plyr)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")

Total_start<-read.csv('../Results/All_Figs/LogisticGrowthData.csv')

a<-read.csv('../Results/Formated/Bae_et_al_2014_processed.csv')
b<-read.csv('../Results/Formated/Bernhard_et_al_2018_processed.csv')
c<-read.csv('../Results/Formated/Blagodatskaya_processed.csv')
d<-read.csv('../Results/Formated/Galarz_et_al_2015_processed.csv')
e<-read.csv('../Results/Formated/Gill_and_DeLacey_1991_processed.csv')
f<-read.csv('../Results/Formated/Heo_et_al_processed.csv')
g<-read.csv('../Results/Formated/Inoue_1977_processed.csv')
h<-read.csv('../Results/Formated/Kapetanakou_processed.csv')
i<-read.csv('../Results/Formated/Kirchman_Rich_1997_processed.csv')
j<-read.csv('../Results/Formated/Koutsoumanis_processed.csv')
k<-read.csv('../Results/Formated/Lee_processed.csv')
l<-read.csv('../Results/Formated/Phillips_and_Griffiths_1987_processed.csv')
# m<-read.csv('../Results/Formated/Ratkowsky***_processed.csv')
n<-read.csv('../Results/Formated/Roth_and_Wheaton_1961_processed.csv')
o<-read.csv('../Results/Formated/Sato-Takabe_processed.csv')
p<-read.csv('../Results/Formated/Silva_et_al_2018_processed.csv')
q<-read.csv('../Results/Formated/Sivonen_1990_processed.csv')
r<-read.csv('../Results/Formated/Stannard_et_al_1985_processed.csv')
s<-read.csv('../Results/Formated/Vankerschaver_processed.csv')
t<-read.csv('../Results/Formated/Wilocx_processed.csv')
u<-read.csv('../Results/Formated/Zwietering_1994_processed.csv')
v<-read.csv('../Results/Formated/Smith_processed.csv')

#bind all together
Total_final <- rbind.fill(a,b,c,d,e,f,g,h,i,j,k,l,n,o,p,q,r,s,t,u,v)

#fix stuff
Total_final <- Total_final %>%
  mutate(Temp = as.numeric(Temp))+
  mutate(PopBio = as.numeric(PopBio))+
  mutate(Time = as.numeric(Time))+
  mutate(Rep = as.numeric(Rep))

write.csv(Total_final,"../Results/Formated/Total_processed.csv", row.names = FALSE)

#total added 1583 lines
