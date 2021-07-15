# install.packages("devtools")
# require(devtools)
# install_version("cowplot", version = "0.9.0", repos = "http://cran.us.r-project.org")
library(tidyverse)
library(broom)
library(growthcurver)
library(cowplot)
library(psych)
library(MASS)
library(fitdistrplus)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
files <- read.csv('../Results/Formated/Total_processed.csv')

data_list <- list()

final_df <- files

#clean it
drops <- c("X", "X.1")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#plotting series
final_df %>%
  ggplot(aes(x=Time,y=PopBio,colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~paste(Species,Medium), scales = "free")

#fit growth curves
df <- final_df %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio")) %>%
  mutate(fits = map(data, ~ SummarizeGrowth(.$Time, .$PopBio))) %>%
  mutate(r = unlist(map(fits,~ .$vals$r)),
         r_p = unlist(map(fits,~ .$vals$r_p)),
         K = unlist(map(fits,~ .$vals$k)),
         K_p = unlist(map(fits,~ .$vals$k_p))) %>%
  #filter(r_p < .1 , K_p < .1) %>%
  mutate(r = ifelse(r <= 0,1,r),
         K = ifelse(K <= 0,1,K),
         a = r/K)

#look at some pretty dots
plot(log(df$r),log(df$K))

#fit boltzmann curves
k <- 8.617e-5

params_df <- df %>%
  dplyr::select(-data,-fits) %>% #remove extra data
  mutate(a = r / K) %>% #calculate a_ii
  gather("Param","Val",r,K,a) %>% #combine r and K to long format
  mutate(Temp = (1 / (k * (Temp+273.15))) - (1 / (k * (15+273.15))) ,Val = log(Val)) %>% #transform for boltz
  #filter(!is.infinite(Val)) %>% #remove inf
  group_by(Citation,Param,Species,Medium,Rep) %>% mutate(n = n()) %>% ungroup() %>%
  filter(n > 2) %>%
  nest(-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,-Param) %>% #nest for boltz fitting
  mutate(fits = map(data, ~lm(.$Val ~ .$Temp)), #fit models
         tidied = map(fits,tidy)) %>% #tidy output
  unnest(tidied) %>% #unnest
  mutate(term = recode(term,`(Intercept)` = "B0", `.$Temp` = "E")) %>% #recode the parameter estimates
  filter(p.value < .5) %>% #cut-off for p-values
  dplyr::select(-std.error,-statistic,-p.value) %>%
  spread(term,estimate)  %>%
  filter(!is.na(B0),!is.na(E))

p1 <- params_df %>%
  ggplot(aes(x=E,group = PopBio_unts, colour = PopBio_unts))+
  geom_density()+
  facet_wrap(~Param)+
  geom_vline(xintercept = 0)

p2 <- params_df %>%
  ggplot(aes(x=B0,group = PopBio_unts, colour = PopBio_unts))+
  geom_density()+
  facet_wrap(~Param)+
  geom_vline(xintercept = 0)

#distributions
plot_grid(p1,p2,nrow = 2)


#lognormal skew
ln_skew <- function(v){
  (exp(v)+2) * sqrt(exp(v)-1)
}

params_df %>%
  group_by(Param,PopBio_unts) %>%
  summarise(n = n(),
            uE = mean(E),vE=var(E),kE=skew(E),
            uB0= mean(B0), vB0 = var(B0),kB0 = skew(B0)) %>%
  gather("measure","value",uE,vE,kE,uB0,vB0,kB0) %>%
  ggplot(aes(x=interaction(measure,Param),y=value,colour=PopBio_unts))+
  geom_point()

#r
x <- params_df %>% filter(Param == "a")

fit_B0 <- fitdistr(x$B0,"normal")
B0_pred = rnorm(1000,fit_B0$estimate["mean"],fit_B0$estimate["sd"])

fit_E <- fitdistr(x$E,"normal")
E_pred = rnorm(1000,fit_E$estimate["mean"],fit_E$estimate["sd"])

p1 <- data.frame(B0 = x$B0) %>%
  ggplot(aes(x=B0))+
  geom_density(colour = "red") +
  geom_density(data = data.frame(B0 = B0_pred),colour = "blue")

p2 <- data.frame(E = x$E) %>%
  ggplot(aes(x=E))+
  geom_density(colour = "red") +
  geom_density(data = data.frame(E = E_pred),colour = "blue")

plot_grid(p1,p2)

descdist(x$E,discrete=FALSE, boot=500)


#specify the model functions:
logistic_model <- function(t, r_max, N_max, N_0){ # The classic logistic equation
  return(N_0 * N_max * exp(r_max * t)/(N_max + N_0 * (exp(r_max * t) - 1)))
}

gompertz_model <- function(t, r_max, N_max, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}

baranyi_model <- function(t, r_max, N_max, N_0, t_lag){  # Baranyi model (Baranyi 1993)
  return(N_max + log10((-1+exp(r_max*t_lag) + exp(r_max*t))/(exp(r_max*t) - 1 + exp(r_max*t_lag) * 10^(N_max-N_0))))
}

buchanan_model <- function(t, r_max, N_max, N_0, t_lag){ # Buchanan model - three phase logistic (Buchanan 1997)
  return(N_0 + (t >= t_lag) * (t <= (t_lag + (N_max - N_0) * log(10)/r_max)) * r_max * (t - t_lag)/log(10) + (t >= t_lag) * (t > (t_lag + (N_max - N_0) * log(10)/r_max)) * (N_max - N_0))
}


