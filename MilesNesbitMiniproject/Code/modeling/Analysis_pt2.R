require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)
library(dplyr)
library(growthcurver)
library(stringr)
library(taxize)
library(nls.multstart)
library(zoo)

rm(list=ls())
graphics.off()

setwd("/home/nesbit/Desktop/Data_RVK/Code")
files <- read.csv('../Results/Formated/Total_processed.csv')

#clean
# drops <- c("X", "X.1")
# files <- files[ , !(names(files) %in% drops)]


#subsetting absorbance
Absorbance <- files %>%
  filter(str_detect(PopBio_unts, 'Abs'))

ABSD <- files %>%
  filter(str_detect(PopBio_unts, 'OD'))

Absorbance <- rbind(Absorbance, ABSD)

#trying through program
ff <- Absorbance %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio")) %>%
  mutate(fits = map(data, ~ SummarizeGrowth(.$Time, .$PopBio))) %>%
  mutate(r = unlist(map(fits,~ .$vals$r)),
         r_p = unlist(map(fits,~ .$vals$r_p)),
         K = unlist(map(fits,~ .$vals$k)),
         K_p = unlist(map(fits,~ .$vals$k_p))) %>%
  filter(r_p < .1 , K_p < .1) %>%
  mutate(r = ifelse(r <= 0,1,r),
         K = ifelse(K <= 0,1,K),
         a = r/K)

#plot
ggplot(ff, aes(x= log(r), y= log(K), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Absorbance")

#try linear?
lin_RVK_Abs <- lm(log(r) ~ log(K), data = ff)
summary(lin_RVK_Abs)

ggplot(ff, aes(x= log(r), y= log(K), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Absorbance")+
  geom_smooth(method='lm', formula= y~x)

#CFU
CFU <- files %>%
  filter(str_detect(PopBio_unts, 'CFU'))

#trying through program
gf <- CFU %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio")) %>%
  mutate(fits = map(data, ~ SummarizeGrowth(.$Time, .$PopBio))) %>%
  mutate(r = unlist(map(fits,~ .$vals$r)),
         r_p = unlist(map(fits,~ .$vals$r_p)),
         K = unlist(map(fits,~ .$vals$k)),
         K_p = unlist(map(fits,~ .$vals$k_p))) %>%
  filter(r_p < .1 , K_p < .1) %>%
  mutate(r = ifelse(r <= 0,1,r),
         K = ifelse(K <= 0,1,K),
         a = r/K)

#plot
ggplot(gf, aes(x= log(r), y= log(K), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Colony Forming Units")

lin_RVK_CFU <- lm(log(r) ~ log(K), data = gf)
summary(lin_RVK_CFU)

length(unique(files$Citation))

length(unique(files$PopBio_unts))



#N
N <- files %>%
  filter(str_detect(PopBio_unts, 'N'))


#growth curver
nf <- N %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio")) %>%
  mutate(fits = map(data, ~ SummarizeGrowth(.$Time, .$PopBio))) %>%
  mutate(r = unlist(map(fits,~ .$vals$r)),
         r_p = unlist(map(fits,~ .$vals$r_p)),
         K = unlist(map(fits,~ .$vals$k)),
         K_p = unlist(map(fits,~ .$vals$k_p))) %>%
  filter(r_p < .1 , K_p < .1) %>%
  mutate(r = ifelse(r <= 0,1,r),
         K = ifelse(K <= 0,1,K),
         a = r/K)

#plot
ggplot(nf, aes(x= log(r), y= log(K), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("N")

#N and CFU = same?
N <- files %>%
  filter(str_detect(PopBio_unts, 'N'))

NCFU <- rbind(N, CFU)

#growth curver
ncfuf <- NCFU %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio")) %>%
  mutate(fits = map(data, ~ SummarizeGrowth(.$Time, .$PopBio))) %>%
  mutate(r = unlist(map(fits,~ .$vals$r)),
         r_p = unlist(map(fits,~ .$vals$r_p)),
         K = unlist(map(fits,~ .$vals$k)),
         K_p = unlist(map(fits,~ .$vals$k_p))) %>%
  filter(r_p < .1 , K_p < .1) %>%
  mutate(r = ifelse(r <= 0,1,r),
         K = ifelse(K <= 0,1,K),
         a = r/K)

#plot
ggplot(ncfuf, aes(x= log(r), y= log(K), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("NCFU")+
  facet_wrap(~Species)



#subset all the way down
smith_cb <-subset(subset(subset(files, Species == 'Serratia'), Rep == 1), Temp == 30)
drops <- c("Citation", "Temp", 'Species', 'Time_units', 'PopBio_unts', 'Medium', 'Rep')
smith_cb <- smith_cb[ , !(names(smith_cb) %in% drops)]


#get models
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

#fit logistic
data <- smith_cb
data <- data %>% 
  rename(
    't'=Time,
   'N' =PopBio
  )


data$LogN <- log(data$N)
data$Log10N <- log10(data$N) # for later

# visualise
ggplot(data, aes(x = t, y = LogN)) + 
  geom_point(size = 3) +
  labs(x = "Time (Hours)", y = "log(Absorbance)") +
  theme_bw(base_size = 16)

#draw a straight line through the linear part on log-transformed data.
lm_growth <- lm(LogN ~ t, data = data[data$t > 250 & data$t < 1000,])
summary(lm_growth)



#nls
logistic_model <- function(t, r_max, N_max, N_0){ # The classic logistic equation
  return(N_0 * N_max * exp(r_max * t)/(N_max + N_0 * (exp(r_max * t) - 1)))
}

# first we need some starting parameters for the model
N_0_start <- min(data$N) # lowest population size
N_max_start <- max(data$N) # highest population size
r_max_start <- max(diff(data$N))/mean(diff(data$t))

fit_logistic <- nlsLM(N ~ logistic_model(t = t, r_max, N_max, N_0), data,
                      list(r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

summary(fit_logistic)

# plot it:
timepoints <- seq(0, 1250, 1)

logistic_points <- logistic_model(t = timepoints, r_max = coef(fit_logistic)["r_max"], N_max = coef(fit_logistic)["N_max"], N_0 = coef(fit_logistic)["N_0"])
df1 <- data.frame(timepoints, logistic_points)
df1$model <- "Logistic"
names(df1) <- c("t", "N", "model")

ggplot(data, aes(x = t, y = N)) +
  geom_point(size = 3) +
  geom_line(data = df1, aes(x = t, y = N, col = model), size = 1) +
  theme_bw(base_size = 16) + # make the background white
  theme(aspect.ratio=1)+ # make the plot square 
  labs(x = "Time", y = "Absorbance")

# Note that we're using Log10N now
N_0_start <- min(data$Log10N)
N_max_start <- max(data$Log10N)
t_lag_start <- data$t[which.max(diff(diff(data$Log10N)))]
r_max_start <- max(diff(data$Log10N))/mean(diff(data$t))


#fit the models
fit_baranyi <- nlsLM(Log10N ~ baranyi_model(t = t, r_max, N_max, N_0, t_lag), data,
                     list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_buchanan <- nlsLM(Log10N ~ buchanan_model(t = t, r_max, N_max, N_0, t_lag), data,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_gompertz <- nlsLM(Log10N ~ gompertz_model(t = t, r_max, N_max, N_0, t_lag), data,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

summary(fit_baranyi)
summary(fit_buchanan)
summary(fit_gompertz)


timepoints <- seq(0, 1250, 1)

baranyi_points <- baranyi_model(t = timepoints, r_max = coef(fit_baranyi)["r_max"], N_max = coef(fit_baranyi)["N_max"], N_0 = coef(fit_baranyi)["N_0"], t_lag = coef(fit_baranyi)["t_lag"])

buchanan_points <- buchanan_model(t = timepoints, r_max = coef(fit_buchanan)["r_max"], N_max = coef(fit_buchanan)["N_max"], N_0 = coef(fit_buchanan)["N_0"], t_lag = coef(fit_buchanan)["t_lag"])

gompertz_points <- gompertz_model(t = timepoints, r_max = coef(fit_gompertz)["r_max"], N_max = coef(fit_gompertz)["N_max"], N_0 = coef(fit_gompertz)["N_0"], t_lag = coef(fit_gompertz)["t_lag"])

df2 <- data.frame(timepoints, baranyi_points)
df2$model <- "Baranyi"
names(df2) <- c("t", "Log10N", "model")

df3 <- data.frame(timepoints, buchanan_points)
df3$model <- "Buchanan"
names(df3) <- c("t", "Log10N", "model")

df4 <- data.frame(timepoints, gompertz_points)
df4$model <- "Gompertz"
names(df4) <- c("t", "Log10N", "model")

model_frame <- rbind(df2, df3, df4)

ggplot(data, aes(x = t, y = Log10N)) +
  geom_point(size = 3) +
  geom_line(data = model_frame, aes(x = t, y = Log10N, col = model), size = 1) +
  theme_bw(base_size  = 16) + 
  theme(aspect.ratio=1) + 
  labs(x = "Time (Hours)", y = "log10(Absorbance)", title = 'Example: Smith, Serratia, 10% LB broth, 30C')


#selection
whoop <- AIC(fit_gompertz, fit_baranyi, fit_buchanan)

#check r and k values
whoop2 <- confint(fit_baranyi)
