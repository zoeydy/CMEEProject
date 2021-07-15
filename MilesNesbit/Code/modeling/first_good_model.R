setwd("/home/nesbit/Desktop/bae_cb_RVK/Code")
require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)



rm(list = ls())
graphics.off()
final_df <- read.csv('../Results/Formated/Total_processed.csv')

#clean it
drops <- c("X", "X.1")
final_df <- final_df[ , !(names(final_df) %in% drops)]

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

#subsetting
bae_cb <-subset(subset(final_df, Temp == 5), Species == 'Chryseobacterium.balustinum')
drops <- c("Citation", "Temp", 'Species', 'Time_units', 'PopBio_unts', 'Medium', 'Rep')
bae_cb <- bae_cb[ , !(names(bae_cb) %in% drops)]

#sanity check
length(unique(bae_cb$Species))

#test plot
bae_cb %>%
  ggplot(aes(x=Time,y=PopBio))+
  geom_point()+
  geom_line()

#starting bae_cb
N_0_start <- min(bae_cb$PopBio)
N_max_start <- max(bae_cb$PopBio)
t_lag_start <- bae_cb$Time[which.max(diff(diff(bae_cb$PopBio)))]
r_max_start <- max(diff(bae_cb$PopBio))/mean(diff(bae_cb$Time))


#fit models
fit_logistic <- nlsLM(PopBio ~ logistic_model(t = Time, r_max, N_max, N_0), bae_cb,
                      list(r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_baranyi <- nlsLM(PopBio ~ baranyi_model(t = Time, r_max, N_max, N_0, t_lag), bae_cb,
                     list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

#didn't work
fit_buchanan <- nlsLM(PopBio ~ buchanan_model(t = Time, r_max, N_max, N_0, t_lag), bae_cb,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_gompertz <- nlsLM(PopBio ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag), bae_cb,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

summary(fit_logistic)
summary(fit_baranyi)
# summary(fit_buchanan) #didn't work
summary(fit_gompertz)

#visualize
timepoints <- seq(0, 6.698795e+02, length.out = 2000)
min(timepoints)
max(timepoints)

logistic_points <- logistic_model(t = timepoints, r_max = coef(fit_logistic)["r_max"], N_max = coef(fit_logistic)["N_max"], N_0 = coef(fit_logistic)["N_0"])

baranyi_points <- baranyi_model(t = timepoints, r_max = coef(fit_baranyi)["r_max"], N_max = coef(fit_baranyi)["N_max"], N_0 = coef(fit_baranyi)["N_0"], t_lag = coef(fit_baranyi)["t_lag"])

buchanan_points <- buchanan_model(t = timepoints, r_max = coef(fit_buchanan)["r_max"], N_max = coef(fit_buchanan)["N_max"], N_0 = coef(fit_buchanan)["N_0"], t_lag = coef(fit_buchanan)["t_lag"])

gompertz_points <- gompertz_model(t = timepoints, r_max = coef(fit_gompertz)["r_max"], N_max = coef(fit_gompertz)["N_max"], N_0 = coef(fit_gompertz)["N_0"], t_lag = coef(fit_gompertz)["t_lag"])

df1 <- data.frame(timepoints, logistic_points)
df1$model <- "Logistic"
names(df1) <- c("Time", "PopBio", "model")

df2 <- data.frame(timepoints, baranyi_points)
df2$model <- "Baranyi"
names(df2) <- c("Time", "PopBio", "model")

# df3 <- data.frame(timepoints, buchanan_points)
# df3$model <- "Buchanan"
# names(df3) <- c("t", "LogN", "model")

df4 <- data.frame(timepoints, gompertz_points)
df4$model <- "Gompertz"
names(df4) <- c("Time", "PopBio", "model")

model_frame <- rbind(df2, df4)

ggplot(bae_cb, aes(x = Time, y = PopBio)) +
  geom_point(size = 1) +
  geom_line(data = model_frame, aes(x = Time, y = PopBio, col = model), size = 1) +
  theme_bw() + # make the background white
  theme(aspect.ratio=1)+ # make the plot square 
  labs(x = "Time", y = "PopBio")

