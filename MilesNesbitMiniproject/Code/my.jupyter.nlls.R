setwd("/home/nesbit/Desktop/R/TheMulQuaBio-master/data")
install.packages("minpack.lm")
require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)




rm(list = ls())
graphics.off()


#write the function
powMod <- function(x, a, b) {
  return(a * x^b)
}


#samraats data
MyData <- read.csv("../data/GenomeSize.csv")

head(MyData)


#subset
Data2Fit <- subset(MyData,Suborder == "Anisoptera")

Data2Fit <- Data2Fit[!is.na(Data2Fit$TotalLength),] # remove NA's

#plot
plot(Data2Fit$TotalLength, Data2Fit$BodyWeight)

PowFit <- nlsLM(BodyWeight ~ powMod(TotalLength, a, b), data = Data2Fit, start = list(a = .1, b = .1))

summary(PowFit)

#visualize the fit
Lengths <- seq(min(Data2Fit$TotalLength),max(Data2Fit$TotalLength),len=200)
coef(PowFit)["a"]
coef(PowFit)["b"]
Predic2PlotPow <- powMod(Lengths,coef(PowFit)["a"],coef(PowFit)["b"])
plot(Data2Fit$TotalLength, Data2Fit$BodyWeight)
lines(Lengths, Predic2PlotPow, col = 'blue', lwd = 2.5)


#check
confint(PowFit)


#practicals


#abundances!!!!!!!!!!
rm(list = ls())
graphics.off()


#Let's first generate some "data" on the number of bacterial cells as a function of time that we can play with:
time <- c(0, 2, 4, 6, 8, 10, 12, 16, 20, 24) # timepoints, in hours
log_cells <- c(3.62, 3.62, 3.63, 4.14, 5.23, 6.27, 7.57, 8.38, 8.70, 8.69) # logged cell counts - more on this below

set.seed(1234) # set seed to ensure you always get the same random sequence if fluctuations  

data <- data.frame(time, log_cells + rnorm(length(time),sd=.1)) # add some random error

names(data) <- c("t", "LogN")

head(data)


#plot data
ggplot(data, aes(x = t, y = LogN)) + geom_point()


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


#generate some starting values for the NLLS fitting
N_0_start <- min(data$LogN)
N_max_start <- max(data$LogN)
t_lag_start <- data$t[which.max(diff(diff(data$LogN)))]
r_max_start <- max(diff(data$LogN))/mean(diff(data$t))


#fit the models
fit_logistic <- nlsLM(LogN ~ logistic_model(t = t, r_max, N_max, N_0), data,
                      list(r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_baranyi <- nlsLM(LogN ~ baranyi_model(t = t, r_max, N_max, N_0, t_lag), data,
                     list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_buchanan <- nlsLM(LogN ~ buchanan_model(t = t, r_max, N_max, N_0, t_lag), data,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

fit_gompertz <- nlsLM(LogN ~ gompertz_model(t = t, r_max, N_max, N_0, t_lag), data,
                      list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))

#model summaries
summary(fit_logistic)
summary(fit_baranyi)
summary(fit_buchanan)
summary(fit_gompertz)


#check fits
timepoints <- seq(0, 24, 0.1)

logistic_points <- logistic_model(t = timepoints, r_max = coef(fit_logistic)["r_max"], N_max = coef(fit_logistic)["N_max"], N_0 = coef(fit_logistic)["N_0"])

baranyi_points <- baranyi_model(t = timepoints, r_max = coef(fit_baranyi)["r_max"], N_max = coef(fit_baranyi)["N_max"], N_0 = coef(fit_baranyi)["N_0"], t_lag = coef(fit_baranyi)["t_lag"])

buchanan_points <- buchanan_model(t = timepoints, r_max = coef(fit_buchanan)["r_max"], N_max = coef(fit_buchanan)["N_max"], N_0 = coef(fit_buchanan)["N_0"], t_lag = coef(fit_buchanan)["t_lag"])

gompertz_points <- gompertz_model(t = timepoints, r_max = coef(fit_gompertz)["r_max"], N_max = coef(fit_gompertz)["N_max"], N_0 = coef(fit_gompertz)["N_0"], t_lag = coef(fit_gompertz)["t_lag"])

df1 <- data.frame(timepoints, logistic_points)
df1$model <- "Logistic"
names(df1) <- c("t", "LogN", "model")

df2 <- data.frame(timepoints, baranyi_points)
df2$model <- "Baranyi"
names(df2) <- c("t", "LogN", "model")

df3 <- data.frame(timepoints, buchanan_points)
df3$model <- "Buchanan"
names(df3) <- c("t", "LogN", "model")

df4 <- data.frame(timepoints, gompertz_points)
df4$model <- "Gompertz"
names(df4) <- c("t", "LogN", "model")

model_frame <- rbind(df1, df2, df3, df4)

ggplot(data, aes(x = t, y = LogN)) +
  geom_point(size = 3) +
  geom_line(data = model_frame, aes(x = t, y = LogN, col = model), size = 1) +
  theme_bw() + # make the background white
  theme(aspect.ratio=1)+ # make the plot square 
  labs(x = "Time", y = "log(Abundance)")

# (a) Calculate the confidence intervals on the parameters of each of the three fitted models, 
#and use model selection (using AIC and/or BIC) as you did before to see if you can determine the best-fitting model among the three.
# 
# (b) Alternatively, for a different random sequence of fluctuations, one or more of the models may fail to fit
#(a singular gradiant matrix error). Try repeating the above fitting with a different random seed (change the integers given to the random.seed( )
#function), or increase the sampling error by increasing the standard deviationand see if it happens. If/when the NLLS optimization
#does fail to converge (the RSS minimum was not found), you can try to fix it by chaning the starting values.
# 
# (c) Repeat the model comparison exercise 1000 times (You will have to write a loop), and determine if/whether one model generally
#wins more often than the others. Note that each run will generate a slightly different dataset, because we are adding a vector of
#random errors every time the "data" are generated. This may result in failure of the NLLS fitting to converge, in which case you
#will need to use the try() or tryCatch functions.
# 
# (d) Repeat (b), but increase the error by increasing the standard deviation of the normal error distributon,
#and see if there are differences in the robustness of the models to sampling/experimental errors. You may also want to
#try changing the distribution of the errors to some non-normal distribution and see what happens



