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
files <- read.csv('../Results/Formated/Total_processedLog.csv')


#create dataframe
df <- data.frame(matrix(ncol = 15, nrow = 0))
df_names <- c('Temp','Time_units', 'PopBio_unts', 'Species', 'Medium', 'Rep', 'Citation',
              'parameter','model', 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
colnames(df) <- df_names


#kill popbio
files <- select(files,-c(PopBio))

#discrete curves
qf <- files %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBioLog"))


#specify the model functions:
gompertz_model <- function(t, r_max, N_max, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}

baranyi_model <- function(t, r_max, N_max, N_0, t_lag){  # Baranyi model (Baranyi 1993)
  return(N_max + log10((-1+exp(r_max*t_lag) + exp(r_max*t))/(exp(r_max*t) - 1 + exp(r_max*t_lag) * 10^(N_max-N_0))))
}

buchanan_model <- function(t, r_max, N_max, N_0, t_lag){ # Buchanan model - three phase logistic (Buchanan 1997)
  return(N_0 + (t >= t_lag) * (t <= (t_lag + (N_max - N_0) * log(10)/r_max)) * r_max
         * (t - t_lag)/log(10) + (t >= t_lag) * (t > (t_lag + (N_max - N_0) * log(10)/r_max)) * (N_max - N_0))
}

#different angle the loop
for(i in 1:876) {  
  tryCatch({
    df_subset <- subset(qf , qf$data == unique(qf$data)[i])$data[[1]]
    
fit_gompertz_multi <- nls_multstart(PopBioLog ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag),
                                    data = df_subset,
                                    start_lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 0),
                                    start_upper = c(t_lag=20, r_max=10, N_0 = 6, N_max = 10),
                                    lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 1),
                                    iter = 500,
                                    supp_errors = "Y")

fit_baranyi_multi <- nls_multstart(PopBioLog ~ baranyi_model(t = Time, r_max, N_max, N_0, t_lag),
                                   data = df_subset,
                                   start_lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 0),
                                   start_upper = c(t_lag=20, r_max=10, N_0 = 6, N_max = 10),
                                   lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 1),
                                   iter = 500,
                                   supp_errors = "Y")

fit_buchanan_multi <- nls_multstart(PopBioLog ~ buchanan_model(t = Time, r_max, N_max, N_0, t_lag),
                                    data = df_subset,
                                    start_lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 0),
                                    start_upper = c(t_lag=20, r_max=10, N_0 = 6, N_max = 10),
                                    lower = c(t_lag=0, r_max=0, N_0 = 0, N_max = 1),
                                    iter = 500,
                                    supp_errors = "Y")

if(!is.null(fit_baranyi_multi)){
  aic_baranyi <- AIC(fit_baranyi_multi)
} else(
  aic_baranyi<- Inf
)


if(!is.null(fit_buchanan_multi)){
  aic_buchanan <- AIC(fit_buchanan_multi)
} else(
  aic_buchanan <- Inf
)

if(!is.null(fit_gompertz_multi)){
  aic_gompertz <- AIC(fit_gompertz_multi)
} else(
  aic_gompertz <- Inf
)


check <- c(aic_baranyi, aic_buchanan, aic_gompertz)

check2 <- which.min(check)

fits <- list(fit_baranyi_multi, fit_buchanan_multi, fit_gompertz_multi)


  

a <- as.data.frame(summary(fits[[check2]])$coeff)
a$parameter <- rownames(a)
a$model <- c('baranyi', 'buchanan', 'gompertz')[check2]
a$discrete_curve <- as.numeric(i)




df <- rbind(df,a)


if (df == '') stop("NA")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

#merge df together
qf$discrete_curve <-  seq.int(nrow(qf))

total <- merge(df,qf,by="discrete_curve")


#save
write.csv(total,"../Results/Formated/Total_processed_Log_RVK_3_models.csv", row.names = FALSE)

