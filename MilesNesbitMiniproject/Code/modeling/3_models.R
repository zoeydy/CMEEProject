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


#kill popbio
files <- select(files,-c(PopBio))

#discrete curves
qf <- files %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBioLog"))


#the loop
for(i in 1:876) {  
  tryCatch({
    df_subset <- subset(qf , qf$data == unique(qf$data)[i])$data[[1]]
    
    #start
    N_0_start <- min(df_subset$PopBioLog)
    N_max_start <- max(df_subset$PopBioLog)
    t_lag_start <- df_subset$Time[which.max(diff(diff(df_subset$PopBioLog)))] #check what exactly this is evaluating optimize lag value
    r_max_start <- max(diff(df_subset$PopBioLog))/mean(diff(df_subset$Time)) #ask Tom how he optimized
    

    fit_baranyi <- nlsLM(PopBioLog ~ baranyi_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                         list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_buchanan <- nlsLM(PopBioLog ~ buchanan_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_gompertz <- nlsLM(PopBioLog ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    #summary
    
 
    if(exists("fit_baranyi")) {  
      y <-as.data.frame(summary(fit_baranyi)$coeff)
      y$parameter <- rownames(y)
      y$model <- "baranyi"
      y$discrete_curve <- as.numeric(i)
    }else {
      y <- NA
      
    }
    
    
    
    if(exists("fit_buchanan"))  {
      z <-as.data.frame(summary(fit_buchanan)$coeff)
      z$parameter <- rownames(z)
      z$model <- "buchanan"
      z$discrete_curve <- as.numeric(i)
    }else {
      z <- NA
      
    }
    
    
    
    if(exists("fit_gompertz")){
      a <- as.data.frame(summary(fit_gompertz)$coeff)
      a$parameter <- rownames(a)
      a$model <- "gompertz"
      a$discrete_curve <- as.numeric(i)
    }else {
      a <- NA
      
    }
    
    
    df <- rbind(df,y,z,a)
    
    if (df == '') stop("NA")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


#tie together new df and old df
qf$discrete_curve <-  seq.int(nrow(qf))

total <- merge(df,qf,by="discrete_curve")

#clean
drops2 <- c('data')
total <- total[ , !(names(total) %in% drops2)]


#sanity check
length(unique(total$discrete_curve))


#subsetting absorbance
Absorbance <- total %>%
  filter(str_detect(PopBio_unts, 'Abs'))

ABSD <- total %>%
  filter(str_detect(PopBio_unts, 'OD'))

Absorbance <- rbind(Absorbance, ABSD)

#plot
ggplot(Absorbance, aes(x= log(r_max), y= log(N_max), color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Absorbance")
