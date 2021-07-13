#install.packages('nls.multstart')
#install.packages('taxize')
# install.packages('zoo')
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

#clean
# drops <- c("X", "X.1")
# files <- files[ , !(names(files) %in% drops)]

#create dataframe
df <- data.frame(matrix(ncol = 15, nrow = 0))
df_names <- c('Temp','Time_units', 'PopBio_unts', 'Species', 'Medium', 'Rep', 'Citation',
              'parameter','model', 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
colnames(df) <- df_names


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
    
    #fit models
    # fit_logistic <- nlsLM(PopBio ~ logistic_model(t = Time, r_max, N_max, N_0), df_subset,
    #                       list(r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_baranyi <- nlsLM(PopBioLog ~ baranyi_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                         list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_buchanan <- nlsLM(PopBioLog ~ buchanan_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_gompertz <- nlsLM(PopBioLog ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    #summary
      
# if(exists("fit_logistic")){
#   x <-as.data.frame(summary(fit_logistic)$coeff)
#   x$parameter <- rownames(x)
#   x$model <- "logistic"
#   x$discrete_curve <- as.numeric(i)
# }else {
#   x <- NA
#   
#     }
#     
   
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

#save
#write.csv(total,"../Results/analysis/initial_analysis.csv")
  
#sanity check
length(unique(total$discrete_curve))
#297/876 curves fit at all


#RVK
RVK <- total[!grepl("N_0", total$parameter), ]
RVK <- RVK[!grepl("t_lag", RVK$parameter), ]


#plot
RVK %>%
  ggplot(aes(x=parameter,y=Estimate))+
  geom_point()+
  facet_wrap(~paste(Species), scales = "free")


#alternatively with growth curver
#fit growth curves
af <- files %>%
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

#percentage_of_fitted_curves_using_growth_curver
percentage_of_fitted_curves_using_growth_curver <- 667/876

#taking values of r and k and putting it for starting values in loop
bf <- files %>%
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

# get just r
bf2 <- subset(bf, select = c('r', 'K'))


#ef
ef <- data.frame(matrix(ncol = 15, nrow = 0))
ef_names <- c('Temp','Time_units', 'PopBio_unts', 'Species', 'Medium', 'Rep', 'Citation',
              'parameter','model', 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
colnames(ef) <- ef_names


#loop pt 2
#the loop
for(i in 1:876) {  
  tryCatch({
    ef_subset <- subset(ef , ef$data == unique(ef$data)[i])$data[[1]]
    
    #start
    N_0_start <- min(df_subset$PopBio)
    N_max_start <- bf2[i,2]
    t_lag_start <- df_subset$Time[which.max(diff(diff(df_subset$PopBio)))] #check what exactly this is evaluating optimize lag value
    r_max_start <- bf2[i,1]
    
    #fit models
    fit_logistic <- nlsLM(PopBio ~ logistic_model(t = Time, r_max, N_max, N_0), df_subset,
                          list(r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_baranyi <- nlsLM(PopBio ~ baranyi_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                         list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_buchanan <- nlsLM(PopBio ~ buchanan_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    
    fit_gompertz <- nlsLM(PopBio ~ gompertz_model(t = Time, r_max, N_max, N_0, t_lag), df_subset,
                          list(t_lag=t_lag_start, r_max=r_max_start, N_0 = N_0_start, N_max = N_max_start))
    #summary
    
    if(exists("fit_logistic")){
      b <-as.data.frame(summary(fit_logistic)$coeff)
      b$parameter <- rownames(b)
      b$model <- "logistic"
      b$discrete_curve <- as.numeric(i)
    }else {files$Species[1]
      b <- NA
      
    }
    
    
    if(exists("fit_baranyi")) {  
      c <-as.data.frame(summary(fit_baranyi)$coeff)
      c$parameter <- rownames(c)
      c$model <- "baranyi"
      c$discrete_curve <- as.numeric(i)
    }else {
      c <- NA
      
    }
    
    
    
    if(exists("fit_buchanan"))  {
      d <-as.data.frame(summary(fit_buchanan)$coeff)
      d$parameter <- rownames(d)
      d$model <- "buchanan"
      d$discrete_curve <- as.numeric(i)
    }else {
      d <- NA
      
    }
    
    
    
    if(exists("fit_gompertz")){
      e <- as.data.frame(summary(fit_gompertz)$coeff)
      e$parameter <- rownames(e)
      e$model <- "gompertz"
      e$discrete_curve <- as.numeric(i)
    }else {
      e <- NA
      
    }
    
    
    ef <- rbind(ef,b,c,d,e)
    
    if (df == 'NA') stop("NA")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


#tie together new df and old df
qf$discrete_curve <-  seq.int(nrow(qf))

total2 <- merge(ef,qf,by="discrete_curve")

#clean
drops2 <- c('data')
total2 <- total2[ , !(names(total2) %in% drops2)]

#save
#write.csv(total,"../Results/analysis/initial_analysis.csv")

#sanity check
length(unique(total2$discrete_curve))


#different angle
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
  ggtitle("Absorbance")+
  facet_wrap(~Temp)

ggplot(ff, aes(x= (r), y= (K)))+
  geom_point(size = 4)

#discrete curves
wf <- Absorbance %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio"))
#length(unique(wf$discrete_curve))

#plot this one too
Absorbance %>%
  ggplot(aes(x=Time,y=PopBio,colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~paste(Species,Medium, Rep), scales = "free")

#subsetting cfu
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
  ggtitle("Colony Forming Units")+
  facet_wrap(~Temp)

plot(ff$r, ff$K)

#discrete curves
xf <- CFU %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,data = c("Time","PopBio"))


#plot this one too
CFU %>%
  ggplot(aes(x=Time,y=PopBio,colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~paste(Species,Medium, Rep), scales = "free")

length(unique(files$Species))


files$Species[1]

#taxonomic data
gnr_resolve("Chryseobacterium.balustinum")

classification("Chryseobacterium balustinum", db = "ncbi")

