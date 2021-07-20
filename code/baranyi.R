
rm(list = ls())
graphics.off()

# 1. read the Data
Data <- read.csv('../data/2_TomData.csv')
Data <- Data[order(Data[,'ID'], Data[,'Time']),]

# 2. define model and GuessStart functions
baranyi_model <- function(t, r_max, K, N_0, t_lag){
  return(N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((K - N_0) * log(10)) + 1)))
}   

GuessStart <- function(t, r_maxS, N_0S, KS, t_lagS, data){
  #browser()
  out <- tryCatch(
    epr <- {
      nlsLM(logN~gompertz_model(t=Time, r_max, N_0, K, t_lag), data,
            start = list(r_max = r_maxS, N_0 = N_0S, K = KS, t_lag = t_lagS))
    },
    # warning = function(w){
    #   #print(w)
    #   return(NA)
    # },
    error = function(e){ 
      #print(e)
      return(NA) 
    }
  )
  if (!is.na(out)){
    AIC <- AIC(out)
    n <- nrow(data)
    k <- 4
    AICc <- AIC(out) + (2*k^2+2*k)/(n-k-1)
    BIC <- BIC(out)
    # calculate the rsq
    RSS <- sum(residuals(out)^2)
    TSS <- sum((data$PopBio-mean(data$PopBio))^2)
    rsq <- 1-(RSS/TSS)
    # calculate the plot datapoints
    time2plot <- seq(min(data$Time), max(data$Time), length=1000)
    predi_gom <- predict(epr, newdata = list(Time = time2plot), interval = "confidence")
    #print(AIC)
    return(list(AIC, AICc, BIC, rsq, r_maxS, N_0S, KS, t_lagS, predi_gom))
  }
  #print(out)
  return(out)
}



