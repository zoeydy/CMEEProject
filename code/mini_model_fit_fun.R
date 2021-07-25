
# model functions
gompertz_model <- function(t, r_max, K, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((K - N_0) * log(10)) + 1)))
} 
baranyi_model <- function(t, r_max, K, N_0, t_lag){  # Baranyi model (Baranyi 1993)
  return(K + log10((-1+exp(r_max*t_lag) + exp(r_max*t))/(exp(r_max*t) - 1 + exp(r_max*t_lag) * 10^(K-N_0))))
}
fun_models <- c("gompertz_model", "baranyi_model")
model_params <- list("Gompertz" = 4, "Baranyi" = 4)


# function to sample inits for one id
sample_inits <- function(data, samp.time, samp.fact){
  # get inits
  N0start <- min(data$logN)
  Kstart <- max(data$logN)
  tlagStart <- data$Time[which.max( diff( diff(data$logN)/diff(data$Time) )/(diff(data$Time[1:nrow(data)-1])) )]
  
  interval <- max(data$logN) - min(data$logN)
  pop_low <- min(data$logN) + 0.2*interval
  pop_high <- min(data$logN) + 0.75*interval
  midle_df <- data[data$logN > pop_low & data$logN < pop_high,]
  
  if (nrow(midle_df)<2) {
    inits.id <- data.frame(N0 = NA, K = NA, tlag = NA, rmax = NA)
  } else {
    lm_fit <- lm(logN~Time, midle_df)
    rStart <- coef(lm_fit)[2]
    
    # sample
    SampTime <- samp.time
    fact <- samp.fact
    SampMin <- runif(SampTime,N0start-abs(N0start)*fact,N0start+abs(N0start)*fact)
    SampMax <- runif(SampTime,Kstart-abs(Kstart)*fact,Kstart+abs(Kstart)*fact)
    Samp_t_lag <- runif(SampTime, tlagStart-abs(tlagStart)*fact, tlagStart+abs(tlagStart)*fact)
    Samp_r_max <- runif(SampTime, 0 , rStart + abs(rStart)*fact)
    
    inits.id <- data.frame(N0 = SampMin, K = SampMax, tlag = Samp_t_lag, rmax = Samp_r_max)
  }
  return(inits.id)
}
#test <- sample_inits(data = data, samp.time = 100,samp.fact = 1)

# function to fit and find the best starting value
GuessStart <- function(data, t, N_0, K, t_lag, r_max, fmla){
  #browser()
  out <- tryCatch(
    epr <- {
      nlsLM(fmla, data, start = list(N_0 = N_0, K = K, t_lag = t_lag,r_max = r_max))
      # nlsLM(fmla, data, start = list(N_0 = samp.inits[10,1], K = samp.inits[10,2], 
      #                                t_lag = samp.inits[10,3], r_max = samp.inits[10,4]))
      # nlsLM(fmla, data, start = list(N_0 = N_0, K = K, 
      #                                t_lag = t_lag, r_max = r_max))
      
      # logN~baranyi_model(t=Time, r_max, N_0, K, t_lag)
    },
    # warning = function(w){
    #   #print(w)
    #   return(NA)
    # },
    error = function(e){ 
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
    predi_bar <- predict(epr, newdata = list(Time = time2plot), interval = "confidence")
    #print(AIC)
    return(list(N_0, K, t_lag, r_max, AIC, AICc, BIC, rsq, predi_bar))
  }
  #print(out)
  return(out)
}
# test.guess <- GuessStart(data, t = Time, N_0 = samp.inits[10,1], K = samp.inits[10,2], 
#           t_lag = samp.inits[10,3], r_max = samp.inits[10,4],fmla)


# function to fit one model to one ID and return successful results
fit_id_model <- function(id, model,fmla){
  fit_result <- list()
  samp.inits <- sample_inits(data, samp.time, samp.fact)
  if (is.na(samp.inits[1,1])) { # check if there is not any starting value
    fit_result <- NULL
  } else{
    for (i in 1:samp.time) {
      # fmla <- as.formula(paste0("logN ~", model.name, "(t=Time, r_max, N_0, K, t_lag)"))
      result <- GuessStart(data,t = Time, N_0 = samp.inits[i,1], K = samp.inits[i,2], 
                           t_lag = samp.inits[i,3], r_max = samp.inits[i,4], fmla)
      if (!is.na(result[1])) {
        fit_result[[i]] <- result
      }
      fit_result <- plyr::compact(fit_result)[1]
    }
  }
  return(fit_result)
}
# test.id.model <- fit_id_model(idname, fun_models[1],fmla)






# # function to fit models to one ID 
# rStars <- function(id, models){
#   fit.id <- list()
#   for (i in 1:length(models)) {
#     model.name <- models[i]
#     fit.id.model <- fit_id_model(idname,model.name)
#     fit.id.model <- plyr::compact(fit.id.model)[[1]]
#     fit.id[[model.name]] <- fit.id.model
#   }
#   return(fit.id)
# }
# # test.id.models <- fit_id_models(idname,fun_models)

