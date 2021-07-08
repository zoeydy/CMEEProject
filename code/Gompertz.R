

rm(list = ls())
graphics.off()


require(minpack.lm)


# 1. define model and starting_value_guess functions
gompertz_model <- function(t, r_max, K, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
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

# 2. read the Data
Data <- read.csv('../data/pop.csv')
Data <- Data[order(Data[,'ID'], Data[,'Time']),]


# 3. writing for loop to get the best value for each ID_data by comparing AIC
start_value_df <- data.frame(ID=NULL, Model=NULL,AIC=NULL,RSq=NULL,r_max=NULL, N_0=NULL, K=NULL, t_lag = NULL)
plot_df <- data.frame()
comp <- list()
for (i in 1:length(unique(Data$ID))){
  #browser()
  # subset the Data by ID 
  id <- unique(Data$ID)[i]
  data <- subset(Data, Data$ID == id)
  
  # get starting values
  N0start <- min(data$logN)
  Kstart <- max(data$logN)
  tlagStart <- data$Time[which.max( diff( diff(data$logN)/diff(data$Time) )/(diff(data$Time[1:nrow(data)-1])) )]
  interval <- max(data$logN) - min(data$logN)
  pop_low <- min(data$logN) + 0.2*interval
  pop_high <- min(data$logN) + 0.75*interval
  if (length(data$Time) <= 5) {
    lm_fit <- lm(logN~Time, data)
  }else{
    if (nrow(data[data$logN > pop_low & data$logN < pop_high,])<2) {
      lm_fit <- lm(logN~Time, data[data$logN > min(data$logN) + 0.05*interval &
                                     data$logN < min(data$logN) + 0.95*interval,])
    }else{
      lm_fit <- lm(logN~Time, data[data$logN > pop_low & data$logN < pop_high,])
    }
  }
  rStart <- coef(lm_fit)[2]
    
  # sample
  SampTime <- 1000
  fact <- 2
  
  SampMin <- runif(SampTime,N0start-abs(N0start)*fact,N0start+abs(N0start)*fact)
  SampMax <- runif(SampTime,Kstart-abs(Kstart)*fact,Kstart+abs(Kstart)*fact)
  Samp_r_max <- runif(SampTime, 0 , rStart + abs(rStart)*fact)
  Samp_t_lag <- runif(SampTime, tlagStart-abs(tlagStart)*fact, tlagStart+abs(tlagStart)*fact)
  
  # for loop in each ID for finding starting value 
  start_id_list <- list()
  plot_id_list <- list()
  for (j in 1:SampTime){
    result <- GuessStart(data,t = Time, r_maxS = Samp_r_max[j], N_0S = SampMin[j], KS = SampMax[j],t_lagS = Samp_t_lag[j])
    start_id_list[[j]] <- result[1:8]
    plot_id_list[[j]] <- result[9]
  }
  
  # write the start_id_list into data frame
  n.obs <- sapply(start_id_list, length)
  seq.max <- seq_len(max(n.obs))
  start_id_mat <- t(sapply(start_id_list, "[", i = seq.max))
  start_id_df <- as.data.frame(start_id_mat)
  start_id_df$index <- seq(1:1000)
  colnames(start_id_df) <- c("AIC","AICc","BIC","RSq","r_max","N_0","K","t_lag","index")
  index_AIC <- start_id_df$index[which.min(start_id_df$AIC)]
  index_AICc <- start_id_df$index[which.min(start_id_df$AICc)]
  index_BIC <- start_id_df$index[which.min(start_id_df$BIC)]
  comp[[i]] <- all(sapply(list(index_AIC,index_BIC,index_AICc), function(x) x == index_AICc))
  
  # write the plot_id_list into data frame
  plot_id_df <- as.data.frame(plot_id_list[index_AICc])
  
  if (nrow(start_id_df[index_AICc,]) == 0) {
    start_id_df <- data.frame(ID=id,Model="gompertz",AIC=NA,AICc=NA,BIC=NA,RSq=NA,r_max=NA,N_0=NA,K=NA,t_lag=NA,index=NA)
    start_value_df <- rbind(start_value_df, start_id_df)

    plot_id_df <- data.frame(pred_gom=NA,ID=id,model="gompertz")
    plot_df <- rbind(plot_df, plot_id_df)
  }else{
    start_id_df <- data.frame(ID=id,model="gompertz",start_id_df[index_AICc,])
    colnames(start_id_df) <- c("ID","Model","AIC","AICc","BIC","RSq","r_max","N_0","K","t_lag","index")
    start_value_df <- rbind(start_value_df, start_id_df)
    
    plot_id_df$ID <- id
    plot_id_df$model <- "gompertz"
    colnames(plot_id_df) <- c("pred_gom", "ID", "model")
    plot_df <- rbind(plot_df, plot_id_df)
  }
}

# 4. save the starting value for each ID into csv
start_value_df <- apply(start_value_df,2,as.character)
write.csv(start_value_df, "../results/gompertz_Starting_Value.csv")
write.csv(plot_df, "../results/gompertz_plot_points.csv")
