rm(list = ls())

# read the data
Data <- read.csv('../data/processed_data.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]

# get starting values
start_value_df <- data.frame(loglin=rep(c(rep('log',4),rep('lin',4)),length(unique(Data$id))), 
                             param=rep(c('N0','Nmax','rmax','tlag'),2*length(unique(Data$id))))
start_value <- data.frame(ID = NULL, value = NULL)
for (i in 1:length(unique(Data$id))) {
  idname <- unique(Data$id)[i]
  dat <- Data[Data$id == idname, ]
  # rmax log
  interval <- max(dat$logPopBio) - min(dat$logPopBio)
  pop_low <- min(dat$logPopBio) + 0.2*interval
  pop_high <- min(dat$logPopBio) + 0.75*interval
  if (length(dat$Time) <= 5) {
    log_lm_fit <- lm(logPopBio~Time, dat)
  }else{
    if (nrow(dat[dat$logPopBio > pop_low & dat$logPopBio < pop_high,])<2) {
      log_lm_fit <- lm(logPopBio~Time, dat[dat$logPopBio > min(dat$logPopBio) + 0.05*interval &
                                             dat$logPopBio < min(dat$logPopBio) + 0.95*interval,])
    }else{
      log_lm_fit <- lm(logPopBio~Time, dat[dat$logPopBio > pop_low & dat$logPopBio < pop_high,])
    }
  }
  # rmax lin
  interval <- max(dat$PopBio) - min(dat$PopBio)
  pop_low <- min(dat$PopBio) + 0.2*interval
  pop_high <- min(dat$PopBio) + 0.75*interval
  if (length(dat$Time) <= 5) {
    lin_lm_fit <- lm(PopBio~Time, dat)
  }else{
    if (nrow(dat[dat$PopBio > pop_low & dat$PopBio < pop_high,])<2) {
      lin_lm_fit <- lm(PopBio~Time, dat[dat$PopBio > min(dat$PopBio) + 0.05*interval &
                                          dat$PopBio < min(dat$PopBio) + 0.95*interval,])
    }else{
      lin_lm_fit <- lm(PopBio~Time, dat[dat$PopBio > pop_low & dat$PopBio < pop_high,])
    }
  }
  # write data
  lin.df <- data.frame(id = rep(idname,4),value = c(min(dat$PopBio),max(dat$PopBio),coef(lin_lm_fit)[2],
                                                    dat$Time[which.max(diff(diff(dat$PopBio)/diff(dat$Time))/(diff(dat$Time[1:nrow(dat)-1])))]))
  log.df <- data.frame(id = rep(idname,4), value = c(min(dat$PopBio),max(dat$PopBio),coef(log_lm_fit)[2],
                                                     dat$Time[which.max( diff( diff(dat$logPopBio)/diff(dat$Time) )/(diff(dat$Time[1:nrow(dat)-1])) )]))
  df <- rbind(log.df, lin.df)
  start_value <- rbind(start_value, df)
}
start_value_df <- cbind(start_value_df,start_value)

write.csv(start_value_df, "../data/starting_value.csv")