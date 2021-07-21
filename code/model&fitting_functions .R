

# load package
require(minpack.lm)

# to fit in log space
Gompertz <- function(t, N_0, N_max, t_lag, r_max){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}
Baranyi <- function(t, N_0, N_max, t_lag, r_max){  # Baranyi model (Baranyi 1993)
  return(N_max + log10((-1+exp(r_max*t_lag) + exp(r_max*t))/(exp(r_max*t) - 1 + exp(r_max*t_lag) * 10^(N_max-N_0))))
}
Buchanan <- function(t, N_0, N_max, t_lag, r_max){ # Buchanan model - three phase logistic (Buchanan 1997)
  return(N_0 + (t >= t_lag) * (t <= (t_lag + (N_max - N_0) * log(10)/r_max)) * r_max * (t - t_lag)/log(10) + (t >= t_lag) * (t > (t_lag + (N_max - N_0) * log(10)/r_max)) * (N_max - N_0))
}

# to fit in linear space
Gompertz.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Gompertz(t, N0, Nmax, r, t_lag)  ))
}

Baranyi.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Baranyi(t, N0, Nmax, r, t_lag)  ))
}

Buchanan.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Buchanan(t, N0, Nmax, r, t_lag)  ))
}

# lists and vectors of model names and functions for fitting and plotting
all_models <- c("Gompertz", "Baranyi", "Buchanan")
all_models.linear <- c("Gompertz.exp", "Baranyi.exp", "Buchanan.exp")

non_linear_models <- c("Logistic", "Gompertz", "Baranyi", "Buchanan")
linear_space_models <- c("Logistic.exp", "Gompertz.exp", "Baranyi.exp", "Buchanan.exp", "Quadratic.exp", "Cubic.exp")

draw_functions <- list("Gompertz" = Gompertz, "Baranyi" = Baranyi, "Buchanan" = Buchanan)
fit_functions <- list("Gompertz" = Gompertz, "Baranyi" = Baranyi, "Buchanan" = Buchanan, "Gompertz.exp" = Gompertz.exp, "Baranyi.exp" = Baranyi.exp, "Buchanan.exp" = Buchanan.exp)
model_params <- list("Gompertz" = 4, "Baranyi" = 4, "Buchanan" = 4, "Gompertz.exp" = 4, "Baranyi.exp" = 4, "Buchanan.exp" = 4)
colours <- c("Gompertz" = 'red', "Baranyi" = 'dark green', "Buchanan" = "blue", "Gompertz.exp" = 'red', "Baranyi.exp" = 'dark green', "Buchanan.exp" = 'blue')



# function to retrieve intial values for a particular id
get_inits <- function(id){
  
  data.id<-dplyr::filter(inits, id == id)
  initvals=list()
  
  for (model in c("log", "lin")){
    initials.both <- dplyr::filter(data.id, loglin == model)$value
    sublist <- setNames(as.list(initials.both) , paste0("x", 1:length(initials.both)))
    initvals[[model]]<-setNames( as.list(initials.both) , paste0("x", 1:length(initials.both)) )
  }
  return(initvals)
}

# plot the models specified by initial parameters to check quality of initial parameters
check_inits <- function(id, models){
  
  id_inits <- get_inits(id)
  
  # filter to get desired time series
  data.id <- dplyr::filter(data, id == id)
  
  # get time range of run
  t_range <- seq(min(data.id$Time),max(data.id$Time),len=200)
  fit_df <- data.frame(t_range)
  for (model_name in models){
    if (model_name != "Logistic"){
      estimates  <- do.call(draw_functions[[model_name]], c(list(t_range), unname(id_inits$log)))
      fit_df[[paste0(model_name)]] <- estimates
    }
    if (model_name == "Logistic"){
      estimates  <- do.call(draw_functions[[model_name]], c(list(t_range), unname(id_inits$lin)[1:3]))
      fit_df[[paste0(model_name)]] <- estimates
    }
  }
  # get information about run
  species <- as.character(data.id$Species[1])
  temp <- as.character(data.id$Temp[1])
  medium <- as.character(data.id$Medium[1])
  repli <- as.character(data.id$Rep[1])
  
  # plot
  p <- ggplot() + geom_point(data = data.id, aes(x = Time, logPopBio) )
  y_m <- melt(fit_df, id.vars=c("t_range"))
  p <- p + geom_line(data = y_m, aes(x = t_range, y = value, group = variable, colour = variable))
  p <- p  + labs(color="Model") + ylab(paste("Log Population Biomass / ", data.id$PopBio_units)) + xlab("Time / Hours")  + theme_bw()
  print(p)
  
}  



# fit specififed models for a specified id
fit_models_multi <- function(id, models){
  # subset data
  subset <- data[data$id == id,]
  # list to store fitted models
  fit_objects <- list(names = models)
  
  for (model_name in models){
    # check if model should be fit in linear space
    linear <- (model_name %in% linear_space_models)
    
    model_fit_function <- fit_functions[[model_name]]
    model_draw_function <- draw_functions[[model_name]]
    
    # check if model can be fit with OLS
    if (!(model_name %in% c(non_linear_models, paste0(non_linear_models, ".exp")))){
      degree <- model_fit_function
      if (linear){
        result <- try(ModelFit <- lm(PopBio ~ poly(Time,degree = degree, raw=T), data=subset), silent = F)
      }
      else {
        result <- try(ModelFit <- lm(logPopBio ~ poly(Time,degree = degree, raw=T), data=subset), silent = F)
      }
    }
    
    # otherwise fit with nlsLM
    else if (model_name %in% c(non_linear_models, linear_space_models)){
      
      # generate list of parameter names of appropriate length
      id_inits <- get_inits(id)
      # if model is logistic, need linear space initial values
      if (model_name == "Logistic" | model_name == "Logistic.exp"){
        # set multistart iterations
        iters <- c(3,3,3)
        # produce param name vector
        params <- paste0("x", 1:3)
        adj_subset <- subset
        # set initial value ranges
        model_inits <- id_inits[["lin"]][-4]
        model_inits.lower <- c(x1 = model_inits[[1]]-0.5*abs(model_inits[[1]]), x2 = model_inits[[2]]-0.5*abs(model_inits[[2]]), x3 = model_inits[[3]]-0.5*abs(model_inits[[3]]))
        model_inits.upper <- c(x1 = model_inits[[1]]+0.5*abs(model_inits[[1]]), x2 = model_inits[[2]]+0.5*abs(model_inits[[2]]), x3 = model_inits[[3]]+0.5*abs(model_inits[[3]]))
        # set bounds
        pop <- subset$PopBio
        lowers = c(-Inf,mean(pop), -Inf)
        uppers = c(mean(pop),max(pop)*1.1, Inf)
      }
      else{
        # param iterations
        iters <- c(3,3,5,5)
        pop <- subset$logPopBio
        tmax <- max(subset$Time)
        # param name vector
        params <- paste0("x", 1:4)
        time_removed <- removed_times[id,"time_removed"]
        # adjust so first data point at t = 0 : found to improve model fit
        adj_subset <- subset
        adj_subset$Time <- adj_subset$Time - time_removed
        # set parameter ranges
        model_inits <- id_inits[["log"]]
        model_inits.lower <- c(x1 = model_inits[[1]]-0.6*abs(model_inits[[1]]), x2 = model_inits[[2]]-0.6*abs(model_inits[[2]]), x3 = model_inits[[3]]-0.8*abs(model_inits[[3]]), x4 = max(model_inits[[4]]-time_removed-0.1*tmax, -time_removed))
        model_inits.upper <- c(x1 = model_inits[[1]]+0.2*abs(model_inits[[1]]), x2 = model_inits[[2]]+0.6*abs(model_inits[[2]]), x3 = model_inits[[3]]+0.8*abs(model_inits[[3]]), x4 = model_inits[[4]]-time_removed+0.1*tmax)
        # set bounds
        min_pop <- min(adj_subset$logPopBio)
        max_pop <- max(adj_subset$logPopBio)
        mean_pop <- mean(adj_subset$logPopBio)
        ran_pop <- max(adj_subset$logPopBio) - min(adj_subset$logPopBio)
        lowers = c(min_pop - ran_pop / 2 - 5,mean_pop,-Inf,-round(time_removed,digits = 10))
        uppers = c(mean_pop, max_pop + ran_pop / 2,Inf,Inf)
      }
      
      # if linear space model, fit to linear PopBio
      if (linear) {
        nlslm_formula <- as.formula(paste("PopBio ~ model_fit_function(Time, ", paste(params, collapse= ","), ")"))
      }
      # if log space model, fit to  logPopBio
      else{
        nlslm_formula <- as.formula(paste("logPopBio ~ model_fit_function(Time, ", paste(params, collapse= ","), ")"))
      }
      # try to fit model with multistart
      try(ModelFit <- nls_multistart(nlslm_formula, d = adj_subset, init.lowers = model_inits.lower, init.uppers = model_inits.upper,lowers = lowers, uppers = uppers ,repeats = iters))
    }
    # if unsuccessful, fit using randomly sampled initial values
    if ( is.null(ModelFit)){
      try(ModelFit <- deep_fit(nlslm_formula, data = subset, inits = unlist(model_inits), lowers = lowers, uppers = uppers) )
      if (is.null(ModelFit)){
        print(paste0("ID ", as.character(id), ": Fit unsuccessful for ", model_name))
      }
    }
    fit_objects[[model_name]] <- ModelFit
    
  }
  return(fit_objects)
}


# fit formula using nlsLM across a grid-search of specified parameter ranges, with optional bounds.
# Return NULL if no models fit
nls_multistart <- function(fmla, d, init.lowers, init.uppers,repeats, lowers=NULL, uppers=NULL){
  
  # get number of parameters in model
  n_param <- length(init.lowers)
  
  # list of parameter samples
  param_ranges <- list()
  for (i in 1:n_param){
    param_ranges[[i]] <- seq(init.lowers[i], init.uppers[i], length.out = repeats[i])
  }
  # array of all possible parameter combinations
  param_combinations <- expand.grid(param_ranges)
  colnames(param_combinations) <- paste0("x", 1:n_param)
  best_M <- NULL
  best_RSS <- Inf
  # attempt to fit with each parameter combination
  for (j in 1:dim(param_combinations)[1]){
    M <- NULL
    try(M <- nlsLM(fmla, data = d, start = param_combinations[j,],lower=lowers, upper=uppers ,  control = list(maxiter = 500)), silent = T)
    if (!is.null(M)){
      RSS <- mean(resid(M)^2)
      # update best model if best RSS is beaten
      if (RSS < best_RSS){
        best_RSS <- RSS
        best_M <- M
      }
    }
  }
  return(best_M)
}



# attempt to fit model to data by randomly varying initial starting values for 10,000 fit attempts
deep_fit <- function(fmla, data, inits, lowers = NULL, uppers = NULL){
  set.seed(4)
  # maximum of 1000 attempts to fit model
  for (k in 1:1000){
    # set initial values to try
    this_inits <- inits*(runif(min = 1,max=3, n = length(inits))^sample(c(-1,1), 4, replace = T))
    this_inits[[4]] <- this_inits[[4]] + runif(1,-0.2,0.2)*max(data$Time)
    M <- NULL
    t <-try( M <- nlsLM(fmla, data = data , start = this_inits, lower = lowers, upper = uppers),silent = T)
    # return model if fit successfully
    if (class(t) != "try-error"){
      return(M)
    }
  }
  return(NULL)
}


######################
# Plotting functions #
######################

# fit and plot specified models for a specified id
plot_fit_multi <- function(id, models, lin = F) {
  subset <- data[data$ID == id,]
  
  # get ID info
  species <- as.character(subset$Species[1])
  temp <- as.character(subset$Temp[1])
  medium <- as.character(subset$Medium[1])
  
  t_range <- seq(min(subset$Time)*0.9, max(subset$Time)*1.1, length.out  = 100)
  
  fitted.models <- fit_models_multi(id, models)
  
  predicted.values <- data.frame(t_range = t_range)
  
  # for each model, if it was fit successfully:
  for (model in names(fitted.models)[-1]){
    if (is.list(fitted.models[[model]])){
      
      # if a removed time needs to be added to t_lag:
      if (model %in% c("Baranyi", "Buchanan", "Gompertz", "Baranyi.exp", "Buchanan.exp", "Gompertz.exp")){
        
        # get params
        params <- unname(coef(fitted.models[[model]]))
        params[[4]] <- params[[4]] + removed_times[id,"time_removed"]
        
      }
      
      else{
        params <- unname(coef(fitted.models[[model]]))
      }
      
      # get model predicitons
      N.pred <- do.call(draw_functions[[model]], c(list(t_range), params))
      
      # check if model prediction in linear space
      if (model %in% linear_space_models){
        predicted.values[[paste0(model)]] <- log10(N.pred)
      }
      else{
        predicted.values[[paste0(model)]] <- N.pred
      }
      
    } 
  }
  
  
  # if non linear plot specified
  if (lin == F){
    # convert to long form and plot
    predicted.values.long <- melt(predicted.values, id.vars=c("t_range"))
    p <- ggplot() + geom_point(data = subset, aes(x = Time, logPopBio) )
    p <- p + geom_line(data = predicted.values.long, aes(x = t_range, y = value, group = variable, colour = variable))
    p <- p + ylab(paste("log(N) / ", subset$PopBio_units)) + xlab("Time / Hours")+ theme_bw() + labs(color="Model")
    print(p)
  }
  else{
    predicted.values.long <- melt(predicted.values, id.vars=c("t_range"))
    # convert model and data to linear space
    subset$logPopBio <- 10^subset$logPopBio
    predicted.values.long$value <- 10^predicted.values.long$value
    p <- ggplot() + geom_point(data = subset, aes(x = Time, logPopBio) )
    p <- p + geom_line(data = predicted.values.long, aes(x = t_range, y = value, group = variable, colour = variable))
    p <- p + ggtitle(paste(c(species, temp, medium), collapse =", " ))  + ylab(paste("Log Population Biomass / ", subset$PopBio_units)) + xlab("Time / Hours") + theme_bw()# + labs(color="Model") + 
    print(p)
  }
  
  
}

# function to compare fitting in log space against fitting in linear space.
plot_fit_multi.compare <- function(id, models) {
  subset <- data[data$ID == id,]
  
  t_range <- seq(min(subset$Time)*0.9, max(subset$Time)*1.1, length.out  = 100)
  
  # log space fits
  fitted.models = fit_models_multi(id, models)
  predicted.values <- data.frame(t_range = t_range)
  
  for (model in names(fitted.models)[-1]){
    if (is.list(fitted.models[[model]])){
      if (model %in% c("Baranyi", "Buchanan", "Gompertz","Baranyi.exp", "Buchanan.exp", "Gompertz.exp")){
        params <- unname(coef(fitted.models[[model]]))
        params[[4]] <- params[[4]] + removed_times[id,"time_removed"]
      }
      else{
        params <- unname(coef(fitted.models[[model]]))
      }
      N.pred <- do.call(draw_functions[[model]], c(list(t_range), params))
      
      if (model %in% linear_space_models){
        predicted.values[[paste0(model)]] <- log10(N.pred)
      }
      else{
        predicted.values[[paste0(model)]] <- N.pred
      }
      
    } 
  }
  
  # fit in log, plot in log
  predicted.values.long <- melt(predicted.values, id.vars=c("t_range"))
  p1 <- ggplot() + geom_point(data = subset, aes(x = Time, logPopBio) )
  p1 <- p1 + geom_line(data = predicted.values.long, aes(x = t_range, y = value, group = variable, colour = variable))
  p1 <- p1 + ggtitle("A" )  + ylab(paste("Log Population Biomass / ", subset$PopBio_units)) + xlab("Time / Hours")+ theme_bw() + theme(legend.position = "none")# + labs(color="Model")
  
  
  # fit in log, plot in linear
  predicted.values.long.linear <- melt(predicted.values, id.vars=c("t_range"))
  subset$logPopBio <- 10^subset$logPopBio
  predicted.values.long.linear$value <- 10^predicted.values.long.linear$value
  p2 <- ggplot() + geom_point(data = subset, aes(x = Time, logPopBio) )
  p2 <- p2 + geom_line(data = predicted.values.long.linear, aes(x = t_range, y = value, group = variable, colour = variable))
  p2 <- p2 + ggtitle("B" )  + ylab(paste("Population Biomass / ", subset$PopBio_units)) + xlab("Time / Hours") + theme_bw() + theme(legend.position = "bottom") + labs(color="Model") 
  
  
  # fit in linear scale
  models.exp <- paste0(models, ".exp")
  fitted.models.exp = fit_models_multi(id, models.exp)
  predicted.values.exp <- data.frame(t_range = t_range)
  
  for (model in names(fitted.models.exp)[-1]){
    if (is.list(fitted.models.exp[[model]])){
      if (model %in% c("Baranyi", "Buchanan", "Gompertz","Baranyi.exp", "Buchanan.exp", "Gompertz.exp")){
        params <- unname(coef(fitted.models.exp[[model]]))
        params[[4]] <- params[[4]] + removed_times[id,"time_removed"]
      }
      else{
        params <- unname(coef(fitted.models.exp[[model]]))
      }
      N.pred <- do.call(draw_functions[[model]], c(list(t_range), params))
      
      if (model %in% linear_space_models){
        predicted.values.exp[[paste0(model)]] <- (N.pred)
      }
      else{
        predicted.values.exp[[paste0(model)]] <- N.pred
      }
      
    } 
  }
  
  # fit in linear, plot in linear
  predicted.values.long.exp <- melt(predicted.values.exp, id.vars=c("t_range"))
  p3 <- ggplot() + geom_point(data = subset, aes(x = Time, logPopBio) )
  p3 <- p3 + geom_line(data = predicted.values.long.exp, aes(x = t_range, y = value, group = variable, colour = variable))
  p3 <- p3 + ggtitle("C" )  + ylab(paste("Population Biomass / ", subset$PopBio_units)) + xlab("Time / Hours")+ theme_bw() + theme(legend.position = "none")# + labs(color="Model")
  
  mylegend <- make.legend(p2)
  
  
  total_plot <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                         p2 + theme(legend.position="none"),
                                         p3 + theme(legend.position="none"),
                                         nrow=1),
                             mylegend, nrow=2,heights=c(10, 1))
  
  return(total_plot)
}

# make legend for grid.arrange plots
make.legend<-function(my_plot){
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

###############################



# function to produce weighted parameter estimates from akaike weight df and parameter estimate df.
weighted_param_est <- function(params, A_weights){
  weighted_params <- list()
  for (id in rownames(params)){
    # count number of models to be averaged
    permitted <- "&"(!is.na(unlist(params[id,])), is.finite(unlist(params[id,])))
    params_id <- unlist(params[id,])[permitted]
    ws <- unlist(A_weights[id,])[permitted]
    # weighted model average
    weighted_params[[id]] <- sum(params_id * ws) / sum(ws)
  }
  return(weighted_params)
}

