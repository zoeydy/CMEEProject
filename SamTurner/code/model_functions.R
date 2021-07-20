## Script: model_functions.R
## Author: Sam Turner sat19@ic.ac.uk
## About: Model specifications. 


# to fit in log space
Quadratic <- function(t, a, b,c){
  return(c*t**2 + b*t + a)
}


Cubic <- function(t, a, b, c, d){
  return(d*t**3 + c*t**2 + b*t + a)
}

Logistic <- function(t, N0, Nmax, r){
  N <- N0 * Nmax / (N0 + (Nmax - N0) * exp(-r*t))
  return( log10(N) )
}

Gompertz <- function(t, N_0, N_max, r_max, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (N_max - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((N_max - N_0) * log(10)) + 1)))
}

Baranyi <- function(t, N_0, N_max, r_max, t_lag){  # Baranyi model (Baranyi 1993)
  return(N_max + log10((-1+exp(r_max*t_lag) + exp(r_max*t))/(exp(r_max*t) - 1 + exp(r_max*t_lag) * 10^(N_max-N_0))))
}

Buchanan <- function(t, N_0, N_max, r_max, t_lag){ # Buchanan model - three phase logistic (Buchanan 1997)
  return(N_0 + (t >= t_lag) * (t <= (t_lag + (N_max - N_0) * log(10)/r_max)) * r_max * (t - t_lag)/log(10) + (t >= t_lag) * (t > (t_lag + (N_max - N_0) * log(10)/r_max)) * (N_max - N_0))
}

# to fit in linear space
Logistic.exp <- function(t, N0, Nmax, r){
  return( 10^(Logistic(t, N0, Nmax, r)  ))
}


Gompertz.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Gompertz(t, N0, Nmax, r, t_lag)  ))
}

Baranyi.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Baranyi(t, N0, Nmax, r, t_lag)  ))
}

Buchanan.exp <- function(t, N0, Nmax, r, t_lag){
  return( 10^(Buchanan(t, N0, Nmax, r, t_lag)  ))
}

Quadratic.exp <- function(t, a, b,c){
  return(c*t**2 + b*t + a)
}


Cubic.exp <- function(t, a, b, c, d){
  return(d*t**3 + c*t**2 + b*t + a)
}

# lists and vectors of model names and functions for fitting and plotting
all_models <- c("Logistic", "Gompertz", "Baranyi", "Buchanan", "Quadratic", "Cubic")
all_models.linear <- c("Logistic.exp", "Gompertz.exp", "Baranyi.exp", "Buchanan.exp", "Quadratic.exp", "Cubic.exp")

non_linear_models <- c("Logistic", "Gompertz", "Baranyi", "Buchanan")
linear_space_models <- c("Logistic.exp", "Gompertz.exp", "Baranyi.exp", "Buchanan.exp", "Quadratic.exp", "Cubic.exp")

draw_functions <- list("Logistic" = Logistic , "Gompertz" = Gompertz, "Baranyi" = Baranyi, "Buchanan" = Buchanan, "Quadratic" = Quadratic, "Cubic" = Cubic,"Logistic.exp" = Logistic.exp , "Gompertz.exp" = Gompertz.exp, "Baranyi.exp" = Baranyi.exp, "Buchanan.exp" = Buchanan.exp, "Quadratic.exp" = Quadratic.exp, "Cubic.exp" = Cubic.exp)
fit_functions <- list("Logistic" = Logistic , "Gompertz" = Gompertz, "Baranyi" = Baranyi, "Buchanan" = Buchanan, "Quadratic" = 2, "Cubic" = 3, "Logistic.exp" = Logistic.exp , "Gompertz.exp" = Gompertz.exp, "Baranyi.exp" = Baranyi.exp, "Buchanan.exp" = Buchanan.exp, "Quadratic.exp" = 2, "Cubic.exp" = 3)
model_params <- list("Logistic" = 3 , "Gompertz" = 4, "Baranyi" = 4, "Buchanan" = 4, "Quadratic" = 3, "Cubic" = 4, "Logistic.exp" = 3 , "Gompertz.exp" = 4, "Baranyi.exp" = 4, "Buchanan.exp" = 4, "Quadratic.exp" = 3, "Cubic.exp" = 4)
colours <- c("Logistic" = 'black' , "Gompertz" = 'red', "Baranyi" = 'dark green', "Buchanan" = 'blue', "Quadratic" = '#66a61e', "Cubic" = '#e6ab02', "Logistic.exp" = 'black' , "Gompertz.exp" = 'red', "Baranyi.exp" = 'dark green', "Buchanan.exp" = 'blue', "Quadratic.exp" = '#66a61e', "Cubic.exp" = '#e6ab02')
 