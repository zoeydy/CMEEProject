## Script: demo_plots.R
## Author: Sam Turner sat19@ic.ac.uk
## About: makes plots to illustrate features of models and model fits. 


# clear environment

rm(list = ls())

# load packages

required_packages <- c("minpack.lm",
                       "ggplot2",
                       "dplyr",
                       "gridExtra"
                       )

for (pkg in required_packages){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

source("model_functions.R")

# demonstrate model parameters

time_range <- seq(0,30,len=600)

df <- data.frame(Time = time_range, Pop = Gompertz(time_range,1,20,5,5))
g <- ggplot(data = df, aes(x = Time, y = Pop)) + geom_line() + labs( x = "Time", y = expression(log[10](N)))
g <- g + geom_hline(yintercept = 1) +  geom_hline(yintercept = 20) + geom_abline( intercept =-10, slope=5/log(10), lty = 2 ) 
g <- g + geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), data = data.frame(x0= 16*log(10)/5, y0 = 6, x1=20*log(10)/5, y1 = 6), lty = 2, color = "red")
g <- g + geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), data = data.frame(x0= 20*log(10)/5, y0 = 6, x1=20*log(10)/5, y1 = 10), lty = 2, color = "red")
g <- g + annotate("text", x = 10.5, y = 7.5, label = "mu[max]", parse = T) + annotate("text", x = 6.5, y = 1.6, label= "t[lag]", parse = T)
g <- g + theme_bw()
ggsave("../results/4_param_demo.pdf", plot = g, width = 5, height = 5)


g

# demonstrate rolling regression
set.seed(1) 
time_range <- seq(0,30,len=20)
pop <- Gompertz(time_range, 1, 20,5,10) + rnorm(20,0,0.5)

g <- ggplot(data = data.frame(Time = time_range, Pop = pop), aes(x = Time, y = Pop)) + geom_point(pch = 20, cex = 2) + labs(x = "Time", y = expression(log[10]~(N)))
maxM <- NA
maxR <- 0
for (i in 1:15){
  M <- lm(pop[i:(i+5)] ~ time_range[i:(i+5)])
  g <-  g + geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), data = data.frame(x0= time_range[i], y0 = summary(M)$coef[1,1] + summary(M)$coef[2,1]*time_range[i], x1=time_range[i+5], y1 = summary(M)$coef[1,1] + summary(M)$coef[2,1]*time_range[i+5]), lty = 2, color = "grey")
  if (summary(M)$coef[2,1] > maxR){
    maxM <- M
    maxR <- summary(M)$coef[2,1]
  }
}
g <-  g + geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), data = data.frame(x0= (min(pop) - summary(maxM)$coef[1,1]) / summary(maxM)$coef[2,1], y0 = min(pop), x1=(max(pop) - summary(maxM)$coef[1,1]) / summary(maxM)$coef[2,1] , y1 = max(pop)), lty = 1, color = "red")
g <- g + geom_hline(yintercept =min(pop), size = .3) + theme_bw()                    
ggsave("../results/rolling_reg.pdf",plot = g, width = 5, height = 5)



# demonstrate parameter bias
time_range <- seq(0,30,len=200)
pop <- Gompertz(time_range, 1, 20,5,10) + rnorm(200,0,0.5)
fitdata <- data.frame(Time = time_range, PopBio = pop)
M.1 <- nlsLM(PopBio~Gompertz(time_range, N0, Nmax, rmax, t_lag), start = c(N0 = 1, Nmax = 20, rmax = 5, t_lag = 10), data = fitdata)
M.2 <- nlsLM(PopBio~Buchanan(time_range, N0, Nmax, rmax, t_lag), start = c(N0 = 1, Nmax = 20, rmax = 5, t_lag = 10), data = fitdata)

c<-unlist(coef(M.1))
d<-unlist(coef(M.2))

g <- ggplot(data = data.frame(Time = time_range, Pop = pop), aes(x= Time, y = Pop)) + geom_point(pch = 20, cex = .2) + labs(x = "Time", y = expression(log[10]~(N)))
g <- g + geom_line(data = data.frame(Time = time_range, Pop = Buchanan(time_range, d[1], d[2], d[3], d[4])), aes(x=Time, y = Pop), color = "red", size = 1)
g <- g + geom_line(data = data.frame(Time = time_range, Pop = Gompertz(time_range, c[1], c[2], c[3], c[4])), aes(x=Time, y = Pop), color = "blue", size = 1)
g <- g + geom_abline(slope=5/log(10),intercept=-21, lty = 2, col = "blue") + theme_bw()

ggsave("../results/fit_difference.pdf", plot = g, width = 5, height  = 5)



# demonstrate poorly constrained parameters
Time <- c(1,2,4,10,17, 19, 21, 23)
pops <- c(1,1,1,6,11,11,11,11) + rnorm(8,0,.1)

g <- ggplot(data = data.frame(Time = Time, Pops = pops), aes(x= Time, y = Pops)) + geom_point() + labs(y=expression(log[10]~(N)))
g <- g + geom_line(data = data.frame(Time = seq(0,23,0.1), Pops = Baranyi(seq(0,23,0.1), 1,11,1*log(10),5) ),aes(x= Time, y = Pops) )
g <- g + geom_line(data = data.frame(Time = seq(0,23,0.1), Pops = Baranyi(seq(0,23,0.1), 1,11,2*log(10),7.5) ),aes(x= Time, y = Pops), color = "blue")
g <- g + geom_line(data = data.frame(Time = seq(0,23,0.1), Pops = Baranyi(seq(0,23,0.1), 1,11,4*log(10),8.75) ),aes(x= Time, y = Pops), color = "red") 
g <- g + theme_bw()
ggsave("../results/constraint_demo.pdf", plot = g, width = 5, height = 5)


################################
# Model parameter demonstrations
################################

# time range to plot model demonstrations
time_range <- seq(0,30,len=600)

# compare shapes of the 4 parameter models
test_plot.fours <- function(fun,Time,...){
  df <- data.frame(Time = Time)
  for (f_name in fun){
    f <- draw_functions[[f_name]]
    df[[f_name]] <- f(Time, ...)
  }
  p <- ggplot(data = df, aes(x = Time))
  cs <- c(colours[["Gompertz"]], colours[["Buchanan"]], colours[["Baranyi"]])
  p <- p + geom_line(aes( y = Buchanan,color = "Buchanan"))
  p <- p + geom_line(aes( y = Baranyi,color = "Baranyi"))
  p <- p + geom_line(aes( y = Gompertz,color = "Gompertz"))
  p <- p + scale_color_manual(values = cs, name = "Model") + theme_bw() + labs(y = "log(N)")
  print(p)
}

# plot a model with specified parameters
test_plot <- function(fun,Time,...){
  Nt <- fun(Time, ...)
  df <- data.frame(Time, Nt)
  p <- ggplot(df, aes(x = Time, y = Nt)) + geom_line() + theme_bw() + labs(y = "log(N)")
  print(p)
}

# plot 4 parameter models with varying rmax
test_plot.range.rmax <- function(fun,Time,N0, Nmax, r,tlag){
  fun <- draw_functions[[fun]]
  r1 <- fun(Time, N0, Nmax, r[1], tlag)
  r2 <- fun(Time, N0, Nmax, r[2], tlag)
  r3 <- fun(Time, N0, Nmax, r[3], tlag)
  
  df <- data.frame(Time, r1=r1,r2=r2,r3=r3)
  cs <- c(colours[["Gompertz"]], colours[["Buchanan"]], colours[["Baranyi"]])
  
  p <- ggplot(data = df, aes(x = Time))
  p <- p + geom_line(aes( y = r1,color = "4"))
  p <- p + geom_line(aes( y = r2,color = "6"))
  p <- p + geom_line(aes( y = r3,color = "8"))
  p <- p + scale_color_manual(values = cs, name = expression(mu[max] )) + theme_bw() + labs(y = "log(N)", title = "A")
  return(p)
}

# plot 4 parameter modesl with varying tlag
test_plot.range.tlag <- function(fun,Time,N0, Nmax, r,tlag){
  fun <- draw_functions[[fun]]
  r1 <- fun(Time, N0, Nmax, r, tlag[1])
  r2 <- fun(Time, N0, Nmax, r, tlag[2])
  r3 <- fun(Time, N0, Nmax, r, tlag[3])
  
  df <- data.frame(Time, r1=r1,r2=r2,r3=r3)
  cs <- c(colours[["Gompertz"]], colours[["Buchanan"]], colours[["Baranyi"]])
  
  p <- ggplot(data = df, aes(x = Time))
  p <- p + geom_line(aes( y = r1,color = "2"))
  p <- p + geom_line(aes( y = r2,color = "3"))
  p <- p + geom_line(aes( y = r3,color = "4"))
  p <- p + scale_color_manual(values = cs, name = expression(t[lag])) + theme_bw() + labs(y = "log(N)", title = "B")
  return(p)
}

# make figure demonstrating model shapes
g1 <- test_plot(Logistic, time_range, 0.1,10, 2) + labs(title = "A")
g2 <- test_plot.fours(c("Gompertz", "Buchanan", "Baranyi"), time_range, 0.1,10, 2, 5) + labs(title = "B")
g.total <- grid.arrange(g1,g2,widths = c(5,6),nrow = 1)
ggsave("../results/model_demo.pdf", plot = g.total, width = 10, height = 4)


# make figure demonstrating effect of parameters
g1 <- test_plot.range.rmax("Baranyi", time_range, 0.1,10, c(1,2,3), 3)
g2 <- test_plot.range.tlag("Baranyi", time_range, 0.1,10, 2, c(4,8,12))
g.total <- grid.arrange(g1,g2,nrow = 1)
ggsave("../results/tlag_rmax_demo.pdf", plot = g.total, width = 8, height = 4)

