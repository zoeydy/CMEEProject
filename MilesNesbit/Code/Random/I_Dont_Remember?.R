rm(list = ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
install.packages("minpack.lm")
require("minpack.lm")
for (j in 1:12){
  df <- read.csv(paste("../Data/All_Figs/All_Figs_Nesbit_", j, ".csv", sep = ""), header = FALSE, stringsAsFactors = F)
}
names(df) = c("Number", "Species", "X", "Y")
str(df)
df<-df[-1,]
XY <- c(3, 4)
plot(df$X, df$Y, main= "All Species All Temperatures Plotted", ylab= "Mystery Number (Biomass?)",
     xlab= "Even More Mysterious Number (Time?)") + theme_bw()
RVK<- function(r, N, K){
  return(rN(1-(N/K)))
}
PowFit <- nlsLM(Y ~ powMod(X, a, b), data = df, start = list(a = .1, b = .1))
Lengths <- seq(min(0), max(350),len=200)
plot(df$X, df$Y, main= "All Species All Temperatures Plotted", ylab= "Mystery Number (Biomass?)",
     xlab= "Even More Mysterious Number (Time?)") lines(Lengths)