setwd("/home/nesbit/Desktop/Data_RVK/Code")
rm(list=ls()) #clears workspace
for (j in 1:12){
df <- read.csv(paste("../Data/All_Figs/All_Figs_Nesbit_", j, ".csv", sep = ""), header = FALSE, stringsAsFactors = F)
}

names(df) = c("Number", "Species", "X", "Y")
str(df)
df<-df[-1,]
XY <- c(3, 4)
plot(df$X, df$Y, main= "All Species All Temperatures Plotted", ylab= "Mystery Number (Biomass?)",
     xlab= "Even More Mysterious Number (Time?)") + theme_bw()

