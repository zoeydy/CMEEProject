setwd("/home/nesbit/Desktop/Data_RVK/Code")
#setwd("../../Data_RVK")

rm(list=ls()) #clears workspace
df_final <- as.data.frame(matrix(0, nrow = 1, ncol = 3)) ## All figs included
names(df_final) = c("Species", "X", "Y")

for (j in 1:12){
## Load data
df <- read.csv(paste("../Data/All_Figs/All_Figs_Nesbit_", j, ".csv", sep = ""), header = FALSE, stringsAsFactors = F)
# df <- df[1:(dim(df)[1]-2), ] #removes the extra two points

## Iterate over all columns
n_data <- ncol(df)/2
df_each_fig <- as.data.frame(matrix(0, nrow = 1, ncol = 3))
names(df_each_fig) = c("Species", "X", "Y")

  for(i in 1:n_data) {
    df_temp <- df[ , ((i*2)-1) : (i*2)]
    ##NEED TO REMOVE EMPTY ROWS
    sp <- df_temp[1,1]
    df_temp <- df_temp[-c(1,2), ]
    df_temp <- cbind(rep(sp, nrow(df_temp)), df_temp)
    names(df_temp) = c("Species", "X", "Y")
    df_final <- rbind(df_final, df_temp)
    df_each_fig <- rbind(df_each_fig, df_temp)
  }
  
  write.csv(df_each_fig, paste("../Results/All_Figs/FigUpdatedNesbit_", j, ".csv", sep = ""))

} 

  write.csv(df_final, "../Results/All_Figs/All_Figs_Nesbit.csv")


# ## Converting one data set into correct format
# new <- df[ ,c(1,2)]
# sps <- new[1,1] #finds the species name
# new <- new[-c(1,2), ] #remove column names
# new <- cbind(Species = rep(sps, dim(new)[1]), new) #adding a column to the data frame with the species name
