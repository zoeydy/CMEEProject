## Script: data_preparation.R
## Author: Sam Turner sat19@ic.ac.uk
## About: Loads and wrangles data set, rescaling time values and removing data from incorrect citation



# lear environment

rm(list=ls())

# load raw data

data <- read.csv("../data/LogisticGrowthData.csv") 
data$Time <- as.numeric(data$Time)

# remove incorrect time series
cit <- "Phillips, J.D. and Griffiths, M.W., 1987. The relation between temperature and growth of bacteria in dairy products. Food Microbiology, 4(2), pp.173-185."
data <- data[data$Citation != cit,]

# initialise experiment ID column
data$ID <- NA

# get all unique combinations of Species, Temperature, Medium, Replicate -> each unique combo is a run
combos <- unique(cbind(as.vector(data$Species), as.vector(data$Temp),as.vector(data$Medium),as.vector(data$Rep)))

# allocate IDs
for (i in 1:dim(combos)[1]){
  combo<-combos[i,]
  species <- combo[1]
  temp <- combo[2]
  medium <- combo[3]
  repli <- combo[4]
  
  data[(data$Species == species) & (data$Temp == temp) & (data$Medium == medium) & (data$Rep == repli),]$ID <- i
}

IDs <- unique(data$ID)

# remove -668 population size data point
data <- data[ data$PopBio > -600, ]



# rescale population sizes

data_rescaled <- data

for (id in IDs){
  if (sum(data[data$ID == id,]$PopBio <= 0) != 0){
    data_rescaled[data_rescaled$ID == id,]$PopBio <- data_rescaled[data_rescaled$ID == id,]$PopBio - min(data_rescaled[data_rescaled$ID == id,]$PopBio)
  }
  else{
    
  }
}

data_rescaled <- data_rescaled[ data_rescaled$PopBio > 0, ]

# add log population sizes
data_log <- data_rescaled
data_log$logPopBio <- log10(data_rescaled$PopBio)

for (id in IDs){
  if (sum(data_log$ID == id) < 5){
    data_log <- data_log[data_log$ID != id,]
    print(id)
  }
}

# save data
write.table(data_log, file = "../data/LogisticGrowthDataLogClean.csv", 
            sep = "\t", row.names=T, col.names=TRUE)

