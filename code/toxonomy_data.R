rm(list = ls())
library(taxize)
library(plyr)

data <- read.csv('../data/pop.csv')


origin_species <- unique(data$Species)
correct_species <- gnr_resolve(names = origin_species)
correct_species <- as.data.frame(correct_species)


temp <- gnr_resolve(names = c('Helianthos annus', 'Homo saapiens'))
temp[,-c(1,4)]
classification(origin_species, db = 'itis')

for (i in origin_species) {
  tryCatch(
    {
      a <- gnr_resolve(i)
      b <- classification(a$matched_name[1], db = 'ncbi')
      c <- b[[1]]
      d <- c[7:10,]
      test <- rbind(test, d)
      print(paste('Value of d :',d))
    }
  )
}
