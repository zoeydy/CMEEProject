
rm(list = ls())
graphics.off()
setwd("~/Documents/Project/code")
require(ggplot2)


# read the data
Data <- read.csv('../data/growth_rate_data.csv')
# id is generated from:
# Data$ID <- paste(Data$Species,"_",Data$Temp,"_",Data$Medium,"_",Data$Rep,"_",Data$Citation)


# # plot the raw data
# for (i in 1:length(unique(Data$id))) {
# 
#   data <- subset(Data, Data$id == i)
#   # plot
#   fileName <- paste0('../results/RawMetaPlot/',i,'.png')
#   png(fileName)
#   ggplot(data , aes(x = Time, y = PopBio)) +
#     geom_point(size = 2) +
#     xlab('Time (Hours)') +
#     ylab('Population Size')
#   ggsave(
#     filename = fileName,
#     width=10,
#     height=10
#   )
#   graphics.off()
# }
# problematic data sets(id) by checking raw plot
# 5,14,17,20,21,104:115,116:119,232,254,260,261,266,279,297,317,330:333,350:353,364:367,414,419,423,424,435,440,544:546,554,570,743,760,762,776:779,784,850,851


# checking problematic data set(id) has negative PopBio
neg <- subset(Data, Data$PopBio < 0)
neg.id <- unique(neg$id)
neg.df <- data.frame()
for (j in neg.id) {
  negdf <- Data[Data$id == j, ]
  neg.df <- rbind(neg.df, negdf)
}
unique(neg.df$id) == unique(neg$id)
check.id <- unique(neg.df$id)
df <- neg.df[neg.df$id == check.id[4],]
ggplot(data = df, aes(x = Time, y = PopBio)) +
  geom_point()
check.id
########## by checking 93 has one negative PopBio value, data set 5 17 21 should be just deleted
# delete
Data <- Data[Data$id != 5,]
Data <- Data[Data$id != 17,]
Data <- Data[Data$id != 21,]
Data <- Data[Data$PopBio >= 0,]
length(unique(Data$id))
# delete point has 0 value of PopBio
Data <- Data[Data$PopBio !=0, ]
# delete problematic data
Data <- Data[Data$id != 544,]
Data <- Data[Data$id != 545,]
Data <- Data[Data$id != 546,]
Data <- Data[Data$id != 762,]
Data <- Data[Data$id != 851,]
# add logPopBio
Data$logPopBio <- log10(Data$PopBio)


# save data
# write.table(Data, file = "../data/processed_data.csv", 
#             sep = "\t", row.names=T, col.names=TRUE)
write.csv(Data, "../data/processed_data.csv", )
colnames(Data)[12] <- 'ID'

# checking data
# str(Data)
# unique(Data$ID)
# length(unique(Data$ID))
# unique(Meta$ID)
# length(unique(Meta$ID))
# unique(unique(Data$PopBio) >0)
# unique(probPopBio$ID)


