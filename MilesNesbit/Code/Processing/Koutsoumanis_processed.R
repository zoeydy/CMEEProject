library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")


#data function
format_data <- function(df,temp,time_units,popBio_unts,species,medium,rep,header){
  #get number of datasets
  n_data <- ncol(df)/2
  data_list <- list()
  for(i in 1:n_data){
    #construct dataframe
    df_temp <- df[ , ((i*2)-1) : (i*2) ] %>%
      filter(row_number() > 1) %>%
      mutate_all(as.numeric)
    
    df_temp[df_temp == ""] <- NA
    df_temp <- df_temp %>% na.omit()
    
    #get header
    header <- colnames(df_temp[1]) %>% gsub("X","",.)
    #add data 
    colnames(df_temp) <- c("Time","PopBio")
    df_temp$header <- header
    df_temp$Temp <- temp
    df_temp$Time_units <- time_units
    df_temp$PopBio_unts <- popBio_unts
    df_temp$Species <- species
    df_temp$Medium <- medium
    df_temp$Rep <- rep
    
    data_list[[i]] <- df_temp
  }
  return(bind_rows(data_list))
}

#gather data
Fig_1 <- read.csv("../Data/Koutsoumanis/Fig_1_beef_0.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"0","Hours","CGU","","beef",1)

Fig_2 <- read.csv("../Data/Koutsoumanis/Fig_2_beef_10.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"10","Hours","CGU","","beef",1)

Fig_3 <- read.csv("../Data/Koutsoumanis/Fig_3_pork_0.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"0","Hours","CGU","","pork",1)

Fig_4 <- read.csv("../Data/Koutsoumanis/Fig_4_pork_10.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"10","Hours","CGU","","pork",1)

#squish
final_df <- bind_rows(Fig_1, Fig_2, Fig_3, Fig_4) %>%
  mutate(Species = header)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#clean
final_df[final_df=="pseodonomads"]<-"pseudonomads"

#check
final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#cite
final_df$Citation <- 'Koutsoumanis, K., Stamatiou, A., Skandamis, P., & Nychas, G. J. E. (2006). Development of a microbial model for the combined effect of temperature and pH on spoilage of ground meat, and validation of the model under dynamic temperature conditions. Applied and Environmental Microbiology, 72(1), 124–134. https://doi.org/10.1128/AEM.72.1.124-134.2006'

#save
write.csv(final_df,"../Results/Formated/Koutsoumanis_processed.csv")
