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
Fig_1 <- read.csv("../Data/Sato-Takabe_et_al_2019/Fig_1.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"","Hours","Total Bacterial Abundance (x10^6 cells/mL)","anoxygenic phototrophic bacteria","seawater",1)

#squish
final_df <- bind_rows(Fig_1) %>%
  mutate(Temp = header)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#check
final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#cite
final_df$Citation <- 'Sato-Takabe, Y., Hamasaki, K., & Suzuki, S. (2019). High temperature accelerates growth of aerobic anoxygenic phototrophic bacteria in seawater. MicrobiologyOpen, 8(5), 1–6. https://doi.org/10.1002/mbo3.710'

#save
write.csv(final_df,"../Results/Formated/Sato-Takabe_processed.csv")
