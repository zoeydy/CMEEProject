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
Fig_1 <- read.csv("../Data/Willocx_et_al_1993/Fig_1.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"8","Hours","CFU/g","P. flourescens","tryptone soya agar/broth",1)

#squish
final_df <- bind_rows(Fig_1) %>%
  mutate(Temp = header)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#temp units
final_df <- final_df %>%
  mutate(Temp = as.numeric(Temp) - 273.15)
         
#check
final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#cite
final_df$Citation <- 'Willocx, F., Mercier, M., Hendrickx, M., & Tobback, P. (1993). Modelling the influence of temperature and carbon dioxide upon the growth of Pseudomonas fluorescens. Food Microbiology, Vol. 10, pp. 159–173. https://doi.org/10.1006/fmic.1993.1016'

#save
write.csv(final_df,"../Results/Formated/Wilocx_processed.csv")
