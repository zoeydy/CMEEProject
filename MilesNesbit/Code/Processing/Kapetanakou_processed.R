library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
read.csv("../Data/Kapetanakou_et_al_2019/Fig_1_10C.csv")

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


#0C
Fig_2_0C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_2_0C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"0C","Hours","CFU/g","","Eruca_satica L.",1)

Fig_4_0C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_4_0C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"0C","Hours","CFU/g","","Eruca_satica L.",1)

#5C
Fig_2_5C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_2_5C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"5C","Hours","CFU/g","","Eruca_satica L.",1)

Fig_4_5C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_4_5C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"5C","Hours","CFU/g","","Eruca_satica L.",1)

#10C
Fig_1_10C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_1_10C.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"10C","Hours","CFU/g","","Eruca_satica L.",1)

Fig_3_10C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_3_10C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"10C","Hours","CFU/g","","Eruca_satica L.",1)

#15C
Fig_1_15C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_1_15C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"15C","Hours","CFU/g","","Eruca_satica L.",1)

Fig_3_15C <- read.csv("../Data/Kapetanakou_et_al_2019/Fig_3_15C.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"15C","Hours","CFU/g","","Eruca_satica L.",1)

#squish it all together
final_df <- bind_rows(Fig_3_15C,Fig_1_15C,Fig_3_10C,Fig_1_10C,Fig_4_5C,Fig_2_5C,Fig_4_0C,Fig_2_0C) %>%
  mutate(Species = header)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#Time
final_df <- final_df %>%
  mutate(Time = Time * 24,
         Time_units = "Hours")
#check
final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#Clean
final_df[final_df=="enterobacteriaceae_sp_5._O2"]<-"enterobacteriaceae_5._O2"
final_df[final_df=="enterobactericeae_5._O2"]<-"enterobacteriaceae_5._O2"

#Check
final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#Cite
final_df$Citation <- "Kapetanakou, A. E., Taoukis, P., & Skandamis, P. N. (2019). Model development for microbial spoilage of packaged fresh‒cut salad products using temperature and in-package CO2 levels as predictor variables. LWT, 113, 108285. https://doi.org/https://doi.org/10.1016/j.lwt.2019.108285"

#clean a little more
final_df[final_df== "0C"] <- "0"
final_df[final_df== "10C"] <- "10"
final_df[final_df== "5C"] <- "5"
final_df[final_df== "15C"] <- "15"

final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#Wrap it up
write.csv(final_df,"../Results/Formated/Kapetanakou_processed.csv")
