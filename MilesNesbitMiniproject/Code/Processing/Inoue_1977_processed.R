library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
read.csv("../Data/Inoue_1977*/Fig_1.csv")
fx<- read.csv("../Data/Inoue_1977*/Fig_1.csv")

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

#Curtobacterium_psychrophilum
Fig_1 <- read.csv("../Data/Inoue_1977*/Fig_1.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Curtobacterium_psychrophilum","Peptone",1)

#Cytophaga_antarctica
Fig_2<- read.csv("../Data/Inoue_1977*/Fig_2.csv",
                       stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Cytophaga_antarctica","Peptone",1)

#Cytophaga_xantha
Fig_3<- read.csv("../Data/Inoue_1977*/Fig_3.csv",
                 stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Cytophaga_xantha","Peptone",1)

#Spirillum_pleomorphum
Fig_4<- read.csv("../Data/Inoue_1977*/Fig_4.csv",
                 stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Spirillum_pleomorphum","Peptone",1)

#Micrococcus_cryophilus
Fig_5<- read.csv("../Data/Inoue_1977*/Fig_5.csv",
                 stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Micrococcus_cryophilus","Peptone",1)

#Pseudomonas_flourescens
Fig_6<- read.csv("../Data/Inoue_1977*/Fig_6.csv",
                 stringsAsFactors = F) %>%
  format_data(.,"header","Hours","Absorbance (660 nm)","Pseudomonas_flourescens","Peptone",1)

#Squish it all together and clean it up
final_df <- bind_rows(Fig_1, Fig_2, Fig_3, Fig_4, Fig_5, Fig_6) %>%
  mutate(Temp = header)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]
final_df[final_df==".0.4"]<-"-0.4"
final_df[final_df==".0.2"]<-"-0.2"

#cite
final_df$Citation <- "Inoue, K. (1977). Effect of temperature on growth of obligately psychrophilic bacteria. The Journal of General and Applied Microbiology, 23(2), 53–63. https://doi.org/10.2323/jgam.23.53"

final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#Wrap it up
write.csv(final_df,"../Results/Formated/Inoue_1977_processed.csv")
