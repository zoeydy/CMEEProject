library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
read.csv("../Data/Heo_et_al_2009*/Fig_1.csv")


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

#5 degrees
Temp_5_df <- read.csv("../Data/Heo_et_al_2009*/Fig_1.csv",
                      stringsAsFactors = F) %>%
  format_data(.,5,"Hours","CFU","header","frankfurters",1)

#15 degrees
Temp_15_df <- read.csv("../Data/Heo_et_al_2009*/Fig_2.csv",
                       stringsAsFactors = F) %>%
  format_data(.,15,"Hours","CFU","header","frankfurters",1)

#25 degrees
Temp_25_df <- read.csv("../Data/Heo_et_al_2009*/Fig_3.csv",
                       stringsAsFactors = F) %>%
  format_data(.,25,"Hours","CFU","header","frankfurters",1)




#get bits of data
final_df <- bind_rows(Temp_5_df,Temp_15_df,Temp_25_df) %>%
  mutate(Species = header) #%>%
  #select(-header)
#Fix time units
final_df <- final_df %>%
  mutate(Time = Time * 24,
         Time_units = "Hours")
#get rid of extra header
drops <- c("header")
Final_df <- final_df[ , !(names(final_df) %in% drops)]

#Cite
Final_df$Citation <-"Heo, C., Choi, Y. S., Kim, C. J., & Paik, H. D. (2009). Estimation of shelf-life of frankfurter using predictive models of spoilage bacterial growth. Korean Journal for Food Science of Animal Resources, Vol. 29, pp. 289–295. https://doi.org/10.5851/kosfa.2009.29.3.289"

Final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

write.csv(Final_df,"../Results/Formated/Heo_et_al_processed.csv")

