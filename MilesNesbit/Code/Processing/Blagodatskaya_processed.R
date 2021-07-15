library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
read.csv("../Data/Blagodatskaya/Blagodatskaya.csv")


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

#GH
GH_df <- read.csv("../Data/Blagodatskaya/Blagodatsk_V2.csv",
                      stringsAsFactors = F) %>%
  format_data(.,"22","Hours","mg*C*(g^-1)","Soil Microbial Community","loamy Luvic Chernozem","header","")


#get temperature
# final_df <- bind_rows(GH_df) %>%
#   mutate(Temp = as.numeric(str_extract(header,regex("[0-9]+(?=_)") )),
#          Rep  = as.numeric(str_extract(header,regex("(?<=_)[0-9]") ))) %>%
#   select(-header)

GH_df$Medium <-GH_df$header
GH_df$Rep <- "1"
GH_df$Citation <- "Blagodatskaya, E. V, Blagodatsky, S. A., Anderson, T.-H., & Kuzyakov, Y. (2007). Priming effects in Chernozem induced by glucose and N in relation to microbial growth strategies. Applied Soil Ecology, 37(1), 95–105. https://doi.org/https://doi.org/10.1016/j.apsoil.2007.05.002"
drops <- c("header")
Final_df <- GH_df[ , !(names(GH_df) %in% drops)]

Final_df %>%
  ggplot(aes(x=Time,y=PopBio, colour = Temp, group = Temp))+
  geom_point()+
  geom_line()+
  facet_wrap(~Species)

#Save
write.csv(Final_df,"../Results/Formated/Blagodatskaya_processed.csv")
