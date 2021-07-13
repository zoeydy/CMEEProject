library(tidyverse)
library(ggplot2)
rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
read.csv("../Data/Ratkovsky/Fig_1.csv")
fx<- read.csv("../Data/Ratkovsky/Fig_1.csv")

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

#Get data
Fig_1 <- read.csv("../Data/Ratkovsky/Fig_1.csv",
                  stringsAsFactors = F) %>%
  format_data(.,"header","Hours"," Logarithm of growth rate constant","Curtobacterium_psychrophilum","Chicken",1)

final_df <- bind_rows(Fig_1) %>%
  mutate(Species = header) %>%
  mutate(Temp = Time) %>%
  mutate(Time = NA)
drops <- c("header")
final_df <- final_df[ , !(names(final_df) %in% drops)]

#give up on this one
final_df %>%
  ggplot(aes(x=Temp,y=PopBio, colour = Species, group = Temp))+
  geom_point()+
  geom_line()

#cite
final_df$Citation <- "Ratkowsky, D. A., Olley, J., McMeekin, T. A., & Ball, A. (1982). Relationship between temperature and growth rate of bacterial cultures. Journal of Bacteriology, 149(1), 1. Retrieved from http://jb.asm.org/content/149/1/1.abstract"

#Save it I guess, but with a caviat
write.csv(final_df,"../Results/Formated/Ratkowsky***_processed.csv")
