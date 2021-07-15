require("minpack.lm")  # for Levenberg-Marquardt nlls fitting
library(tidyverse)
library(ggplot2)
library(dplyr)
library(growthcurver)
library(stringr)
library(taxize)
library(nls.multstart)
library(zoo)
library(broom)
library(growthcurver)
library(cowplot)
library(psych)
library(MASS)
library(fitdistrplus)
library(directlabels)
library(lmodel2)

rm(list=ls())
graphics.off()
setwd("/home/nesbit/Desktop/Data_RVK/Code")
files <- read.csv("../Results/Formated/5Total_processed_Log_RVK_3_models.csv")
# fa <- read.csv('../Results/Formated/Total_processed.csv')
# 
# specieslist<- unique(fa$Species)
# write.csv(specieslist, '../Results/Formated/UnformattedSpeciesList.csv')

#clean the data
files$r <- NA
files$K <- NA
files$t_lag <- NA
files$N_0 <- NA
files$r.t.value <- NA
files$K.t.value <- NA
files$t_lag.t.value <- NA
files$N_0.t.value <- NA

for (i in unique(files$discrete_curve)){
  
  # Subset data by curve number
  ff_sub <- subset(files, discrete_curve == i)
  
  # Find estimated values
  r <- subset(ff_sub, parameter == "r_max")$Estimate[1]
  r_tval <- subset(ff_sub, parameter == "r_max")$t.value[1]
  Nmax <- subset(ff_sub, parameter == "N_max")$Estimate[1]
  Nmax_tval <- subset(ff_sub, parameter == "N_max")$t.value[1]
  N0 <- subset(ff_sub, parameter == "N_0")$Estimate[1]
  N0_tval <- subset(ff_sub, parameter == "N_0")$t.value[1]
  tlag <- subset(ff_sub, parameter == "t_lag")$Estimate[1]
  tlag_tval <- subset(ff_sub, parameter == "t_lag")$t.value[1]
  
  # Add data 
  index <- which(files$discrete_curve == i)
  for (j in index){
    files$r[j] <- r
    files$r.t.value[j] <- r_tval
    files$K[j] <- Nmax
    files$K.t.value[j] <- Nmax_tval
    files$N_0[j] <- N0
    files$N_0.t.value[j] <- N0_tval
    files$t_lag[j] <- tlag
    files$t_lag.t.value[j] <- tlag_tval
  }
  
}

#nest
ff <- files %>%
  nest(-Temp,-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts, -Estimate,
       -t.value, -discrete_curve, -model, -parameter, -r, -K, -t_lag, -N_0,-r.t.value,
       -K.t.value,-N_0.t.value,-t_lag.t.value,-Nutrient.Content,- spp_species,- spp_genus,
       - spp_genus,-spp_class,-spp_family, data = c("Time","PopBio"))

# Tidy ff
ff_new <- subset(ff, parameter == "r_max")


#drops
drops <- c("parameter", "t.value", "Estimate")
df <- ff_new[ , !(names(ff_new) %in% drops)]


#fit boltzmann curves
k <- 8.617e-5


params_df <- df %>%
  dplyr::select(-data) %>% #remove extra data
  mutate(a = r / K) %>% #calculate a_ii
  gather("Param","Val",r,K,a) %>% #combine r and K to long format
  mutate(Temp = (1 / (k * (Temp+273.15))) - (1 / (k * (15+273.15))) ,Val = log(Val)) %>% #transform for boltz
  filter(!is.infinite(Val)) %>% #remove inf
  group_by(Citation,Param,Species,Medium,Rep) %>% mutate(n = n()) %>% ungroup() %>%
  filter(n > 2) %>%
  nest(-Species,-Medium,-Rep,-Citation,-Time_units,-PopBio_unts,-Param,-Nutrient.Content, -spp_class, -spp_species) %>% #nest for boltz fitting
  mutate(fits = map(data, ~lm(.$Val ~ .$Temp)), #fit models
         tidied = map(fits,tidy)) %>% #tidy output
  unnest(tidied) %>% #unnest
  mutate(term = recode(term,`(Intercept)` = "B0", `.$Temp` = "E")) %>% #recode the parameter estimates
  #filter(p.value < .05) %>% #cut-off for p-values
  dplyr::select(-std.error,-statistic,-p.value) %>%
  spread(term,estimate)  %>%
  filter(!is.na(B0),!is.na(E))




p1 <- params_df %>%
  ggplot(aes(x=E,group = PopBio_unts, colour = PopBio_unts))+
  geom_density()+
  facet_wrap(~Param)+
  geom_vline(xintercept = 0)

p2 <- params_df %>%
  ggplot(aes(x=B0,group = PopBio_unts, colour = PopBio_unts))+
  geom_density()+
  facet_wrap(~Param)+
  geom_vline(xintercept = 0)

#distributions
plot_grid(p1,p2,nrow = 2)


#CFU B0
CFU <-  filter(params_df, PopBio_unts =='N' |PopBio_unts =='CFU'|PopBio_unts =='CFU/g')


CFU2 <-CFU %>% unnest(data)

#plot b0 CFU
CFU4 <- CFU %>% 
  filter(Param != "a") %>%
  dplyr::select(-data,-fits,-E) %>%
  spread(Param,B0) %>%
  filter(!is.na(K) & !is.na(r))

#rma
lmodel1 <- lmodel2(K~r, data = CFU4, range.y = 'interval', range.x = 'interval')
rma<- list()
T_rma1 = seq(min(CFU4$r), max(CFU4$r), length.out = 100)
rma1 = data.frame(rlog =T_rma1,
                 rmapred =(0.2734124	*T_rma1) +		1.880166 )

  ggplot(CFU4, aes(x=r,y=K))+
  geom_point(size= 3)+
  # geom_smooth(method='lm', formula= y~x, se= F)+
  ggtitle('Log r and K by B0')+
  xlab('Log of r')+
  ylab('Log of K')+
  geom_line(data = rma1, aes(x= rlog, y = (rmapred)), color = 'red', size = 3)
  
  
  

#lognormal skew
ln_skew <- function(v){
  (exp(v)+2) * sqrt(exp(v)-1)
}

params_df %>%
  group_by(Param,PopBio_unts) %>%
  summarise(n = n(),
            uE = mean(E),vE=var(E),kE=skew(E),
            uB0= mean(B0), vB0 = var(B0),kB0 = skew(B0)) %>%
  gather("measure","value",uE,vE,kE,uB0,vB0,kB0) %>%
  ggplot(aes(x=interaction(measure,Param),y=value,colour=PopBio_unts))+
  geom_point()

#r
x <- params_df %>% filter(Param == "a")

fit_B0 <- fitdistr(x$B0,"normal")
B0_pred = rnorm(1000,fit_B0$estimate["mean"],fit_B0$estimate["sd"])

fit_E <- fitdistr(x$E,"normal")
E_pred = rnorm(1000,fit_E$estimate["mean"],fit_E$estimate["sd"])

p1 <- data.frame(B0 = x$B0) %>%
  ggplot(aes(x=B0))+
  geom_density(colour = "red") +
  geom_density(data = data.frame(B0 = B0_pred),colour = "blue")

p2 <- data.frame(E = x$E) %>%
  ggplot(aes(x=E))+
  geom_density(colour = "red") +
  geom_density(data = data.frame(E = E_pred),colour = "blue")

plot_grid(p1,p2)

descdist(x$E,discrete=FALSE, boot=500)




#filter out 0
df2 <- filter(df, K > 0)

drops2 <- c("data")
df2 <- df2[ , !(names(df2) %in% drops2)]

515/876

#try the plot
Absorbance <- df2 %>%
  filter(str_detect(PopBio_unts, 'Abs'))

ABSD <- df2 %>%
  filter(str_detect(PopBio_unts, 'OD'))

Absorbance <- rbind(Absorbance, ABSD)

# Absorbance <- filter(Absorbance, K < 150)
# 
# Absorbance <- filter(Absorbance, r < 20)
  

#create log
Absorbance$Klog <- log(Absorbance$K)
Absorbance$rlog <- log(Absorbance$r)

#filter
Absorbance <- filter(Absorbance, Klog >  -2)
Absorbance <- filter(Absorbance, rlog >  -20)

#plot
ggplot(Absorbance, aes(x= rlog, y= Klog, color = Temp))+
  geom_point(size = 2)+ 
  ggtitle("Absorbance")+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~spp_species, scales = 'free')


#CFU
CFU <- df2 %>%
  filter(str_detect(PopBio_unts, 'CFU'))
CFU2 <- df2 %>%
  filter(str_detect(PopBio_unts, 'N'))
CFU3 <- df2 %>%
  filter(str_detect(PopBio_unts, 'CGU'))

CFU <- rbind(CFU, CFU2, CFU3)
#create log
CFU$Klog <- log(CFU$K)
CFU$rlog <- log(CFU$r)

#filter
CFU <- filter(CFU, Klog >  -4)
CFU <- filter(CFU, rlog >  -20)
CFU2 <- CFU %>%
  filter(!str_detect(spp_class, 'NA'))

#rma
lmodel <- lmodel2(Klog~rlog, data = CFU2, range.y = 'interval', range.x = 'interval')
rma<- list()
T_rma = seq(min(CFU2$rlog), max(CFU2$rlog), length.out = 100)
rma = data.frame(rlog =T_rma,
                 rmapred =(0.3543012*T_rma) +	2.354129 )


#plot
ggplot(CFU2, aes(x= rlog, y= Klog, color = Temp))+
  geom_point(size = 5)+ 
  scale_colour_gradientn(colours=rainbow(4))+
  ggtitle("CFU")+
  xlab('log of r')+
  ylab('log of K')+
  # geom_smooth(method='lm', formula= y~x, se = F, size = 3)+
  #facet_wrap(~spp_species, scales = 'free')
  geom_line(data = rma, aes(y = (rmapred)), color = 'red', size = 3)


##r changes with temp but not K
##real analysis
TPC_fits <- params_df %>% 
  dplyr::select(-data,-fits) %>%
  pivot_wider(names_from = Param, values_from = c(B0,E)) %>%
  filter(PopBio_unts =='N' |PopBio_unts =='CFU'|PopBio_unts =='CFU/g') %>%
  group_by(spp_species) %>%
  summarise(B0_r = mean(B0_r), B0_K = mean(B0_K), E_r = mean(E_r))

T_vec = seq(min(CFU$Temp), max(CFU$Temp), length.out = 100)

res <- list()
for(i in 1:nrow(TPC_fits)){
  data <- TPC_fits[i,]
  res[[i]] <- data.frame(Temp = T_vec,
                 a0 = exp(data$B0_r - data$B0_K),
                 aE = 2*data$E_r,
                 log_a_pred = (data$B0_r - data$B0_K) + (2*data$E_r * ((1 / (k * (T_vec+273.15))) - (1 / (k * (15+273.15))))) ,
             spp_species = as.character(data$spp_species))
}

CFU %>%
  mutate(a = r / K)  %>%
  ggplot(aes(Temp,log(a))) +
        geom_point(size= 3)+
        facet_wrap(~spp_species,scales="free")+
       # geom_hline(aes(yintercept = a0), data = bind_rows(res))
        geom_line(data = bind_rows(res), aes(y = (log_a_pred)), color = 'red')+
        geom_smooth(method='lm', formula= y~x)+  
  ggtitle("Comparison of Actual Line by Predicted Line")+
  xlab('Temperature')+
  ylab('log of a')
  

#example graphs
fakedata <- read.csv("../Results/Formated/FakeData1.csv")
fakedata$Fake_Data_Y <- as.character(fakedata$Fake_Data_Y)
fakedata$Fake_Data_Y <- as.integer(fakedata$Fake_Data_Y)
fakedata2 <- fakedata %>%
  filter(Hypothesis == 'MTE Tradeoff'| Hypothesis == 'K not temperature dependent'| Hypothesis == 'No tradeoff')

#plot
ggplot(fakedata2, aes(x= Fake_Data_X, y= Fake_Data_Y, color = Style))+
  # geom_point(size = 2)+ 
  facet_wrap(~Hypothesis)+
  ggtitle("Predicted Graphs from Different Hypothesis")+
  xlab(' ')+
  ylab(' ')+
  geom_line(aes(x= Fake_Data_X, y= Fake_Data_Y, color = r.K, size = 2))+
  guides(size = 'none')+
  geom_dl(aes(label = r.K),color = 'black', method = list(dl.combine("last.points"), cex = 2 ))+
  theme(plot.margin = unit(c(1,3,1,1), "lines"))+
  theme(legend.position = 'none')

hjust = -0.5
fakedata3 <- read.csv("../Results/Formated/FakeData2.csv")
fakedata3$Fake_Data_Y <- as.character(fakedata3$Fake_Data_Y)
fakedata3$Fake_Data_Y <- as.integer(fakedata3$Fake_Data_Y)

#plot
ggplot(fakedata3, aes(x= Fake_Data_X, y= Fake_Data_Y, color = r.K))+
  geom_point(size = 2)+ 
  facet_wrap(~Hypothesis)+
  ggtitle("Predicted Graphs from Different Hypothesis")+
  xlab(' ')+
  ylab(' ')+
  geom_smooth(method='lm', formula= y~x, se =F, size =5)+ 
  # geom_line(aes(x= Fake_Data_X, y= Fake_Data_Y, color = r.K, size = 5))+
  guides(size = 'none')

#analysis chart
##real analysis
TPC_fits <- params_df %>% 
  dplyr::select(-data,-fits) %>%
  pivot_wider(names_from = Param, values_from = c(B0,E)) %>%
  filter(PopBio_unts =='N' |PopBio_unts =='CFU'|PopBio_unts =='CFU/g') %>%
  group_by(spp_species) %>%
  summarise(B0_r = mean(B0_r), B0_K = mean(B0_K), E_r = mean(E_r))

rely <- data.frame(T_min=as.numeric(c()), T_max=as.numeric(c()), spp_species=as.character(c()), stringsAsFactors = FALSE)
for (j in unique(CFU$spp_species)) {
  CFI <- subset(CFU, spp_species == j)
  T_vec <- data.frame(T_min = min(CFI$Temp), T_max = max(CFI$Temp), spp_species = j)
  rely = rbind(T_vec, rely)
}


TPC_fits <- TPC_fits[ !grepl("NA yeast", TPC_fits$spp_species) , ]
res <- list()
for(i in 1:nrow(TPC_fits)){
  data <- TPC_fits[i,]
  minmax = subset(rely, spp_species == unique(as.character(data$spp_species)))
  print(as.character(minmax$spp_species))
  T_min = minmax$T_min[1]
  T_max = minmax$T_max[1]
  
  T_seq = seq(T_min, T_max, length.out=100) 
  res[[i]] <- data.frame(Temp = T_seq,
                         a0 = exp(data$B0_r - data$B0_K),
                         aE = 2*data$E_r,
                         log_a_pred = (data$B0_r - data$B0_K) + (2*data$E_r * ((1 / (k * (T_seq+273.15))) - (1 / (k * (15+273.15))))) ,
                         spp_species = as.character(data$spp_species))
  print(unique(as.character(data$spp_species)))
}



estimated_tt <- data.frame(spp_species=as.factor(c()),

                           
CFU %>%
  mutate(a = r / K)  %>%
  ggplot(aes(Temp,log(a))) +
  loga=as.numeric(c()),
  Temp2=as.numeric(c()))
for (spp_df in res){
  Temp2 <- (1 / (k * (spp_df$Temp+273.15))) - (1 / (k * (15+273.15)))
  loga <- log(spp_df$a0)
  estimated_tt <- rbind(estimated_tt, data.frame(Temp2, loga, spp_species=spp_df$spp_species))
}

CFU <- CFU[ !grepl("NA PCT", CFU$spp_species) , ]
CFU <- CFU[ !grepl("Brochothrix thermosphacta", CFU$spp_species) , ]
CFU <- CFU[ !grepl("NA yeast", CFU$spp_species) , ]
CFU %>%
  mutate(a = r / K)  %>%
  ggplot(aes(Temp,log(a))) +
  geom_point(size= 3)+
  facet_wrap(~spp_species,scales="free")+
  # geom_hline(aes(yintercept = a0), data = bind_rows(res))
  geom_line(data = bind_rows(res), aes(y = (log_a_pred)), color = 'red', size = 1.5)+
  geom_smooth(method = 'lm', se = FALSE, size = 1)+
  ggtitle("Comparison of Actual Line by Predicted Line")+
  xlab('Temperature')+
  ylab('log of a')


tt <- CFU %>%
  mutate(a = r / K) %>%
  mutate(loga = log(a)) %>%
  mutate(Temp2 = (1 / (k * (Temp+273.15))) - (1 / (k * (15+273.15))))

#analysis chart
z <- list()

species_name_list <- unique(tt$spp_species)
for (i in species_name_list){
  fsummary = lm(loga ~ Temp2, data = tt, subset = spp_species == i)
  gsummary <- as.data.frame(summary(fsummary)$coeff)
  gsummary$species = i
  z <- rbind(gsummary, z)
}


# Get df of log(a) and normalised 1/kt & spp values
tt_2 <- tt %>% dplyr::select(spp_species, Temp2, loga)

# Lm to get slope
observed_data_fits <- list()
for (i in unique(tt_2$spp_species)){
  fsummary = lm(loga ~ Temp2, data = tt_2, subset = spp_species == i)
  gsummary <- as.data.frame(summary(fsummary)$coeff[2])
  gsummary$species = i
  observed_data_fits <- rbind(gsummary, observed_data_fits)
}
colnames(observed_data_fits) <- c("slope_obs", "species")
# end up with df of slope & spp

#get slope difference
final_df <- bind_rows(res) %>%
  group_by(spp_species) %>%
  summarise(slope_exp = mean(aE)) %>%
  rename(species = spp_species) %>%
  dplyr::left_join(., observed_data_fits , by="species") %>%
  mutate(slope_diff = slope_obs - slope_exp)



ggplot(final_df, aes(x= slope_diff))+
  geom_density(color= 'red')+
  ggtitle('Observed - Expected E Distribution')+
  xlab('Observed -Expected E')+
  geom_vline(xintercept=c(0), color = 'blue')



#done
