
rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")

require(ggplot2)
library(ggpubr)


#################
# read the data #
#################
Data <- read.csv('../data/pop.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]
IDs <- unique(Data$id)

# plot data
plot.df <- read.csv("../data/gomp.plot.csv")

# info
infos <- read.csv("../data/gomp.info.csv")

info0 <- subset(infos, infos$rmax > 0 & infos$tlag > 0) 
info0[, "temp_c"] <- round(info0[, "temp_c"], digits = 2)

info0$log.one.t <- log(1/info0$tlag)
info0$log.r <- log(info0$rmax)




# set boltzmann constant 8.617*10^(-5)
K = 8.617*10^(-5)

save_plot <- function(filepath, plo){
  pdf(paste0(filepath,".pdf"))
  print(plo)
  graphics.off()
}
library(RColorBrewer)
# gene_col <- function(n){
#   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#   return(as.character(sample(col_vector, n)))
# }



# devide the temperature into 4 groups ((-5:10, 10:20, 20:30, 30:40))
df1 <- subset(info0, info0$temp_c >= -5 & info0$temp_c < 10)
df2 <- subset(info0, info0$temp_c >= 10 & info0$temp_c < 20)
df3 <- subset(info0, info0$temp_c >= 20 & info0$temp_c < 30)
df4 <- subset(info0, info0$temp_c >= 30 & info0$temp_c < 40)
df1$temp_c_group <- "-5~10 °C"
df2$temp_c_group <- "10~20 °C"
df3$temp_c_group <- "20~30 °C"
df4$temp_c_group <- "30~40 °C"
info4 <- rbind(df1,df2,df3,df4)
info4$temp_c_group <- factor(info4$temp_c_group, levels = unique(info4$temp_c_group))
# devide the temperature into 5 groups ((-5:5, 5:10, 10:20, 20:30, 30:40))
dm1 <- subset(info0, info0$temp_c >= -5 & info0$temp_c < 5)
dm2 <- subset(info0, info0$temp_c >= 5 & info0$temp_c < 10)
dm3 <- subset(info0, info0$temp_c >= 10 & info0$temp_c < 20)
dm4 <- subset(info0, info0$temp_c >= 20 & info0$temp_c < 30)
dm5 <- subset(info0, info0$temp_c >= 30 & info0$temp_c < 40)
dm1$temp_c_group <- "-5~5 °C"
dm2$temp_c_group <- "5~10 °C"
dm3$temp_c_group <- "10~20 °C"
dm4$temp_c_group <- "20~30 °C"
dm5$temp_c_group <- "30~40 °C"
info5 <- rbind(dm1,dm2,dm3,dm4,dm5)
info5$temp_c_group <- factor(info5$temp_c_group, levels = unique(info5$temp_c_group))

# # color
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# correlation between rmax and 1/tlag
p <- ggplot(data = info0, aes(x = log.one.t, y = log.r, colour = temp_c)) +
  geom_point(size = 1) +
  theme_set(theme_bw()) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = NA), panel.border = element_blank(),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17)) +
  stat_smooth(method = 'loess', formula = 'y~x')
  # stat_smooth(formula = y~x, method = lm, fullrange= TRUE,se = TRUE)
save.plot("../results/rt_plot/log_rt_col_temp",p)
p <- ggplot(data = info4, aes(x = log.one.t, y = log.r)) +
  geom_point(size = 1) +
  stat_smooth(formula = y~x, method = lm, se = TRUE)
save.plot("../results/rt_plot/log_rt_temp",p)
p <- ggplot(data = info4, aes(x = log.one.t, y = log.r, colour = temp_c_group)) +
  geom_point(size = 1) +
  stat_smooth(formula = y~x, method = lm, se = TRUE) 
  # annotate(label = paste0("temperature = ", infos$temp_group[1]), geom = "text",
  #          x = min(infos$tlag,na.rm = TRUE), y = max(infos$rmax, na.rm = TRUE),
  #          hjust = 0, vjust = 0,
  #          size = 2.5)
save.plot("../results/rt_plot/log_rt_temp_group",p)


# plot rmax and 1/tlag V.S. temperature
p <- ggplot(data = info4, aes(x = temp_c, y = log.r, colour = temp_c_group)) +
  geom_point(size = 1) +
  stat_smooth(formula = y~x, method = lm, se = TRUE)
save.plot("../results/rt_plot/log_rmax_temp_group",p)
#1/tlag
p <- ggplot(data = info4, aes(x = temp_c, y = log.one.t, colour = temp_c_group)) +
  geom_point(size = 1) +
  stat_smooth(formula = y~x, method = lm, se = TRUE) 
# stat_regline_equation(label.y = -7, size = 2) 
# annotate(label = paste0("temperature = ", info$temp_group[1]), geom = "text",
#          x = min(info$Nmax,na.rm = TRUE), y = max(info$tlag, na.rm = TRUE),
#          hjust = 0, vjust = 0,
#          size = 2.5)
save.plot("../results/rt_plot/log_onetlag_temp_group",p)




















####################################
# define plot function
plotrt <- function(dat){
  ggplot(data = dat, aes(x = log(tlag), y = log(rmax))
         ,xlim(min(log(tlag)), max(log(tlag)))
  ) +
    geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
    geom_point()+
    stat_cor(label.y = max(dat$rmax, na.rm = TRUE)+1, size = 2)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
    stat_regline_equation(label.y = max(dat$rmax, na.rm = TRUE)+2.3, size = 2) +
    annotate(label = paste0("temperature = ", dat$temp_group[1]), geom = "text",
             x = min(dat$tlag,na.rm = TRUE), y = max(dat$rmax, na.rm = TRUE),
             hjust = 0, vjust = 0,
             size = 2.5)
}
plotrk <- function(dat){
  ggplot(data = dat, aes(x = log(Nmax), y = log(rmax))
         ,xlim(min(log(Nmax)), max(log(Nmax)))
  ) +
    geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
    geom_point()+
    stat_cor(label.y = max(dat$rmax, na.rm = TRUE)+1, size = 2)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
    stat_regline_equation(label.y = max(dat$rmax, na.rm = TRUE)+2.3, size = 2) +
    annotate(label = paste0("temperature = ", dat$temp_group[1]), geom = "text",
             x = min(dat$Nmax,na.rm = TRUE), y = max(dat$rmax, na.rm = TRUE),
             hjust = 0, vjust = 0,
             size = 2.5)
}
plott_temp <- function(dat){
  ggplot(data = dat, aes(log(dat$tlag))) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
    geom_density(alpha=0.6) +
    # labs(title=paste("log(tlag) histogram plot, temp = ", dat$temp_group),x="log(tlag)", y = "Count") +
    annotate(label = paste0("temperature = ", dat$temp_group[1]), geom = "text", 
             x = mean(dat$tlag), y = 0.4, hjust = 1, vjust = 0, size = 2.5)
}
plotr_temp <- function(dat){
  ggplot(data = dat, aes(log(dat$rmax))) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
    geom_density(alpha=0.6) +
    annotate(label = paste0("temperature = ", dat$temp_group[1]), geom = "text", 
             x = mean(dat$rmax), y = 0.4, hjust = 1, vjust = 0, size = 2.5)
  # labs(title=paste("log(rmax) histogram plot, temp = ", dat$temp_group),x="log(rmax)", y = "Count")
}
plotk_temp <- function(dat){
  ggplot(data = dat, aes(log(dat$Nmax))) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
    geom_density(alpha=0.6) +
    annotate(label = paste0("temperature = ", dat$temp_group[1]), geom = "text", 
             x = mean(dat$Nmax), y = 0.4, hjust = 1, vjust = 0, size = 2.5)
}
# plotrtemp.log <- function(dat){
#   ggplot(data = dat, aes(x = temp, y = log(rmax), colour = temp_group)) +
#     geom_point(size = 1) +
#     stat_smooth(formula = y~x, method = lm, se = TRUE)
# }


#############################################
# get complete information and plot fitting #
#############################################
info <- data.frame()

for (i in 1:length(IDs)){
  # subset info and data by idname
  idname <- IDs[i]
  data <- Data[Data$id == idname,]
  plot.id <- plot.df[plot.df$id == idname, ]
  info.id <- infos[infos$id == idname, ]
  # get info about temperature, species and medium
  info.id$temp <- data$Temp[1]
  info.id$Species <- data$Species[1]
  info.id$Medium <- data$Medium[1]
  # add kelvin temperature
  info.id$kelvin_temp <- info.id$temp + 273
  # rbind the data frame
  info <- rbind(info, info.id)
  
  
  # plot the fitting
  if (is.na(unique(plot.id$plot.point)[1])){
    print(paste("Fit was not successful for ID:",idname))
  }else{

    legend <- paste0('AICc = ', info.id$AICc, '\nAIC = ', info.id$AIC,
                  '\nBIC = ', info.id$BIC, '\nRsquare = ', info.id$rsq)

    # plot fit plot
    FileName <- paste0("../results/fit_gomp/plot_", idname,".png")
    png(file = FileName)
    p <- ggplot(data, aes(x = Time, y = logN)) +
      geom_point(size = 1) +
      labs(x = "Time (h)", y = "Logarithm of the population size (logN)") +
      ggtitle("Gompertz model comparison plot") +
      geom_line(data = plot.id, aes(x = time, y = plot.point), size=1) +
      # theme(legend.position = 'bottom') +
      annotate('text', label = legend, x = min(data$Time), y = min(data$logN), hjust = -.5, vjust = 0)
    # stat_smooth(method = lm, level = 0.95, aes(colour="Cubic")) +
    # scale_colour_manual(name="Model", values=c("darkblue", "darkred", "darkgreen"))
    print(p)
    graphics.off()
  }
}

###############################
# Plot Typical Growth Pattern #
###############################
# 128 765 813 635
typ.data <- subset(Data, Data$id == 635)
plot.id <- plot.df[plot.df$id == 635, ]
info.id <- infos[infos$id == 635, ]
#legend <- paste0('AICc = ', info.id$AICc, '\nAIC = ', info.id$AIC,
#                 '\nBIC = ', info.id$BIC, '\nRsquare = ', info.id$rsq)
legend <- "AICc = -330.10\nAIC = -331.01\nBIC = -321.55" #\nR^2 = 0.95
pdf(file = "../results/typ_gomp_fit.pdf")
p <- ggplot(typ.data, aes(x = Time, y = logN)) +
  geom_point() +
  labs(x = "Time (h)", y = "Log of Population Size (log(N))") +
  # ggtitle("Typical Gompertz Model Fit") +
  geom_line(data = plot.id, aes(x = time, y = plot.point), size=1) +
  theme_set(theme_bw()) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17)) +
#legend.position=c(3.5,0.39) 
  annotate('text', label = legend, x = 17, y = -2.3, hjust = 0, vjust = 0, size = 4.6)

print(p)
graphics.off()



###################################
# get info for Linear Mixed Model #
###################################
# library(lme4)
# lmm <- lmer(rmax ~ temp+(1|Species)+(1|Medium), data = info)
# ori.lm <- lm(rmax ~ temp, data = info)
# log.lm <- lm(log(rmax) ~ temp, data = info)
# summary(lmm)
# summary(ori.lm)
# summary(log.lm)

###################################
# r_max V.S. t_lag by temperature #
###################################
# 1. plot tlag VS rmax grouped by temp
pdf(file = "../results/rt_plot/log_rt_temp")
ggarrange(plotrt(df1), plotrt(df2), plotrt(df3), plotrt(df4),
          #labels = paste(df1$temp_group,df2$temp_group,df3$temp_group,df4$temp_group),
          ncol = 1, nrow = 4)
graphics.off()
# 2. plot tlag hist gourp by temperature
pdf(file = "../results/rt_plot/tlag_temp")
ggarrange(plott_temp(df1), plott_temp(df2), plott_temp(df3), plott_temp(df4),
          #labels = paste(df1$temp_group,df2$temp_group,df3$temp_group,df4$temp_group),
          ncol = 4, nrow = 1)
graphics.off()
# 3. plot rmax hist gourp by temperature

pdf(file = "../results/rt_plot/rmax_temp")
ggarrange(plotr_temp(df1), plotr_temp(df2), plotr_temp(df3), plotr_temp(df4),
          #labels = paste(df1$temp_group,df2$temp_group,df3$temp_group,df4$temp_group),
          ncol = 4, nrow = 1)
graphics.off()

# 4. plot k hist gourp by temperature

pdf(file = "../results/rt_plot/k_temp")
ggarrange(plotk_temp(df1), plotk_temp(df2), plotk_temp(df3), plotk_temp(df4),
          #labels = paste(df1$temp_group,df2$temp_group,df3$temp_group,df4$temp_group),
          ncol = 4, nrow = 1)
graphics.off()



##################################


ggplot(data = info, aes(x = log(1/tlag), y = log(rmax))) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm, se = TRUE)

ggplot(data = info, aes(x = log(rmax), y = log(1/tlag))) +
  geom_point() +
  stat_smooth(formula = y~x, method = lm, se = TRUE)

infos$rt <- infos$rmax*infos$tlag
ggplot(data = info, aes(x = log(1/tlag), y = log(rmax), colour = temp_group)) +
  geom_point() 
  stat_smooth(formula = y~x, method = lm, se = TRUE)





##################################
# r_max V.S. K by temperature#
##################################
# 1. plot K VS rmax grouped by temp
pdf(file = "../results/rt_plot/log_rk_temp")
ggarrange(plotrk(df1), plotrk(df2), plotrk(df3), plotrk(df4),
          #labels = paste(df1$temp_group,df2$temp_group,df3$temp_group,df4$temp_group),
          ncol = 1, nrow = 4)
graphics.off()


###############################
# linear mixed model analyses #
###############################


# ####################
# # r_max V.S. t_lag #
# ####################
# 
# # 1. plot tlag VS rmax grouped by temp
# info.plot <- subset(info, info$tlag > 0 & info$rmax > 0)
# png(filename = "../results/rt_plot/tr_plot")
# p <- ggplot(data = info.plot, aes(x = tlag, y = rmax, colour = temp_group)) +
#   geom_point(size = 1) +
#   geom_smooth(method='lm')
#   # scale_color_gradientn(colours = rainbow(33))
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/log_tr_plot")
# # yvar <- log(info.plot$rmax)
# # xvar <- log(info.plot$tlag)
# # fmla <-  yvar~ xvar
# p <- ggplot(data = info.plot, aes(x = log(tlag), y = log(rmax), colour = temp_group)) +
#   geom_point(size = 1) +
#   # stat_smooth(method = "lm", formula = log(rmax)~log(tlag), se = TRUE)
#   # stat_smooth(aes(fill= temp_group, color = temp_group), method='lm', se = TRUE, formula = fmla) +
#   # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~~')), formula = fmla)
#   geom_smooth(method='lm', se = TRUE)
#   # stat_poly_eq(formula = fmla, 
#   #              aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~')),
#   #              parse = TRUE) 
# 
# 
# 
# print(p)
# graphics.off()
# # library(ggpmisc)
# # library(ggpubr)
# 
# # 2. plot tlag hist gourp by temperature
# png(filename = "../results/rt_plot/tlag_hist")
# p <- ggplot(info.plot, aes(x = tlag, colour = temp_group, fill = temp_group)) +
#   geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
#   geom_density(alpha=0.6) +
#   geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
#              linetype="dashed") +
#   labs(title="tlag histogram plot",x="tlag", y = "count") +
#   theme_classic()
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/tlag_hist_facet")
# p <- ggplot(info.plot, aes(x = tlag, color = temp_group, fill = temp_group)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   # scale_fill_viridis(discrete = TRUE) +
#   # scale_color_viridis(discrete = TRUE) +
#   geom_density(alpha=0.6) +
#   geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
#              linetype="dashed") +
#   labs(title="tlag histogram plot",x="tlag", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   # xlab("tlag") +
#   # ylab("Count") +
#   facet_wrap(~temp_group)
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/log_tlag_hist")
# ggplot(info.plot, aes(x = log(tlag), colour = temp_group, fill = temp_group)) +
#   geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
#   geom_density(alpha=0.6) +
#   # geom_vline(data=info.plot, aes(xintercept=mean(log(info.plot$temp)), 
#                                  # color=temp_group), linetype="dashed") +
#   labs(title="log(tlag) histogram plot",x="tlag", y = "count") +
#   theme_classic()
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/log_tlag_hist_facet")
# p <- ggplot(info.plot, aes(x = log(tlag), color = temp_group, fill = temp_group)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   # scale_fill_viridis(discrete = TRUE) +
#   # scale_color_viridis(discrete = TRUE) +
#   geom_density(alpha=0.6) +
#   # geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
#   #            linetype="dashed") +
#   labs(title="tlag histogram plot",x="log(tlag)", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   # xlab("tlag") +
#   # ylab("Count") +
#   facet_wrap(~temp_group)
# print(p)
# graphics.off()
# 
# # 3. plot rmax hist gourp by temperature
# png(filename = "../results/rt_plot/rmax_hist")
# p <- ggplot(info.plot, aes(x = rmax, colour = temp_group, fill = temp_group)) +
#   geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
#   geom_density(alpha=0.6) +
#   geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
#              linetype="dashed") +
#   labs(title="rmax histogram plot",x="rmax", y = "count") +
#   theme_classic()
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/rmax_hist_facet")
# p <- ggplot(info.plot, aes(x = rmax, color = temp_group, fill = temp_group)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   geom_density(alpha=0.6) +
#   labs(title="rmax histogram plot",x="rmax", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   facet_wrap(~temp_group)
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/log_rmax_hist")
# p <- ggplot(info.plot, aes(x = log(rmax), colour = temp_group, fill = temp_group)) +
#   geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
#   geom_density(alpha=0.6) +
#   # geom_vline(data=info.plot, aes(xintercept=mean(log(info.plot$temp)), 
#   # color=temp_group), linetype="dashed") +
#   labs(title="log(rmax) histogram plot",x="rmax", y = "count") +
#   theme_classic()
# print(p)
# graphics.off()
# 
# png(filename = "../results/rt_plot/log_rmax_hist_facet")
# p <- ggplot(info.plot, aes(x = log(rmax), color = temp_group, fill = temp_group)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   geom_density(alpha=0.6) +
#   labs(title="rmax histogram plot",x="log(rmax)", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   facet_wrap(~temp_group)
# print(p)
# graphics.off()
# 
# 
# # plot by species
# png(filename = "../results/rt_plot/log_rmax_hist_facet_spe")
# p <- ggplot(info.plot, aes(x = log(rmax), color = Species, fill = Species)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   geom_density(alpha=0.6) +
#   labs(title="log(rmax) histogram plot group by species",x="log(rmax)", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   facet_wrap(~Species)
# print(p)
# graphics.off()
# 
# 
# png(filename = "../results/rt_plot/log_tlag_hist_facet_spe")
# p <- ggplot(info.plot, aes(x = log(tlag), color = Species, fill = Species)) +
#   geom_histogram(alpha=0.6, binwidth = 5) +
#   geom_density(alpha=0.6) +
#   labs(title="log(tlag) histogram plot group by species",x="log(tlag)", y = "count") +
#   theme_classic() +
#   theme(
#     legend.position = 'none',
#     panel.spacing = unit(0.1, 'lines'),
#     strip.text.x = element_text(size = 9)
#   ) +
#   facet_wrap(~Species)
# print(p)
# graphics.off()

################ check list
# color palette