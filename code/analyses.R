
rm(list = ls())
graphics.off()

setwd("~/Documents/CMEEProject/code")

require(ggplot2)

# 1. read the data, starting value and compare models
infos <- read.csv("../data/gomp.info.csv")
plot.df <- read.csv("../data/gomp.plot.csv")

Data <- read.csv('../data/mini_pop.csv')
Data <- Data[order(Data[,'id'], Data[,'Time']),]

IDs <- unique(Data$id)

################
# plot fitting #
################
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
  # rbind the data frame
  info <- rbind(info, info.id)
  
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
      # ggtitle("Model comparison plot") +
      geom_line(data = plot.id, aes(x = time, y = plot.point), size=1) +
      # theme(legend.position = 'bottom') +
      annotate('text', label = legend, x = min(data$Time), y = min(data$logN), hjust = -.5, vjust = 0) 
    # stat_smooth(method = lm, level = 0.95, aes(colour="Cubic")) +
    # scale_colour_manual(name="Model", values=c("darkblue", "darkred", "darkgreen"))
    print(p)
    graphics.off()
  }
}

# devide the temperature into 4 groups ((-5:10, 10:20, 20:30, 30:40))
df1 <- subset(info, info$temp >= -5 & info$temp < 10)
df2 <- subset(info, info$temp >= 10 & info$temp < 20)
df3 <- subset(info, info$temp >= 20 & info$temp < 30)
df4 <- subset(info, info$temp >= 30 & info$temp < 40)
df1$temp_group <- "-5:10"
df2$temp_group <- "10:20"
df3$temp_group <- "20:30"
df4$temp_group <- "30:40"
info <- rbind(df1,df2,df3,df4)

####################
# r_max V.S. t_lag #
####################

# 1. plot tlag VS rmax grouped by temp
info.plot <- subset(info, info$tlag > 0 & info$rmax > 0)
png(filename = "../results/rt_plot/tr_plot")
p <- ggplot(data = info.plot, aes(x = tlag, y = rmax, colour = temp_group)) +
  geom_point(size = 1) +
  geom_smooth(method='lm')
  # scale_color_gradientn(colours = rainbow(33))
print(p)
graphics.off()

png(filename = "../results/rt_plot/log_tr_plot")
# yvar <- log(info.plot$rmax)
# xvar <- log(info.plot$tlag)
# fmla <-  yvar~ xvar
p <- ggplot(data = info.plot, aes(x = log(tlag), y = log(rmax), colour = temp_group)) +
  geom_point(size = 1) +
  # stat_smooth(method = "lm", formula = log(rmax)~log(tlag), se = TRUE)
  # stat_smooth(aes(fill= temp_group, color = temp_group), method='lm', se = TRUE, formula = fmla) +
  # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~~')), formula = fmla)
  geom_smooth(method='lm', se = TRUE)
  # stat_poly_eq(formula = fmla, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~')),
  #              parse = TRUE) 



print(p)
graphics.off()
# library(ggpmisc)
# library(ggpubr)

# 2. plot tlag hist gourp by temperature
png(filename = "../results/rt_plot/tlag_hist")
p <- ggplot(info.plot, aes(x = tlag, colour = temp_group, fill = temp_group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6) +
  geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
             linetype="dashed") +
  labs(title="tlag histogram plot",x="tlag", y = "count") +
  theme_classic()
print(p)
graphics.off()

png(filename = "../results/rt_plot/tlag_hist_facet")
p <- ggplot(info.plot, aes(x = tlag, color = temp_group, fill = temp_group)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  # scale_fill_viridis(discrete = TRUE) +
  # scale_color_viridis(discrete = TRUE) +
  geom_density(alpha=0.6) +
  geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
             linetype="dashed") +
  labs(title="tlag histogram plot",x="tlag", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  # xlab("tlag") +
  # ylab("Count") +
  facet_wrap(~temp_group)
print(p)
graphics.off()

png(filename = "../results/rt_plot/log_tlag_hist")
ggplot(info.plot, aes(x = log(tlag), colour = temp_group, fill = temp_group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
  geom_density(alpha=0.6) +
  # geom_vline(data=info.plot, aes(xintercept=mean(log(info.plot$temp)), 
                                 # color=temp_group), linetype="dashed") +
  labs(title="log(tlag) histogram plot",x="tlag", y = "count") +
  theme_classic()
print(p)
graphics.off()

png(filename = "../results/rt_plot/log_tlag_hist_facet")
p <- ggplot(info.plot, aes(x = log(tlag), color = temp_group, fill = temp_group)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  # scale_fill_viridis(discrete = TRUE) +
  # scale_color_viridis(discrete = TRUE) +
  geom_density(alpha=0.6) +
  # geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
  #            linetype="dashed") +
  labs(title="tlag histogram plot",x="log(tlag)", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  # xlab("tlag") +
  # ylab("Count") +
  facet_wrap(~temp_group)
print(p)
graphics.off()

# 3. plot rmax hist gourp by temperature
png(filename = "../results/rt_plot/rmax_hist")
p <- ggplot(info.plot, aes(x = rmax, colour = temp_group, fill = temp_group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6) +
  geom_vline(data=info.plot, aes(xintercept=mean(info.plot$temp), color=temp_group),
             linetype="dashed") +
  labs(title="rmax histogram plot",x="rmax", y = "count") +
  theme_classic()
print(p)
graphics.off()

png(filename = "../results/rt_plot/rmax_hist_facet")
p <- ggplot(info.plot, aes(x = rmax, color = temp_group, fill = temp_group)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  geom_density(alpha=0.6) +
  labs(title="rmax histogram plot",x="rmax", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~temp_group)
print(p)
graphics.off()

png(filename = "../results/rt_plot/log_rmax_hist")
p <- ggplot(info.plot, aes(x = log(rmax), colour = temp_group, fill = temp_group)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
  geom_density(alpha=0.6) +
  # geom_vline(data=info.plot, aes(xintercept=mean(log(info.plot$temp)), 
  # color=temp_group), linetype="dashed") +
  labs(title="log(rmax) histogram plot",x="rmax", y = "count") +
  theme_classic()
print(p)
graphics.off()

png(filename = "../results/rt_plot/log_rmax_hist_facet")
p <- ggplot(info.plot, aes(x = log(rmax), color = temp_group, fill = temp_group)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  geom_density(alpha=0.6) +
  labs(title="rmax histogram plot",x="log(rmax)", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~temp_group)
print(p)
graphics.off()


# plot by species
png(filename = "../results/rt_plot/log_rmax_hist_facet_spe")
p <- ggplot(info.plot, aes(x = log(rmax), color = Species, fill = Species)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  geom_density(alpha=0.6) +
  labs(title="log(rmax) histogram plot group by species",x="log(rmax)", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~Species)
print(p)
graphics.off()


png(filename = "../results/rt_plot/log_tlag_hist_facet_spe")
p <- ggplot(info.plot, aes(x = log(tlag), color = Species, fill = Species)) +
  geom_histogram(alpha=0.6, binwidth = 5) +
  geom_density(alpha=0.6) +
  labs(title="log(tlag) histogram plot group by species",x="log(tlag)", y = "count") +
  theme_classic() +
  theme(
    legend.position = 'none',
    panel.spacing = unit(0.1, 'lines'),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~Species)
print(p)
graphics.off()

################ check list
# color palette