library(reshape2)
library(ggplot2)
tasks = c("Literature Review", "Refining\nProject Design","Preliminary Simulations", "Model Fitting &\nModel Selection", "Thesis Writing")
tasks = factor(tasks, levels = tasks)
start.date = as.Date(c("2021-01-10", "2021-01-10", "2021-03-01", "2021-03-15", "2021-03-15"))
end.date = as.Date(c("2021-02-10", "2021-03-15", "2021-04-01", "2021-08-01", "2021-09-15"))
chart = data.frame(tasks, start.date,end.date)
chart = melt(chart, measure.vars = c("start.date", "end.date"))
#pdf("Gantt_chart.pdf", width = 9, height = 3)

ggplot(chart, aes(value, tasks)) + 
  geom_line(size = 6) + theme_bw()+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 months") +
  xlab(NULL) +
  ylab(NULL)
# graphics.off()
