library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(ggsignif)

setwd('/Users/liupeiwen/BC/New_analysis/1.Toomanycells_result/Tcell/')

data <- read.csv('data/lineplot/t_CD8-CXCL13lineplot.csv',row.names = 1)

subdata <- filter(data, data$Patient %in% c('P019','P012','P002','P005'))

ggplot(subdata, aes(x=Timeline, y=fraction, color=Response, group=Response))+
  scale_x_discrete(limits = c('Pre_treatment','Post_treatment'))+
  theme_classic()+
  ylab('propotion of \n non-responde related cluste \n in CD8-CXCL13 cells (%)')+
  geom_point(size=3) +
  
  geom_path(group=subdata$Patient)

ggplot(data,aes(x=Response,y=fraction)) +
  geom_boxplot() +
  geom_jitter(width=0.4) +
  theme_classic()

