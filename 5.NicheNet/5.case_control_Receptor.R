library(ggplot2)
library(dplyr)
library(reshape2)

setwd('/Users/liupeiwen/BC/New_analysis/CD8-CXCL13/5.NicheNet/')

i = 'TBX21'
data <- read.csv(paste('./result/',i,'_top10_Receptor_mean_expr_regscore.csv',sep = ''))


ggplot(data, aes(x=group, y=reorder(receptor,regulatory_score))) +
  geom_point(aes(size=detection_rate, fill=mean_expression),
             shape=21) +
  scale_fill_distiller(palette = "Reds",direction = 1) +
  theme_classic() +
  ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab('Receptors')




