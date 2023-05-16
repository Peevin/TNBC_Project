library(ggplot2)
library(ggvenn)
library(dplyr)

setwd('/Users/liupeiwen/BC/New_analysis/0.Signature/')

adata <- readRDS('./data/CD8-CXCL13_sig.rds')
names(adata) <- paste0('Cluster',names(adata))
ggvenn(adata,show_percentage = F,
       fill_color = c('#63b2ee','#ffa510','#f89588'),
       fill_alpha = 0.7,
       text_size = 8
       )
  





df <- data.frame(cluster=names(adata))


for (j in 1:nrow(df)) {
  df$count[j] <- nrow(adata[[df$cluster[j]]])
}

df$cluster <- factor(df$cluster,levels = df$cluster)

ggplot(df,aes(x=cluster,y=count)) +
  geom_col(color='black',fill='white') +
  theme_classic() +
  ylab('signature count')


