library(ggplot2)

setwd('/Users/liupeiwen/BC/New_analysis/1.Toomanycells_result/Tcell/')

cluster <- c('t_CD4-CXCL13', 't_CD4_Tcm-LMNA', 't_CD4_Tn-LEF1',
             't_CD4_Treg-FOXP3', 't_CD8-CXCL13', 't_CD8_MAIT-KLRB1',
             't_CD8_Teff-GNLY', 't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683',
             't_Tact-IFI6', 't_Tact-XIST', 't_Tprf-MKI67')

for (j in cluster) {
  adata <- read.csv(paste0('data/Cellcycle_count/',j,'_cellcycle_count.csv'))
  adata <- adata[order(adata$cluster),]
  
  df <- data.frame(row.names = unique(adata$cluster))
  df$cluster <- paste0('Cluster',rownames(df))
  
  for (i in 1:nrow(df)) {
    data <- adata[adata$cluster==rownames(df)[i],]
    df$cluster_count[i] <- nrow(data)
    df$G1_count[i] <- nrow(data[data$phase=='G1',])
    df$S_count[i] <- nrow(data[data$phase=='S',])
    df$G2M_count[i] <- nrow(data[data$phase=='G2M',])
    df$G1_Prop <- df$G1_count / df$cluster_count
    df$S_Prop <- df$S_count / df$cluster_count
    df$G2M_Prop <- df$G2M_count / df$cluster_count
  }
  
  df$cluster <- factor(df$cluster,levels = df$cluster)
  
  
  g1 <- ggplot(df,aes(x=cluster)) +
    geom_line(aes(y=S_Prop,group=1,color='lightgreen')) +
    geom_line(aes(y=G1_Prop,group=2,color='red')) +
    geom_line(aes(y=G2M_Prop,group=3,color='blue')) +
    theme_classic() +
    theme(axis.text.x = element_text (angle = 90, hjust = 1),plot.title = element_text (hjust = 0.5)) +
    ggtitle(j)
  
  ggsave(paste0('figure/Cellcycle/',j,'_cellcycle_line.pdf'),width = 10,height = 6,
         dpi=600)
}



adata <- read.csv('data/Cellcycle_count/t_CD8-CXCL13_cellcycle_count.csv')
adata <- adata[order(adata$cluster),]


df <- data.frame(row.names = unique(adata$cluster))
df$cluster <- paste('Cluster',rownames(df),sep='')


for (i in 1:nrow(df)) {
  data <- adata[adata$cluster==rownames(df)[i],]
  df$cluster_count[i] <- nrow(data)
  df$G1_count[i] <- nrow(data[data$phase=='G1',])
  df$S_count[i] <- nrow(data[data$phase=='S',])
  df$G2M_count[i] <- nrow(data[data$phase=='G2M',])
  df$G1_Prop <- df$G1_count / df$cluster_count
  df$S_Prop <- df$S_count / df$cluster_count
  df$G2M_Prop <- df$G2M_count / df$cluster_count

}


df$cluster <- factor(df$cluster,levels = df$cluster)

ggplot(df,aes(x=cluster)) +
  geom_line(aes(y=S_Prop,group=1,color='lightgreen')) +
  geom_line(aes(y=G1_Prop,group=2,color='red')) +
  geom_line(aes(y=G2M_Prop,group=3,color='blue')) +
  theme_classic() +
  theme(axis.text.x = element_text (angle = 90, hjust = 1),plot.title = element_text (hjust = 0.5)) +
  ggtitle('CD8-CXCL13')


ggplot(df,aes(x=cluster)) +
  # geom_col(aes(y=S_Prop,color='lightgreen')) +
  geom_col(aes(y=G1_Prop,color='red')) +
  theme_classic() +
  geom_col(aes(y=G2M_Prop,color='blue')) +
  theme(axis.text.x = element_text (angle = 90, hjust = 1),plot.title = element_text (hjust = 0.5))
  



