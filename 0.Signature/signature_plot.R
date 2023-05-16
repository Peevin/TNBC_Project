library(ggridges)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggsignif)

setwd('/Users/liupeiwen/BC/New_analysis/0.Signature/')


adata <- LoadH5Seurat('../1.Toomanycells_result/Tcell/data/Tcell_PDL1.h5seurat')
# ann <- read.csv('./data/CD8_EX_SCC_obs.csv',row.names = 1)
# adata <- AddMetaData(adata,ann)
adata <- adata[,adata$Sub_Cluster=='t_CD4_Treg-FOXP3']
meta <- read.csv('../1.Toomanycells_result/Tcell/result/cluster/t_CD4_Treg-FOXP3_cluster.csv',row.names = 1)
adata <- AddMetaData(adata,meta)
# adata <- adata[,adata$cellSubType=='CD8_EX']
adata <- NormalizeData(adata,scale.factor = 1e6)

matrix <- as.data.frame(GetAssayData(adata, slot = "data"))
matrix <- as.data.frame(t(matrix))

b <- read.csv('./data/CD4_Treg-FOXP3_sig_weight.csv')
weight <- b$scale_w
names(weight) <- b$gene


matrix <- matrix[,colnames(matrix) %in% names(weight)]
adata$score <- 0
for (i in rownames(matrix)) {
  c <- matrix[i,]
  score <- 0
  for (j in colnames(matrix)) {
    score <- c[[j]] * weight[[j]] + score
  }
  adata$score[i] <- score
}




metabolism_score <- adata@meta.data %>% 
  group_by(cluster) %>% 
  summarize(sum_value = mean(score))


df <- adata@meta.data[c('cluster','score')]
colnames(df) <- c('cluster','score')
df <- df[order(df$cluster,decreasing = F),]
df$cluster <- as.character(df$cluster)
df$cluster <- paste0('Cluster',df$cluster)
df$cluster <- factor(df$cluster, levels = unique(df$cluster))
df$color <- ifelse(df$cluster %in% c('Cluster4','Cluster8','Cluster9',
                                     'Cluster11','Cluster16','Cluster18','Cluster19'), 
                   'NR-E','non-NR-E')


g1 <- ggplot(df,aes(cluster,score)) +
  stat_boxplot(geom="errorbar",width=0.5) +
  geom_boxplot(aes(fill=color),outlier.shape = NA) +

  scale_fill_manual(values = c('NR-E'="#4A7CB3","non-NR-E"="white")) +
  stat_summary(fun.y='mean',geom='point',shape=23,size=2,fill='white')+
  theme_classic() +
  labs(title = 'CD4_Treg-FOXP3') +
  ylab('NR-E score') +
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle = 90,size=15),
        axis.text.y = element_text(size=15),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        legend.position = 'none'
        ) +
  geom_hline(yintercept=mean(df$score),linetype=2)

g1


# compared <- list(c('16','15'),c('16','22'))
# 
# 
# g1 + 
#   geom_signif(comparisons = compared, 
#                  step_increase = 0.1,
#                  test =t.test ,map_signif_level = T)
# 
# 
# ggsave('figure/CD4_Treg-FOXP3_sig_score.pdf',width = 8,height = 5,dpi = 400)

# t.test(df[df$cluster=='5',]$score,df$score)
# t.test(df[df$cluster=='12',]$score,df$score)
# t.test(df[df$cluster=='13',]$score,df$score)$p.value
metabolism_score$cluster <- paste0('Cluster',metabolism_score$cluster)

for (i in 1:nrow(metabolism_score)) {
  metabolism_score$P.Val[i] <- t.test(df[df$cluster==metabolism_score$cluster[i],]$score,df$score,alternative = 'greater')$p.value

}
metabolism_score$Sig <- ifelse(metabolism_score$P.Val<0.05,'Sig','N.S.')
metabolism_score$adj.P <- p.adjust(metabolism_score$P.Val,method = 'BH')
metabolism_score$`-log(adj.P)` <- -log10(metabolism_score$adj.P)
mean(df$score)
metabolism_score$FC <- metabolism_score$sum_value / mean(df$score)
metabolism_score$log2FC <-log2(metabolism_score$FC)
metabolism_score <- metabolism_score[order(metabolism_score$FC,decreasing = T),]
write.csv(metabolism_score,'./result/CD4_Treg-FOXP3_signature_score.csv',quote = F,row.names = F)
#-------------------------------------------

