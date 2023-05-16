library(ggplot2)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

setwd('/Users/liupeiwen/BC/New_analysis/1.Toomanycells_result/Tcell/')

adata <- LoadH5Seurat('data/Tcell_PDL1.h5seurat')
adata <- subset(adata,Sub_Cluster=='t_CD8-CXCL13')
meta <- read.csv('result/cluster/t_CD8-CXCL13_cluster.csv',row.names = 1)
adata <- AddMetaData(adata,meta)



adata$cluster <- paste0('cluster',adata$cluster)
adata$cluster <- factor(adata$cluster,levels = unique(adata$cluster))

df <- adata@meta.data
df <- df[order(df$cluster),]
df$cluster <- paste0('Cluster',df$cluster)
df$cluster <- factor(df$cluster,levels = unique(df$cluster))


ggplot(df,aes(x=cluster,fill=Sample)) +
  geom_bar(position = 'fill') +
  theme_classic() +
  scale_fill_manual(values = c("Post_P002_T"="#D1352B","Post_P005_T"="#826F9D",
                               "Post_P012_T"="#599592","Post_P017_T"="#70A268",
                               "Post_P019_T"="#8B5F99","Pre_P002_T"="#D2766A",
                               "Pre_P004_T"="#F5BE48","Pre_P005_T"="#F8F65E",
                               "Pre_P007_T"="#B38D3E","Pre_P012_T"="#CA7592",
                               "Pre_P016_T"="#D68CB4","Pre_P019_T"="#999999")
                    ) +
  theme(axis.text.x = element_text(angle = 90,size=10),
        axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        legend.position = 'none'
        ) +
  ylab("Fraction")

ggsave('figure/CD8-CXCL13_sample_fraction.pdf',width = 8,height = 6,dpi=600)


ori <- nrow(df[df$Response=='SD',]) / nrow(df)
ggplot(df,aes(cluster,fill=Response)) +
  geom_bar(position = 'fill') +
  theme_classic() +
  scale_fill_manual(values = c("PR"="#D1352B","SD"="#4A7CB3")
                    ) +
  theme(axis.text.x = element_text(angle = 90,size=10),
        axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        legend.position = 'none'
  ) +
  geom_hline(yintercept=ori,linetype=2, size=0.2) +
  ylab("Fraction")
ggsave('figure/CD8-CXCL13_Response_fraction.pdf',width = 8,height = 6,dpi=600)
