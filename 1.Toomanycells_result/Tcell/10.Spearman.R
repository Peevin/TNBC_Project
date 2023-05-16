library(ggplot2)
library(Hmisc)
library(tidyverse)
library(reshape2)
library(ComplexHeatmap)

setwd('/Users/liupeiwen/BC/New_analysis/1.Toomanycells_result/Tcell/')

cluster <- c('t_CD4-CXCL13', 't_CD4_Tcm-LMNA', 't_CD4_Tn-LEF1',
             't_CD4_Treg-FOXP3', 't_CD8-CXCL13', 't_CD8_MAIT-KLRB1',
             't_CD8_Teff-GNLY', 't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683',
             't_Tact-IFI6', 't_Tact-XIST', 't_Tprf-MKI67')

df <- data.frame(row.names = cluster)
df$Cluster <- rownames(df)

for (j in cluster) {
  adata <- read.csv(paste0('data/spearman_data/',j,'_spearman_data.csv'),)
  
  adata$Cluster = paste('Cluster',adata$Cluster,sep = '')
  adata$Cluster <- factor(adata$Cluster,levels = as.vector(adata$Cluster))
  result1 <- cor.test(adata$G1_Prop,adata$exhuast_score,method = 'spearman')
  result2 <- cor.test(adata$S_Prop,adata$exhuast_score,method = 'spearman')
  result3 <- cor.test(adata$G2M_Prop,adata$exhuast_score,method = 'spearman')
  
  df[j,'G1_spearman_cor'] <- result1$estimate
  df[j,'S_spearman_cor'] <- result2$estimate
  df[j,'G2M_spearman_cor'] <- result3$estimate
  
  df[j,'G1_p.Val'] <- result1$p.value
  df[j,'S_p.Val'] <- result2$p.value
  df[j,'G2M_p.Val'] <- result3$p.value
  
}

df <- df[-(nrow(df)),]
data <- as.data.frame(df[,2:4]) %>% #将矩阵转换成数据框
  mutate(x=rownames(df)) %>%  #新建一列x，是11种属性变量
  melt(id='x') %>%                   #将宽数据转换成长数据，更适合ggplot2绘图
  rename('y'='variable','Corr'='value')  #将variable列名改成y



ggplot(data,aes(y,x,fill=Corr)) +
  geom_tile()+  #色块函数
  scale_fill_gradient2(low = 'Blue',mid = 'white',high ='red',
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL)+
  theme_()


df_heat <- df[,2:4]

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(df_heat[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 10
    ))
}

Heatmap(df_heat,
        border = "black",
        cluster_columns =F,
        cluster_rows = F,
        row_names_side = 'left',
        column_title = 'exhuast-cellcycle spearman corr',
        column_title_side = 'top', 
        column_title_gp = gpar(fontsize = 15,
                               fontface = 'bold', 
                               # fill = 'red',
                               # col = 'white',
                               border = F
                              ),
        name = 'Spearman_corr',
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = cell_fun
        )




# ----------------code test-----------------
adata <- read.csv('data/spearman_data/t_CD4_Tcm-LMNA_spearman_data.csv')

adata$Cluster = paste('Cluster',adata$Cluster,sep = '')
adata$Cluster <- factor(adata$Cluster,levels = as.vector(adata$Cluster))

cor(adata$S_Prop,adata$exhuast_score,method = 'spearman')
cor.test(adata$S_Prop,adata$exhuast_score,method = 'spearman')

cor(adata$G1_Prop,adata$exhuast_score,method = 'spearman')
cor.test(adata$G1_Prop,adata$exhuast_score,method = 'spearman')

cor(adata$G2M_Prop,adata$exhuast_score,method = 'spearman')
a <- cor.test(adata$G2M_Prop,adata$exhuast_score,method = 'spearman')
a$p.value
a$estimate
