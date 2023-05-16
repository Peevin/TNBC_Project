library(readxl)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)

setwd('/Users/liupeiwen/BC/New_analysis/6.NATMI/')

files = list.files("data/lig_rec/", pattern = "xlsx", full.names = TRUE)
df = map_dfr(files, read_xlsx) 

lig_nichenet <- read.csv('../5.NicheNet_rec_target/data/CD8_Trm-ZNF683_ligand_activity.csv')
lig_p <- lig_nichenet %>%
  filter(pearson>0)

lig_rec_nich_p <- df %>%
  filter(ligand %in% lig_p$test_ligand)

lig_count <- lig_rec_nich_p %>% 
  distinct(`sending cell`,ligand) %>% 
  group_by(`sending cell`) %>% 
  summarise(n_distinct(ligand))
colnames(lig_count)<- c('type','lig_count')
lig_count <- data.frame(row.names = lig_count$type,lig_count=lig_count$lig_count)

# receptor count by cluster
rec_count <- lig_rec_nich_p %>% 
  distinct(`target cell`,receptor) %>% 
  group_by(`target cell`) %>% 
  summarise(n_distinct(receptor))
colnames(rec_count) <- c('type','rec_count')
rec_count <- data.frame(row.names = rec_count$type,rec_count=rec_count$rec_count)


lig_rec_count <- data.frame(row.names = rownames(lig_count))

for (i in rownames(lig_rec_count)) {
  df <- lig_rec_nich_p %>% filter(`sending cell`==i)
  for (j in rownames(rec_count)) {
    lig_rec_count[i,j] <- nrow(df %>% filter(`target cell`==j))
    
  }
}


col <- colorRampPalette(c("white", "#FDF5EB","#FAE7D2","#F5D3A8","#F2B176","#ED924E","#DF7031","#C85223","#993F18","#722D10"))( 10 )
# lig_col <- colorRamp2(breaks = seq(0,1500,length.out=10), colors = c("white", "#FDF5EB","#FAE7D2","#F5D3A8","#F2B176","#ED924E","#DF7031","#C85223","#993F18","#722D10"))


lig_ann <- rowAnnotation(
  lig_count=anno_barplot(lig_count$lig_count,
                         width = unit(2, "cm"),
                         gp = gpar('white'))
  
)

rec_ann <- columnAnnotation(
  rec_count=anno_barplot(rec_count$rec_count,
                         height = unit(2, "cm"),
                         gp = gpar('white'))
  
)

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(lig_rec_count[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 10
    ))
}


Heatmap(
  lig_rec_count,
  # cluster_rows = F,
  # cluster_columns = F,
  show_column_dend = F,
  show_row_dend = F,
  name = 'lig_rec_count',
  col = col,
  row_names_side = 'left',
  # right_annotation = lig_ann,
  # top_annotation = rec_ann,
  cell_fun = cell_fun,
  # column_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2))),
  # row_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2)))
)

Heatmap(
  lig_rec_count,
  # cluster_rows = F,
  # cluster_columns = F,
  show_column_dend = F,
  show_row_dend = F,
  col = col,
  row_names_side = 'left',
  # right_annotation = lig_ann,
  name = 'lig_rec_count',
  # top_annotation = rec_ann,
  cell_fun = cell_fun,
  # column_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2))),
  # row_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2)))
)


a <- lig_rec_nich_p %>% filter(`target cell`=='t_CD8-CXCL13')
b <- a[grepl('macro',a$`sending cell`),]
b$`ligand mean expression` <- log10(b$`ligand mean expression`)
b$`ligand mean expression` <- scale(b$`ligand mean expression`)
b$`receptor mean expression` <- log10(b$`receptor mean expression`)
b$`receptor mean expression` <- scale(b$`receptor mean expression`)

ggplot(b,aes(x=`sending cell`, y=ligand,color=`ligand mean expression`)) +
  geom_point() +
  scale_color_gradient2(low = "blue", # 下限颜色
                        mid = "white", # 中值颜色
                        high = "red3", # 上限颜色
                        midpoint = 0
  ) +
  theme_classic()

ggplot(b,aes(x=`sending cell`, y=ligand,color=`ligand mean expression`)) +
  geom_point() +
  scale_color_gradient2(low = "blue", # 下限颜色
                        mid = "white", # 中值颜色
                        high = "red3", # 上限颜色
                        midpoint = 0
                        ) +
  theme_bw()+
  geom_point(aes(# size=`value`,
    color=`ligand mean expression`))+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  geom_vline(xintercept=c(-3,3.5,8.5,11.5,15.5),size=.8)

p1<-ggplot(b,aes(x=`target cell`, y=receptor)) #热图绘制
p2 <- p1+
  scale_color_gradient2(low = "blue", # 下限颜色
                        mid = "white", # 中值颜色
                        high = "red3", # 上限颜色
                        midpoint = 0
  )+
  # scale_color_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))+
  theme_bw()+
  geom_point(aes(# size=`value`,
                 color=`receptor mean expression`))+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  geom_vline(xintercept=c(-3,3.5,8.5,11.5,15.5),size=.8)
p2
