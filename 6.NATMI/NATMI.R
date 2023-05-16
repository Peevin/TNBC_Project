library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(openxlsx)

setwd('/Users/liupeiwen/BC/New_analysis/6.NATMI/')

lig_rec_count <- read.xlsx('./data/c2c_source_matrix.xlsx',rowNames = T)


lig_rec <- lig_rec %>% 
  filter(Ligand.detection.rate>=0.2) %>% 
  filter(Receptor.detection.rate>=0.2)

lig_nichenet <- read.csv('./data/ligand_activity.csv')
lig_p <- lig_nichenet %>%
  filter(pearson>0)
# lig_p <- lig_nichenet %>%
#   filter(pearson<0)

lig_rec_nich_p <- lig_rec

lig_rec_nich_p <- lig_rec %>%
  filter(Ligand.symbol %in% lig_p$test_ligand)

# ligand count by cluster
lig_count <- lig_rec_nich_p %>% 
  distinct(Sending.cluster,Ligand.symbol) %>% 
  group_by(Sending.cluster) %>% 
  summarise(n_distinct(Ligand.symbol))
colnames(lig_count)<- c('type','lig_count')
lig_count <- data.frame(row.names = lig_count$type,lig_count=lig_count$lig_count)

# receptor count by cluster
rec_count <- lig_rec_nich_p %>% 
  distinct(Target.cluster,Receptor.symbol) %>% 
  group_by(Target.cluster) %>% 
  summarise(n_distinct(Receptor.symbol))
colnames(rec_count) <- c('type','rec_count')
rec_count <- data.frame(row.names = rec_count$type,rec_count=rec_count$rec_count)

aggregate(lig_rec_nich_p$Target.cluster,by=list(type=lig_rec_nich_p$Target.cluster),length)

lig_rec_count <- data.frame(row.names = rownames(lig_count))

for (i in rownames(lig_rec_count)) {
  df <- lig_rec_nich_p %>% filter(Sending.cluster==i)
  for (j in rownames(rec_count)) {
    lig_rec_count[i,j] <- nrow(df %>% filter(Target.cluster==j))
    
  }
}


col <- colorRampPalette(c("white", "#FDF5EB","#FAE7D2","#F5D3A8","#F2B176","#ED924E","#DF7031","#C85223","#993F18","#722D10"))( 10 )
lig_col <- colorRamp2(breaks = seq(0,1500,length.out=10), colors = c("white", "#FDF5EB","#FAE7D2","#F5D3A8","#F2B176","#ED924E","#DF7031","#C85223","#993F18","#722D10"))


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
  right_annotation = lig_ann,
  top_annotation = rec_ann,
  cell_fun = cell_fun,
  column_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2))),
  row_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2)))
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


a <- lig_rec_nich_p %>% filter(Target.cluster=='CD8_EX Tcell')

