library(openxlsx)
library(dplyr)
library(ComplexHeatmap)

setwd('/Users/liupeiwen/BC/New_analysis/6.NATMI/')

lig_rec_count <- read.xlsx('./data/c2c_source_matrix_NATMI_NicheNet_subcluster.xlsx',rowNames = T)




col <- colorRampPalette(c("white", "#FDF5EB","#FAE7D2","#F5D3A8","#F2B176","#ED924E","#DF7031","#C85223","#993F18","#722D10"))(10)
# col <- colorRampPalette(c("white", "#FDF5EB","#F5D3A8","#ED924E","#DF7031","#993F18"))(6)

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
  column_title = 'Cell-type expressing receptor',
  column_title_side = 'bottom',
  row_title = 'Cell-type expressing ligand',
  # column_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2))),
  # row_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2)))
)

