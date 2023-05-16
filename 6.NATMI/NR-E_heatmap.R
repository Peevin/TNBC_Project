library(readxl)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)

setwd('/Users/liupeiwen/BC/New_analysis/6.NATMI/')

files = list.files("data/lig_rec/", pattern = "xlsx", full.names = TRUE)
df = map_dfr(files, read_xlsx) 

cluster <- c('CD4-CXCL13', 'CD4_Tcm-LMNA', 'CD4_Tn-LEF1',
             'CD4_Treg-FOXP3','CD8-CXCL13',
             'CD8_Tem-GZMK', 'CD8_Trm-ZNF683','Bmem-CD27',
             'pB-IGHG1'
)
# sub <- c('CD4-CXCL13')

NR_E <- data.frame(row.names = unique(df$`sending cell`))

# Heatmap data
for (sub in cluster) {
  lig_nichenet <- read.csv(paste0('../5.NicheNet_rec_target/data/',sub,'_ligand_activity.csv'))
  lig_p <- lig_nichenet %>%
    filter(pearson>0)
  lig_rec_nich_p <- df %>%
    filter(ligand %in% lig_p$test_ligand)
  
  lig_rec_count <- data.frame(row.names = unique(df$`sending cell`))
  
  # df1 <- lig_rec_nich_p %>% filter(`target cell`==paste0('t_',sub))
  for (i in rownames(lig_rec_count)) {
    lig_rec_count[i,sub] <- length(unique(lig_rec_nich_p %>% filter(`sending cell`==i) %>% .$ligand))
  }
  
  NR_E <- cbind(NR_E, lig_rec_count)
}


# Heatmap(NR_E)
col <- colorRampPalette(c("white","#2878B5"))( 10)
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(NR_E[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 10
    ))
}

NR_E <- t(NR_E)

colnames(NR_E) <- substring(colnames(NR_E),3,nchar(colnames(NR_E)))

Heatmap(
  NR_E,
  # cluster_rows = F,
  # cluster_columns = F,
  show_column_dend = F,
  show_row_dend = F,
  name = 'ligand_count',
  col = col,
  row_names_side = 'left',
  # right_annotation = lig_ann,
  # top_annotation = rec_ann,
  cell_fun = cell_fun,
  # column_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2))),
  # row_names_gp = gpar(col=c('#32B897','black','black',rep('#C82423',12),rep('#32B897',2),rep('black',2),'#C82423','#32B897',rep('#2878B5',10),'black','#32B897',rep('#C82423',2),rep('#32B897',2),'#2878B5',rep('#C82423',2)))
)



