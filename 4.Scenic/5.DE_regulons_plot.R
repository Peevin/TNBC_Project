# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk",force = TRUE)


library(SeuratDisk)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)


library(tidyverse)
library(magrittr)
library(viridis)

library(network)
library(igraph)

cluster <- 't_CD8-CXCL13'
TMC <- c(12,13,16)
setwd('/Users/liupeiwen/BC/New_analysis/4.Scenic/')
#读取h5ad数据。h5ad是python的Scanpy读取文件格式，对其进行格式转换，并得到.h5seurat格式的文件
# Convert('./data/t_CD8-CXCL13.h5ad', "h5seurat",
#         overwrite = TRUE,assay = "RNA")
adata <- LoadH5Seurat("./data/Tcell_PDL1.h5seurat")
adata <- subset(adata,Sub_Cluster==cluster)
ann <- read.csv(paste0('../1.Toomanycells_result/Tcell/result/cluster/',cluster,'_cluster.csv'),
                row.names = 1)
adata <- AddMetaData(adata,ann)

adata$divide <- ifelse(adata$cluster %in% TMC,'NR_Enrich','non-NR_Enrich')

Idents(adata) <- 'divide'
Idents(adata)

# 读取scenic数据
scenicLoomPath <- file.path(paste0('./data/',cluster,'/pyscenic_output.loom'))
motifEnrichmentFile <- file.path(paste0('./data/',cluster,'/reg.csv'))

file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)

loom <- open_loom(scenicLoomPath)
# Read information from loom file:
# exprMat <- get_dgem(loom)
# exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
# embeddings <- get_embeddings(loom)
close_loom(loom)

AUCmat <- AUCell::getAUC(regulonAUC)

adj_regulonAucTresholds <- as.data.frame(regulonAucThresholds)
adj_regulonAucTresholds$regulonAucThreshold <- rownames((adj_regulonAucTresholds))
rownames(adj_regulonAucTresholds) <- adj_regulonAucTresholds$regulonAucThresholds
regulonAucThresholds <- as.character(adj_regulonAucTresholds[, 2])
names(regulonAucThresholds) <- rownames(adj_regulonAucTresholds)
regulonAucThresholds


binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}

binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)

adata[['AUC']] <- CreateAssayObject(data = AUCmat)
adata[['AUCBinary']] <- CreateAssayObject(data = binaryRegulonActivity)


#  Identify differential regulons based AUC scores
DefaultAssay(adata) <- 'AUC'

deg <- FindMarkers(adata, ident.1 = 'NR_Enrich', ident.2 = 'non-NR_Enrich', 
                   logfc.threshold = 0, test.use = 't')

deg.up <- deg[which(deg$avg_log2FC > 0 & deg$p_val_adj < 0.01), ]
deg.up$TF <- rownames(deg.up)
deg.dn <- deg[which(deg$avg_log2FC < -0 & deg$p_val_adj < 0.01), ]
deg.ls <- list(deg.up, deg.dn)
names(deg.ls) <- c('NR_Enrich','non-NR_Enrich')

# Visu check differential regulons
adata <- ScaleData(adata, assay = 'AUC', features = rownames(AUCmat))

c_deup <- as.vector(rownames(deg.up))
c_dedn <- as.vector(rownames(deg.dn))

c_de <- c(c_deup,c_dedn)


Idents(adata) <- factor(Idents(adata),levels = c('NR_Enrich','non-NR_Enrich'))

g1 <- DoHeatmap(adata, features = c_deup, slot = 'scale.data', 
                raster = F,group.colors = c('#4A7CB3','#D1352B'))+ 
  scale_fill_gradient2(
      low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
      mid = "white",
      high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
      midpoint = 0,
      guide = "colourbar",
      aesthetics = "fill"
    )
g1
g2 <- DoHeatmap(adata, features = c_dedn, slot = 'scale.data', raster = F)+ scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 0,
  guide = "colourbar",
  aesthetics = "fill"
)
g2

g3 <- DoHeatmap(adata, features = c_de, slot = 'scale.data', raster = F)+ scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 0,
  guide = "colourbar",
  aesthetics = "fill"
)
g3

ggsave('figure/t_CD8-CXCL13_DE_regulons_up.pdf',g1,width = 10,height = 8, dpi=1000)







meta <- read.csv('./output/t_CD8-CXCL13_DE_regulons_up_sig_phyper_test.csv',row.names = 1)
meta_sig <- deg.up[deg.up$TF %in% meta$regulon_name,]
meta_sig <- meta_sig[meta_sig$avg_log2FC>0,]
c_deup_sig <- as.vector(rownames(meta_sig))
g3 <- DoHeatmap(adata, features = c_deup_sig, slot = 'scale.data', raster = F,
                group.colors=c('blue','red'))+ 
  scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 0,
  guide = "colourbar",
  aesthetics = "fill"
)
g3
ggsave('figure/t_CD8-CXCL13_DE_regulons_up_sig.pdf',g3,width = 10,height = 8, dpi=700)



# g2 <- DoHeatmap(adata, features = c_dedn, slot = 'scale.data', raster = F)+ scale_fill_gradient2(
#       low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
#       mid = "white",
#       high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
#       midpoint = 0,
#       guide = "colourbar",
#       aesthetics = "fill"
#     )
# g2
# ggsave('figure/DE_regulons_down.png',g2,width = 10,height = 8, dpi=700)

# DoHeatmap(adata, features = c_combine, slot = 'scale.data', raster = F) + scale_fill_viridis()
# DoHeatmap(adata, features = c_combine, slot = 'scale.data', raster = F) + scale_fill_gradient2(
#   low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
#   mid = "white",
#   high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
#   midpoint = 0,
#   guide = "colourbar",
#   aesthetics = "fill"
# )

# 输出差异regulons文件
deg_out <- regulons[names(regulons) %in% c_deup_sig]
##这里snps是一个大list
c <- file("./output/t_CD8-CXCL13_up_sig_regulons.txt", "wa")
for(i in 1:length(deg_out)){
  o <- as.character(c(names(deg_out[i]), deg_out[[i]]))
  writeLines(o,  c, sep = "\t")
  writeLines("\n", c, sep = "" )
}
close(c)

deg_fc <- deg[rownames(deg) %in% c_deup,c(2,5)]
deg_fc <- deg_fc[order(deg_fc$avg_log2FC,decreasing = T),]
write.table(deg_fc, './output/t_CD8-CXCL13_DE_up_regulons_fc_qval.csv',row.names=T,quote = F)
