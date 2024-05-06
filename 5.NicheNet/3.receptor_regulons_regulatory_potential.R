library(dplyr)
library(ggplot2)
library(nichenetr)

setwd('/Users/liupeiwen/BC/New_analysis/CD8-CXCL13/5.NicheNet/')

# CD8-CXCL13
# read regulons
CD8_CXCL13_regulons <- readLines('../4.Scenic/output/DE_regulons.txt')
res <- strsplit(CD8_CXCL13_regulons, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -1)
CD8_CXCL13_regulons <- res
CD8_CXCL13_regulons <- CD8_CXCL13_regulons[c('JUND(+)','RARA(+)','E2F3(+)',
                                             'GABPA(+)','ZFX(+)','FOXJ2(+)',
                                             'MTF2(+)','SMAD3(+)','TBX21(+)','MAFF(+)')]

# read up-regulated metabolic genes
CD8_CXCL13_metgene <- read.csv('../4.Scenic//data/t_CD8-CXCL13_nebula_result.csv')
CD8_CXCL13_metgene <- CD8_CXCL13_metgene[CD8_CXCL13_metgene$expression=='Up-regulated',]
CD8_CXCL13_metgene <- as.character(CD8_CXCL13_metgene$gene)

# consturct metabolic regulons
CD8_CXCL13_met_regulons <- list()
for (i in names(CD8_CXCL13_regulons)) {
  CD8_CXCL13_met_regulons[[i]] <- intersect(CD8_CXCL13_regulons[[i]], CD8_CXCL13_metgene)
}

# macro-CCL2
# read regulons
macro_CCL2_regulons <- readLines('macro-CCL2/Scenic/output/DE_regulons.txt')
res <- strsplit(macro_CCL2_regulons, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -1)
macro_CCL2_regulons <- res
macro_CCL2_regulons <- macro_CCL2_regulons['CREB1(+)']

# read up-regulated metabolic genes
macro_CCL2_metgene <- read.csv('macro-CCL2/Scenic/data/t_macro-CCL2_nebula_result.csv')
macro_CCL2_metgene <- macro_CCL2_metgene[macro_CCL2_metgene$expression=='Up-regulated',]
macro_CCL2_metgene <- as.character(macro_CCL2_metgene$gene)

# consturct metabolic regulons
macro_CCL2_met_regulons <- list()
for (i in names(macro_CCL2_regulons)) {
  macro_CCL2_met_regulons[[i]] <- intersect(macro_CCL2_regulons[[i]], macro_CCL2_metgene)
}


# mDC-LAMP3
# read regulons
mDC_LAMP3_regulons <- readLines('mDC-LAMP3/Scenic/output/DE_regulons.txt')
res <- strsplit(mDC_LAMP3_regulons, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -1)
mDC_LAMP3_regulons <- res
mDC_LAMP3_regulons <- mDC_LAMP3_regulons[c('NFKB1(+)','FOS(+)',
                                           'JUN(+)','CREM(+)','YY1(+)',
                                           'KLF6(+)','SREBF2(+)','JUNB(+)',
                                           'RUNX3(+)')]

# read up-regulated metabolic genes
mDC_LAMP3_metgene <- read.csv('mDC-LAMP3/Scenic/data/t_mDC-LAMP3_nebula_result.csv')
mDC_LAMP3_metgene <- mDC_LAMP3_metgene[mDC_LAMP3_metgene$expression=='Up-regulated',]
mDC_LAMP3_metgene <- as.character(mDC_LAMP3_metgene$gene)

# consturct metabolic regulons
mDC_LAMP3_met_regulons <- list()
for (i in names(mDC_LAMP3_regulons)) {
  mDC_LAMP3_met_regulons[[i]] <- intersect(mDC_LAMP3_regulons[[i]], mDC_LAMP3_metgene)
}

# bind 3 metabolic regulons
met_regulons <- c(CD8_CXCL13_met_regulons,macro_CCL2_met_regulons,mDC_LAMP3_met_regulons)

met_regulons <- CD8_CXCL13_met_regulons

setwd('/Users/liupeiwen/NicheNet/')

Receptor_target <- readRDS('OmniNetworks_NNformat/receptor_TF_matrixWithweights.rds')

TF <- read.table('OmniNetworks_NNformat/allTFs_hg38.txt',col.names = 'TF')


Receptor_TF <- Receptor_target[rownames(Receptor_target) %in% TF$TF,]

up_reg <- c('JUND','RARA','E2F3',
            'GABPA','ZFX','FOXJ2',
            'MTF2','SMAD3','TBX21','MAFF')

Receptor_TF_sub <- as.data.frame(Receptor_TF[rownames(Receptor_TF) %in% up_reg,])

heat_plot <- list()

getColNumber<-function(data,colname){
  CN <- colnames(data)
  for (j in 1:ncol(data)){
    if (CN[j] == colname) {
      return (j)
    }
  }
  return (FALSE)
}

moveColToTop<-function(data,colname){
  colNumber <- getColNumber(data,colname)
  data <- data[,c(colNumber,2:(colNumber-1),(colNumber+1):ncol(data))]
  
  return (data)
}

for (i in up_reg) {
  TF <- c(i,met_regulons[[paste(i,'(+)',sep = '')]])
  TF_sub <- Receptor_TF_sub[i,]
  TF_sub <- t(TF_sub)
  TF_sub <- as.data.frame(TF_sub)
  TF_sub$receptor <- rownames(TF_sub)
  TF_sub <- TF_sub[order(TF_sub[,i],decreasing = T),]
  TF_rec <- head(TF_sub,10)
  colnames(TF_rec)[1] <- 'regulatory_score'
  write.csv(TF_rec,paste('/Users/liupeiwen/BC/New_analysis/CD8-CXCL13/5.NicheNet/result//',i,'_top10_Receptor_mean_expr_regscore.csv',sep = ''),
            quote = F, row.names = F)
  TF_rec <- TF_rec$receptor
  
  TF_regulons_rec_tf <- Receptor_target[rownames(Receptor_target) %in% TF,TF_rec]
  TF_regulons_rec_tf <- t(TF_regulons_rec_tf)
  TF_regulons_rec_tf <- as.data.frame(TF_regulons_rec_tf)
  TF_regulons_rec_tf <- TF_regulons_rec_tf[order(TF_regulons_rec_tf[,i]),]
  
  # TF_regulons_rec_tf <- moveColToTop(TF_regulons_rec_tf,i)
  TF_regulons_rec_tf <- TF_regulons_rec_tf %>% select(i,everything())  
  TF_regulons_rec_tf <- as.matrix(TF_regulons_rec_tf)
  
  
  
  p_receptor_regulon_network <- TF_regulons_rec_tf %>% 
    make_heatmap_ggplot("Receptors",paste(i,"_regulons",sep = ''), 
                        color = "blue",legend_position = "top", x_axis_position = "top",
                        legend_title = "Regulatory potential") + 
    scale_fill_gradient2() +  
    theme(axis.text.x = element_text(face = "italic"))
  p_receptor_regulon_network
  heat_plot[[i]] <- p_receptor_regulon_network
}


heat_plot$JUND
heat_plot$RARA
heat_plot$E2F3
heat_plot$GABPA
heat_plot$ZFX
heat_plot$FOXJ2
heat_plot$MTF2
heat_plot$SMAD3
heat_plot$TBX21
heat_plot$MAFF












i <- 'RUNX3'
TF <- c(i,met_regulons[[paste(i,'(+)',sep = '')]])
TF_sub <- Receptor_TF_sub[i,]
TF_sub <- t(TF_sub)
TF_sub <- as.data.frame(TF_sub)
TF_sub$receptor <- rownames(TF_sub)
TF_sub <- TF_sub[order(TF_sub[,i],decreasing = T),]
TF_rec <- head(TF_sub,10)
TF_rec <- TF_rec$receptor

TF_regulons_rec_tf <- Receptor_target[rownames(Receptor_target) %in% TF,TF_rec]
TF_regulons_rec_tf <- t(TF_regulons_rec_tf)
TF_regulons_rec_tf <- as.data.frame(TF_regulons_rec_tf)
TF_regulons_rec_tf <- TF_regulons_rec_tf[order(TF_regulons_rec_tf[,i]),]
TF_regulons_rec_tf <- TF_regulons_rec_tf %>% select(i,everything())
TF_regulons_rec_tf <- as.matrix(TF_regulons_rec_tf)



p_receptor_regulon_network <- TF_regulons_rec_tf %>% 
  make_heatmap_ggplot("Receptors",paste(i,"_regulons",sep = ''), 
                      color = "blue",legend_position = "top", x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradient2() +  
  theme(axis.text.x = element_text(face = "italic"))
  



p_receptor_regulon_network




