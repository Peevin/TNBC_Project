library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(ggplot2)

setwd('/Users/liupeiwen/BC/New_analysis/4.Scenic/')
getwd()
packageVersion("SCENIC")

scenicLoomPath <- file.path('./data/t_CD8-CXCL13_pyscenic_output.loom')
motifEnrichmentFile <- file.path('t_CD8-CXCL13_reg.csv')
file.exists(scenicLoomPath)

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





# read metabolic gene set
mb_gene <- read.table('./data/metabolism_genelist_total.csv')
mb_gene <- as.character(mb_gene$V1)

# filter regulons with metabolic genes
mb_regulons <- list()
for (i in 1:length(regulons)){
  for (j in mb_gene){
    if (j %in% regulons[[i]]){
      mb_regulons[[names(regulons[i])]] <- regulons[[i]]
      break
    }
    
  }
  
}
head(mb_regulons)

# calculate regulons size
regulons_count <- data.frame(row.names = names(mb_regulons))
for (i in 1:length(mb_regulons)){
  regulons_count$geneset[i] <- length(mb_regulons[[i]])
}

# read metabolic gene
mb_gene <- read.table('./data/metabolism_genelist_total.csv')
# read NEBULA DE data
DE_gene <- read.csv('./data/t_CD8-CXCL13_nebula_result_mb.csv')

# calculate the count of metabolic genes in each regulon
for (i in 1:length(mb_regulons)){
  regulons_count$mb_gene[i] <- sum(mb_gene$V1 %in% mb_regulons[[i]])
}

unique(DE_gene$expression)
gene_up_reg <- nrow(DE_gene[DE_gene$expression=='Up-regulated',])
# gene_down_reg <- nrow(DE_gene[DE_gene$expression=='Down-regulated',])
# calculate the up/down regulated genes in each regulon
for (i in 1:length(mb_regulons)){
  regulons_count$gene_up_reg[i] <- sum(DE_gene[DE_gene$expression == 'Up-regulated', 'gene'] %in% mb_regulons[[i]])
  # regulons_count$gene_down_reg[i] <- sum(DE_gene[DE_gene$expression == 'Down-regulated', 'gene'] %in% mb_regulons[[i]])
  
}

# calculate the p.Val: Hypergeometric test
for (i in 1:nrow(regulons_count)){
  regulons_count$up_P.Val[i] <- phyper(regulons_count$gene_up_reg[i]-1,
                                       regulons_count$mb_gene[i],
                                       3175-regulons_count$mb_gene[i],
                                       gene_up_reg, lower.tail=F)
  # regulons_count$down_P.Val[i] <- phyper(regulons_count$gene_down_reg[i]-1,
  #                                      regulons_count$mb_gene[i],
  #                                      3170-regulons_count$mb_gene[i],
  #                                      gene_down_reg, lower.tail=F)
}

# calculate adj.P
regulons_count$adj.P <- p.adjust(regulons_count$up_P.Val,method = 'BH')



# calculate the -log10(P.Val)
regulons_count$`-Log10(adj.P)` <- -log10(regulons_count$adj.P)
# regulons_count$`-Log10(p.Val)_down` <- -log10(regulons_count$down_P.Val)


regulons_count$regulon_name <- rownames(regulons_count)

regulons_enrich_up <- regulons_count[order(regulons_count$adj.P,decreasing = F),]
# regulons_enrich_down <- regulons_count[order(regulons_count$down_P.Val,decreasing = F),]
regulons_enrich_up_sig <- regulons_enrich_up[regulons_enrich_up$adj.P<0.01,]



ggplot(regulons_enrich_up_sig,
       aes(x=reorder(regulon_name,`-Log10(adj.P)`),y=`-Log10(adj.P)`,fill=`-Log10(adj.P)`)) +
  theme_classic() +
  geom_col() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  xlab(label = 'Regulons') +
  labs(title = 'CD8-CXCL13') +
  theme(plot.title=element_text(hjust=0.5)) +
  # scale_fill_manual(values = c("Sig.E"='#FF0000',"Non-Sig.E"='pink',
  #                              "Sig.NE"='blue',"Non-Sig.NE"='lightblue'))+
  coord_flip() 

ggsave('figure/t_CD8-CXCL13_phyper.01.pdf',width = 6,height = 8,dpi = 500)
 

regulons_count <- regulons_count[order((regulons_count$`-Log10(adj.P)`),decreasing = T),]

 
write.csv(regulons_count,'./output/t_CD8-CXCL13_DE_regulons_phyper_test.csv',
          quote = F)
write.csv(regulons_enrich_up_sig,'./output/t_CD8-CXCL13_DE_regulons_up_sig_phyper_test.csv',quote = F)











