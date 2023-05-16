library(SCopeLoomR)
library(SCENIC)
library(dplyr)


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

# read NEBULA DE data
DE_gene <- read.csv('../2.Nebula/result/t_CD8-CXCL13_nebula_result_mb.csv')
up <- DE_gene[DE_gene$expression=='Up-regulated',]$gene

sig_reg <- read.csv('./output/t_CD8-CXCL13_DE_regulons_up_sig_phyper_test.csv',row.names = 1)

mb_regulons <- mb_regulons[names(mb_regulons) %in% sig_reg$regulon_name]

df <- as.data.frame(x=mb_regulons$`BATF(+)`)
colnames(df) <- 'Target'
df$TF <- 'BATF'
df <- df[,c(2,1)]

ls <- list()
for (j in names(mb_regulons)) {
  df <- as.data.frame(mb_regulons[[j]])
  colnames(df) <- 'Target'
  df$TF <- j
  df <- df[,c(2,1)]
  ls[[j]] <- df
}

for (i in names(ls)) {
  write.csv(ls[[i]],paste0('output/TF_Target/regulon_',i,'.csv'),quote = F,row.names = F,col.names = F)
}


ls <- do.call(rbind,ls)
rownames(ls) <- NULL
ls <- ls[ls$Target %in% up,]
ls$TF <- substring(ls$TF,1,nchar(ls$TF)-3)
ls$TF_Target <- paste(ls$TF,ls$Target,sep = '_')

adj <- read.csv('./output/t_CD8-CXCL13_adj.tsv',sep = '\t')
adj$TF_Target <- paste(adj$TF,adj$target,sep = '_')

adj <- adj %>% 
  filter(adj$TF_Target %in% ls$TF_Target)

colnames(adj) <- c('source','target','importance','TF_Target')
adj$relation <- adj$source

a <- as.data.frame(table(adj$target))
a <- a[a$Freq>=5,]

adj <- adj[adj$target %in% a$Var1,]



# dir.create('./output/TF_Target/')
write.csv(adj,'./output/TF_Target/CD8-CXCL13_TF_Target_5.csv',quote = F,row.names = F)


TF <- data.frame(unique(adj$source))
colnames(TF) <- 'Id'
TF$Label <- TF$Id
TF$meta <- 'TF'

Target <- data.frame(unique(adj$target))
colnames(Target) <- 'Id'
Target$Label <- Target$Id
Target$meta <- 'Target'


df2 <- rbind(TF,Target)




write.csv(df2,'./output/TF_Target/CD8-CXCL13_label_5.csv',row.names = F,quote = F)
