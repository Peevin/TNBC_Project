library(Seurat)
library(SeuratDisk)
library(dplyr)
library(nebula)




setwd('/Users/liupeiwen/BC/New_analysis/0.Signature/')

# Convert('../1.Toomanycells_result/Tcell/data/Tcell_PDL1.h5ad', "h5seurat",
#         overwrite = TRUE,assay = "RNA")
adata <- LoadH5Seurat("../1.Toomanycells_result/Tcell/data/Tcell_PDL1.h5seurat")
adata <- adata[,adata$Sub_Cluster=='t_CD4_Treg-FOXP3']
ann <- read.csv('../1.Toomanycells_result/Tcell/result/cluster/t_CD4_Treg-FOXP3_cluster.csv',row.names = 1)
adata <- AddMetaData(adata,ann)
# load metabolic genes
mb <- read.table('./data/metabolism_genelist_total.csv')
mb <- as.vector(mb$V1)


cluster <- c(4,8,9,11,16,18,19)
cluster <- as.character(cluster)
adata$cluster <- as.character(adata$cluster)
# add metadata for nebula
for (j in cluster) {
  adata@meta.data[paste0('label',j)] <- ifelse(adata$cluster==j,'case','control')
}


# choose the metabolic gene expression matrix
mb_adata <- adata[rownames(adata) %in% mb]

# get normalized expression matrix
adata <- NormalizeData(adata,scale.factor = 1e6)
adata <- adata[rownames(adata) %in% mb]


# dataframe of the metabolic gene expression matrix
counts <- as.data.frame(mb_adata@assays$RNA@data)
count <- as.matrix(counts)

batch <- mb_adata@meta.data['Patient']

# select genes with detection rate >= 0.2
det <- list()
for (j in cluster) {
  adata$label <- adata@meta.data[paste0('label',j)]
  data <- subset(adata,label=='case')
  matrix <- as.data.frame(GetAssayData(data, slot = "data"))
  df <- data.frame(row.names = rownames(matrix))

  df$det.rate <- rowSums(matrix != 0) / ncol(matrix)
  df$mean <- rowMeans(matrix)
  signature <- rownames(df[df$det.rate>=0.2,])
  det[[j]] <- signature
}


sig <- list()
# nebula
for (j in cluster) {
  mb_adata$label <- mb_adata@meta.data[paste0('label',j)]
  label <-  mb_adata@meta.data['label']
  
  df = model.matrix(~label, data=label)
  data_g = group_cell(count=count,id=label$label,pred=df)
  
  re = nebula(data_g$count,data_g$id,pred=data_g$pred)
  
  result <- re$summary
  
  result <- result[,c(8,2,6)]
  
  colnames(result) <- c('gene','log2FC','pVal')
  rownames(result) <- result$gene
  result$adj.p <- p.adjust(result$pVal,method = 'BH')
  
  result$'-Log10(adj.P)' <- -log10(result$adj.p)
  # 由于这里差异分析是control比case，所以取个负值
  result$log2FC <- -result$log2FC
  
  result <-  result[order(result$log2FC,decreasing = T),]
  
  result$expression <- ifelse(result$log2FC >= 1 & result$adj.p < 0.05,"Up-regulated",
                              ifelse(result$log2FC <= -1 & result$adj.p < 0.05,"Down-regulated","NS."))
  
  # 删除log2FC很高，但NS的代谢基因
  result <- result %>%
    filter(!(abs(result$log2FC) > 10 & result$adj.p > 0.05))

  write.csv(result,paste0('./result/CD4_Treg-FOXP3_',j,'.csv'),quote = F,row.names = F)

  # result <- result[result$expression=='Up-regulated',]
  # result <- result %>% filter(log2FC>0) %>% filter(adj.p<0.05)
  # 
  # signature <- result[result$gene %in% det[[j]],]
  # signature <- signature %>% head(20)
  # signature <-  signature %>% mutate(order = row_number())
  # 
  # adata$label <- adata@meta.data[paste0('label',j)]
  # data <- subset(adata,label=='case')
  # matrix <- as.data.frame(GetAssayData(data, slot = "data"))
  # df <- data.frame(row.names = rownames(matrix))
  # 
  # df$det.rate <- rowSums(matrix != 0) / ncol(matrix)
  # df$mean <- rowMeans(matrix)
  # df$gene <- rownames(df)
  # df <- df[rownames(df) %in% signature$gene,]
  # 
  # df2 <- merge(signature,df,by='gene')
  # df2 <- df2[order(df2$log2FC,decreasing = T),]
  # 
  # sig[[j]] <- df2
  
}


for (j in cluster) {
  result <- read.csv(paste0('result/CD4_Treg-FOXP3_',j,'_signature.csv'))

  signature <- result[result$gene %in% det[[j]],]
  # signature <- result
  # signature <- signature %>% head(20)
  signature <-  signature %>% mutate(order = row_number())

  adata$label <- adata@meta.data[paste0('label',j)]
  data <- subset(adata,label=='case')
  matrix <- as.data.frame(GetAssayData(data, slot = "data"))
  df <- data.frame(row.names = rownames(matrix))

  df$det.rate <- rowSums(matrix != 0) / ncol(matrix)
  df$mean <- rowMeans(matrix)
  df$gene <- rownames(df)
  df <- df[rownames(df) %in% signature$gene,]

  df2 <- merge(signature,df,by='gene')
  df2 <- df2[order(df2$log2FC,decreasing = T),]

  sig[[j]] <- df2
}
Reduce(union,sig)

saveRDS(sig,'./data/CD4_Treg-FOXP3_sig.rds')

a <- data.frame()
for (i in cluster) {
  a <- rbind(a,sig[[i]][c('gene','order')])
}

a$weight <- 1 / a$order

weight <- a %>% group_by(gene) %>% summarize(sum_value = sum(weight))
weight$scale_w <- (weight$sum_value - min(weight$sum_value)) / (max(weight$sum_value) - min(weight$sum_value))
b <- as.vector(weight$scale_w)
names(b) <- weight$gene

write.csv(weight[c('gene','scale_w')],'data/CD4_Treg-FOXP3_sig_weight.csv',quote = F,row.names = F)

# ---------------------------------------------------------------


matrix <- as.data.frame(GetAssayData(adata, slot = "data"))
matrix <- as.data.frame(t(matrix))

adata$score <- 0
for (i in rownames(matrix)) {
  c <- matrix[i,]
  score <- 0
  for (j in weight$gene) {
    score <- c[[j]] * b[[j]] + score
  }
  adata$score[i] <- score
}

z <- adata@meta.data %>% group_by(cluster) %>% summarize(sum_value = mean(score))
adata@meta.data %>% group_by(cluster) %>% summarize(sum_value = mean(score))




Convert('./data/CD8_EX_SCC.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
adata2 <- LoadH5Seurat("./data/CD8_EX_SCC.h5seurat",meta.data=T)
ann <- read.csv('data/CD8_EX_SCC_obs.csv',row.names = 1)
adata2 <- AddMetaData(adata2,ann)
adata2 <- NormalizeData(adata2,scale.factor = 1e6)

matrix2 <- as.data.frame(GetAssayData(adata2, slot = "data"))
matrix2 <- as.data.frame(t(matrix2))

matrix2 <- matrix2[,colnames(matrix2) %in% weight$gene]

adata2$score <- 0
for (i in rownames(matrix2)) {
  c <- matrix2[i,]
  score <- 0
  for (j in colnames(matrix2)) {
    score <- c[[j]] * b[[j]] + score
  }
  adata2$score[i] <- score
}
z2 <- adata2@meta.data %>% group_by(cluster_TMC) %>% summarize(sum_value = mean(score))

n <- t(adata2[['score']])
write.csv(n,'data/CD8_EX_sig_score.csv',quote = F,col.names = F)
adata2$score