library(nebula)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(showtext)
#> Warning: package 'nebula' was built under R version 4.1.2
# input raw count matrix，note:do not normalize the data
setwd('/Users/liupeiwen/BC/New_analysis/2.Nebula/')

cluster <- c('t_CD4-CXCL13', 't_CD4_Tcm-LMNA', 't_CD4_Tn-LEF1',
             't_CD4_Treg-FOXP3', 't_CD8-CXCL13', 't_Tact-XIST',
             't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683','t_Bmem-CD27',
             't_pB-IGHG1','t_macro-MMP9'
             )
cluster <- c('t_Tact-XIST','t_macro-MMP9')

for (j in cluster) {
  counts <- read.csv(paste0('./data/',j,'_nebula_mb.csv'), row.names = 1)
  # mb <- read.csv('./data/metabolism_genelist.csv',header = F)
  # counts <- counts[rownames(counts) %in% mb$V1,]
  
  count <- as.matrix(counts)
  
  batch <- read.csv(paste0('./data/',j,'_nebula_batch.csv'), header = F)
  
  label <- read.csv(paste0('./data/',j,'_nebula_label3.csv'), header = F)
  
  df = model.matrix(~V1, data=label)
  head(df)
  
  # data_g = group_cell(count=count,id=batch$V1,pred=df)
  
  re = nebula(count,batch$V1,pred=df)
  
  result <- re$summary
  
  result <- result[,c(8,2,6)]
  
  colnames(result) <- c('gene','logFC','pVal')
  rownames(result) <- result$gene
  result$adj.P <- p.adjust(result$pVal,method = 'BH')
  
  result$'-Log10(adj.P)' <- -log10(result$adj.P)
  # 由于这里差异分析是control比case，所以取个负值
  result$logFC <- -result$logFC
  
  result <-  result[order(result$logFC,decreasing = T),]
  
  
  result$expression <- ifelse(result$logFC >= 1 & result$adj.P < 0.05,"Up-regulated",
                              ifelse(result$logFC <= -1 & result$adj.P < 0.05,"Down-regulated","NS."))
  
  # 删除logFC很高，但NS的代谢基因
  result <- result %>%
    filter(!(abs(result$logFC) > 10 & result$adj.P > 0.05))
  
  
  # g <- ggplot(result, aes(logFC,`-Log10(adj.P)`,color=expression))+
  #   geom_point(size=0.1)+
  #   theme_classic()+
  #   scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
  #   geom_hline(yintercept=-log10(0.05),linetype=4, size=0.2)+
  #   geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
  #   labs(title = j)+
  #   theme(plot.title=element_text(hjust=0.5))
  # 
  # ggsave(paste0('figure/',j,'_volcano.pdf'),width = 8,height = 8,dpi=600)
  
  # write.csv(result,quote = F, row.names = F,
  #           file = paste0('./result/',j,'_nebula_result_total.csv'),
  # )
  
  # mb <- read.csv('data/metabolism_genelist_total.csv',header = F)
  # mb_result <- result[result$gene %in% mb$V1,]
  mb_result <- result
  
  g1 <- ggplot(mb_result, aes(logFC,`-Log10(adj.P)`,color=expression))+
    geom_point(size=0.3)+
    theme_classic()+
    scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
    geom_hline(yintercept=-log10(0.05),linetype=4, size=0.2)+
    geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
    theme(plot.title=element_text(hjust=0.5),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.y=element_text(size=20),
          axis.title.x=element_text(size=20),
          legend.position = 'none'
    ) +
    labs(title = j)+
    theme(plot.title=element_text(hjust=0.5))
  
  
  # 给5个上调和下调的基因标签
  adata <- bind_rows(
    mb_result %>%
      filter(expression == 'Down-regulated') %>%
      arrange(pVal,desc(abs(logFC)))%>%
      head(10),
    mb_result %>%
      filter(expression == 'Down-regulated') %>%
      tail(10)
    # mb_result %>%
    #   filter(expression == 'Down-regulated') %>%
    #   arrange(pVal,desc(abs(logFC)))%>%
    #   head(10),
    # mb_result %>%
    #   filter(logFC < 0) %>%
    #   tail(10)
  )
  
  g1+geom_text_repel(data = adata,
                     aes(logFC, `-Log10(adj.P)`,label=gene),
                     size = 5)  
  
  ggsave(paste0('figure/',j,'_volcano_mb.pdf'),width = 8,height = 8,dpi=600)
  
  write.csv(mb_result,quote = F, row.names = F,
            file = paste0('./result/',j,'_nebula_result_mb.csv'),
  )
  
  # 导出上调和下调的差异代谢基因
  up <- mb_result[mb_result$expression=='Up-regulated',]$gene
  up <- as.data.frame(up)
  write.table(up,paste0('./result/',j,'_up.csv'), col.names = F, row.names = F, quote = F)
  
}



counts <- read.csv('./data/t_CD8-CXCL13_nebula_total.csv', row.names = 1)
count <- as.matrix(counts)

batch <- read.csv('./data/t_CD8-CXCL13_nebula_batch.csv', header = F)

label <- read.csv('./data/t_CD8-CXCL13_nebula_label_cluster16.csv', header = F)

df = model.matrix(~V1, data=label)
head(df)

data_g = group_cell(count=count,id=batch$V1,pred=df)

re = nebula(count,batch$V1,pred=df)

result <- re$summary


result <- result[,c(8,2,6)]

colnames(result) <- c('gene','logFC','pVal')
rownames(result) <- result$gene

result$adj.P <- p.adjust(result$pVal,method = 'BH')

result$'-Log10(adj.P)' <- -log10(result$adj.P)
# 由于这里差异分析是control比case，所以取个负值
result$logFC <- -result$logFC



result <-  result[order(result$logFC,decreasing = T),]


ggplot(result, aes(logFC,`-Log10(adj.P)`))+
  geom_point()

result$expression <- ifelse(result$logFC >= 1 & result$adj.P < 0.05,"Up-regulated",
                          ifelse(result$logFC <= -1 & result$adj.P < 0.05,"Down-regulated","NS."))

# 删除logFC很高，但NS的代谢基因
result <- result %>%
  filter(!(abs(result$logFC) > 10 & result$pVal > 0.05))




g1 <- ggplot(result, aes(logFC,`-Log10(adj.P)`,color=expression))+
  geom_point(size=0.1)+
  theme_classic()+
  scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
  geom_hline(yintercept=-log10(0.05),linetype=4, size=0.2)+
  geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
  labs(title = 'CD8-CXCL13')+
  theme(plot.title=element_text(hjust=0.5))
g1




# result <- read.csv('result/t_CD8-CXCL13_nebula_result_total.csv')
# 导出gsea富集分析input数据
gsea_input <- result[, c(1,2)]
gsea_input = gsea_input[order(gsea_input$logFC,decreasing = T),]
write.table(gsea_input,quote = F, col.names = FALSE, row.names = F,
            file = 'result/CD8-CXCL13_gsea_data_cluster16_total.rnk',
            sep = '\t'
)

mb <- read.csv('data/metabolism_genelist_total.csv',header = F)
mb_result <- result[result$gene %in% mb$V1,]
g1 <- ggplot(mb_result, aes(logFC,`-Log10(P.Val)`,color=expression))+
  geom_point(size=0.1)+
  theme_classic()+
  scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
  geom_hline(yintercept=-log10(0.05),linetype=4, size=0.2)+
  geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
  labs(title = 'CD8-CXCL13')+
  theme(plot.title=element_text(hjust=0.5))
g1



# 给5个上调和下调的基因标签
adata <- bind_rows(
  result %>%
    filter(expression == 'Up-regulated') %>%
    arrange(pVal,desc(abs(logFC)))%>%
    head(5),
  result %>%
    head(10),
  # result %>%
  #   filter(logFC < 0) %>%
  #   tail(5),
  result %>%
    filter(expression == 'Down-regulated') %>%
    arrange(pVal,desc(abs(logFC)))%>%
    head(10)
)

g1+geom_text_repel(data = adata,
                   aes(logFC, `-Log10(P.Val)`,label=gene),
                   size = 2,show.legend = F
                   )  


# 导出上调和下调的差异代谢基因
up <- result[result$expression=='Up-regulated',]$gene
up <- as.data.frame(up)
# down <- rownames(result[result$expression=='Down-regulated',])
# down <- as.data.frame(down)
write.csv(result,quote = F, row.names = F,
          file = './data/t_CD8-CXCL13_nebula_result.csv',
)
write.table(up,'./result/t_CD8-CXCL13_up.csv', col.names = F, row.names = F, quote = F)
# write.table(down,'down_regulated_mbgene_CD8-CXCL13.csv', col.names = F, row.names = F, quote = F)
down <- result[result$logFC < 0,]$gene
down <- as.data.frame(down)
write.table(down,'./result/t_CD8-CXCL13_down.csv', col.names = F, row.names = F, quote = F)
gene <- result[result$logFC<0,]$gene %>% tail(100)


write.table(gene,'./result/CD8_CXCL13_logFC.csv', col.names = F, row.names = F, quote = F)



cluster <- c('t_CD4-CXCL13', 't_CD4_Tcm-LMNA', 't_CD4_Tn-LEF1',
             't_CD4_Treg-FOXP3', 't_CD8-CXCL13', 't_Tact-XIST',
             't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683')
for (j in cluster) {
  data <- read.csv(paste0('./result/',j,'_nebula_result_mb.csv'))
  up <- data[data$expression=='Up-regulated',]$gene
  write.table(up,paste0('./result/',j,'_up.txt'),quote = F,col.names = F,row.names = F)
}

sample_data$offset
data(sample_data)
sample_data$count[1:5,1:5]
sample_data$pred
