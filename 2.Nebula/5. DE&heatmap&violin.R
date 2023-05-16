library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggpubr)


setwd('/Users/liupeiwen/BC/New_analysis/2.Nebula/')
## 差异分析
cluster <- 't_CD8-CXCL13'
# 读取对数化的表达数据
data<-read.csv(paste0('data/',cluster,'_nebula_log.csv'),
             header = 1,row.names = 1,sep=',')

# 读取case组和contral组的标签文件
group<-read.csv(paste0('data/',cluster,'_nebula_label.csv'),
                header = F)
group$V2 <- ifelse(group$V1=='ICI-opposite','ICIopposite','ICIsupport')
group <- group$V2


suppressMessages(library(limma))
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)


contrast.matrix<-makeContrasts("ICIopposite-ICIsupport",levels=design)



##step1
fit <- lmFit(data,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

nrDEG = nrDEG[order(nrDEG$logFC,decreasing = T),]
data = nrDEG

rnk = data.frame(row.names = rownames(nrDEG),nrDEG[,1])

#
write.table(rnk,quote = F, col.names = FALSE,
          file = paste('/Users/liupeiwen/toomanycells/analysis_result', globalc, treatment, subc, node, gseadata, sep='/'),
          sep = '\t'
          )


write.table(rownames(adata[adata$expression=='Up-regulated',]),
            file = paste('/Users/liupeiwen/toomanycells/analysis_result', globalc, treatment, subc, node, godata, sep='/'),
            sep = ' ',
            quote = F,
            row.names = F,
            col.names = F
            )

data$feature <- rownames(nrDEG)


data$'-log10(adj.P.Val)' <- -log10(data$adj.P.Val)
colnames(data)[1] = 'log2FC'
ggplot(data, aes_string('log2FC','-log10(adj.P.Val)'))+
  geom_point()

data$expression <- ifelse(data$log2FC >= 1 & data$adj.P.Val < 0.05,"Up-regulated",
                          ifelse(data$log2FC <=-1 & data$adj.P.Val < 0.05,"Down-regulated","NS."))

g1 <- ggplot(data, aes(log2FC,-log10(adj.P.Val),color=expression))+
  geom_point(size=0.5)+
  theme_classic()+
  scale_color_manual(values = c("Down-regulated"='blue',"NS."='black',"Up-regulated"='red'))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  labs(title = 't_CD8-CXCL13')+
  theme(plot.title=element_text(hjust=0.5))
g1


# 给5个上调和下调的基因标签
adata <- bind_rows(
  data %>%
    filter(expression == 'Up-regulated') %>%
    arrange(adj.P.Val,desc(abs(log2FC)))%>%
    head(5),
  data %>%
    filter(expression == 'Down-regulated') %>%
    arrange(adj.P.Val,desc(abs(log2FC))) %>%
    head(5)
)

g1+geom_text_repel(data = adata,
                   aes(log2FC, -log10(adj.P.Val),label=feature),
                   size = 2)  



## heatmap plot

# 读取scale之后的分组的平均表达矩阵
em <- read.csv('/Users/liupeiwen/toomanycells/analysis_result/Tcell/ICB/t_CD8_Tem-GZMK_df/mean.scale.em.ICB.csv',header = 1,row.names = 1)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))


up <- c('GLA','SESN1','FUS','UBB','DDIT4')


GLA <- c('SPI1','USF2','USF1','NFE2','MYC','FOS','ETS1','STAT1')

em1 <- em[,colnames(em) %in% up]
em1 <- em1[, up]


pheatmap::pheatmap(em1,cluster_cols = F,
                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   cluster_rows = F,angle_col = 90,
                   cellwidth = 20,
                   cellheight = 35,
                   fontsize = 8,
                   #gaps_col = c(5,7,8,11),
                   gaps_row = 2,
                   breaks = bk,
                   #labels_col = 'ASS1',
                   main = 't_CD8_Tem-GZMK',
                   
)



## violin plot

# 读取log之后的表达数据（原始样本）
me <- read.csv('/Users/liupeiwen/toomanycells/analysis_result/Tcell/ICB/t_CD8_Tem-GZMK_df/em.log.ICB.csv',header = 1,row.names = 1)

setwd('/Users/liupeiwen/toomanycells/analysis_result/Tcell/ICB/t_CD8_Tem-GZMK_df/node5/figure')

# 设置显著性检验的组别
compaired <- list(c('SD_Pre','SD_Post'),c('SD_Pre','PR_Pre'))


for(i in up){
  g1 <- ggplot(me,aes_string(x='group',y=i, fill='group'))+
    geom_violin()+
    theme_classic()+
    theme(legend.position = 'none', plot.title=element_text(hjust=0.5))+
    labs(title = paste('anti-PDL1+Chemo',i,sep='_'))+
    scale_x_discrete(limits=c('SD_Pre','SD_Post','PR_Pre','PR_Post'))+
    #scale_y_continuous(limits=c(0,11))+
    geom_signif(comparisons = compaired, 
                step_increase = 0.1,
                test =wilcox.test ,map_signif_level = T)
  
  
  ggsave(file=paste('violin_ICB_',i,'.svg',sep=''), plot=g1, width=10, height=6)
}







