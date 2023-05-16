library(ggplot2)
library(dplyr)
library(ggrepel)


setwd('/Users/liupeiwen/BC/New_analysis/1.Toomanycells_result/Tcell/')

cluster <- c('t_CD4-CXCL13', 't_CD4_Tcm-LMNA', 't_CD4_Tn-LEF1',
             't_CD4_Treg-FOXP3', 't_CD8-CXCL13', 't_CD8_MAIT-KLRB1',
             't_CD8_Teff-GNLY', 't_CD8_Tem-GZMK', 't_CD8_Trm-ZNF683',
             't_Tact-IFI6', 't_Tact-XIST', 't_Tprf-MKI67')
cluster <- c('t_CD8-CXCL13')
for(j in cluster){
  
  # read count data
  adata <- read.csv(paste('data/cluster_count/',j,'_cluster_count.csv',sep = ''))
  
  # prop.test
  for (i in rownames(adata)){
    if(adata[i,'SD_Prop'] == adata[1,'SD_Prop']){
      a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                     c(adata[i,'Count'],adata[1,'Count']),alternative = 'two.sided')
      adata[i,'p_Val'] <- a$p.value
      adata[i,'alternative'] <- 'two.sided'
    }
    else if(adata[i,'SD_Prop'] < adata[1,'SD_Prop']){
      a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                     c(adata[i,'Count'],adata[1,'Count']),alternative = 'less')
      adata[i,'p_Val'] <- a$p.value
      adata[i,'alternative'] <- 'less'
    }
    else{
      a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                     c(adata[i,'Count'],adata[1,'Count']),alternative = 'greater')
      adata[i,'p_Val'] <- a$p.value
      adata[i,'alternative'] <- 'greater'
    }
  }
  
  adata$`-log10pVal` <- -log10(adata$p_Val)
  adata$Cluster = paste('Cluster',adata$Cluster,sep = '')
  adata$Cluster <- factor(adata$Cluster,levels = as.vector(adata$Cluster))
  
  
  
  for (i in 1:nrow(adata)) {
    adata$Fold_Change[i] = adata[i,'SD_Prop'] / adata[1,'SD_Prop']
    adata$`log2FC`[i] = log2(adata$Fold_Change[i])
  }
  
  
  # volcano plot
  # g1 <- ggplot(adata[-1,], aes(log2FC,`-log10pVal`,fill=Count))+
  #   geom_point(size=2,shape=21)+
  #   scale_fill_distiller(palette = "Reds",direction = 1) +
  #   theme_classic()+
  #   # scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
  #   geom_hline(yintercept=-log10(1e-20),linetype=4, size=0.2)+
  #   geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
  #   labs(title = j)+
  #   theme(plot.title=element_text(hjust=0.5))
  # label = adata %>% 
  #   filter(log2FC > 1) %>% 
  #   filter(p_Val < 1e-20)
  # 
  # g2 <- g1+geom_text_repel(data = label,
  #                    aes(log2FC, `-log10pVal`,label=Cluster),
  #                    size = 4)  
  # # 选择火山图的保存路径
  # ggsave(paste('figure/propotion/',j,'_volcano.pdf',sep=''),
  #        width=8,height=8,dpi=600,device='pdf')
  
  adata$label <- ifelse(adata$SD_Prop >= adata[1,'SD_Prop']+0.25 & adata$p_Val < 1e-4,"NR_enrich",
                        ifelse(adata$SD_Prop <= adata[1,'SD_Prop']-0.25 & adata$p_Val < 1e-4,"R_enrich","None_enrich"))
  colnames(adata)[colnames(adata)=='SD_Prop'] <- 'NR_Prop'
  # 火山图2
  g3 <- ggplot(adata[-1,], aes(NR_Prop,`-log10pVal`,fill=label))+
    
    theme_classic()+
    geom_hline(yintercept=-log10(1e-3),linetype=4, size=0.2)+
    geom_vline(xintercept=c(max(adata[1,'NR_Prop']-0.25,0),min(adata[1,'NR_Prop']+0.25,1)),
               linetype=4, size=0.2)+
    scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
    geom_vline(xintercept=adata[1,'NR_Prop'],linetype=2) +
    geom_rect(aes(xmin=max(adata[1,'NR_Prop']-0.25,0), xmax=min(adata[1,'NR_Prop']+0.25,1), 
                  ymin=-Inf, ymax=Inf),fill='grey',alpha = 0.01)+
    labs(title = j)+
    theme(plot.title=element_text(hjust=0.5)) +
    geom_point(size=5,shape=21,alpha=.6)+
    scale_fill_manual(values = c("R_enrich"='red',"None_enrich"='darkgrey',"NR_enrich"='blue'))
    
  
  label = adata %>% 
    filter(NR_Prop > adata[1,'NR_Prop']+0.25) %>% 
    filter(p_Val < 1e-3)
  
  g4 <- g3+geom_text_repel(data = label,
                           aes(NR_Prop, `-log10pVal`,label=Cluster),
                           size = 5)

  ggsave(paste('figure/volcano/',j,'_threshold.pdf',sep=''),
         width=8,height=8,dpi=900,device='pdf')
  
  
# 
  # plot
  g <- ggplot(adata[-1,],aes(x=Cluster,y=NR_Prop,fill=Count)) +
    geom_segment(aes(x=Cluster,xend=Cluster,y=adata[1,'NR_Prop'],yend=NR_Prop),color='grey',linetype='dotted') +
    geom_point(aes(size=Count),shape=21,color='black') +
    # scale_color_manual(values = c("grey","black")) +
    scale_fill_distiller(palette = "Reds",direction = 1) +
    theme_classic() +
    theme(axis.text.x = element_text (angle = 90, hjust = 1),plot.title = element_text (hjust = 0.5)) +
    ggtitle(j) +
    geom_hline(yintercept=adata[1,'NR_Prop'],linetype=2, size=0.2)

  # 选择比例图的保存路径
  ggsave(paste('figure/propotion/',j,'_Prop.pdf',sep=''),
         width=8,height=,dpi=600,device='pdf')

  # write.csv(adata,paste('./result/cluster/',j,'_cluster_count.csv',sep=''),quote = F,row.names = F)
  
}

g3
# ---------------------------以下为测试代码-----------------------

adata <- read.csv('data/cluster_count/t_CD4_Treg-FOXP3_cluster_count.csv')

for (i in rownames(adata)){
  if(adata[i,'SD_Prop'] == adata[1,'SD_Prop']){
    a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                   c(adata[i,'Count'],adata[1,'Count']),alternative = 'two.sided')
    adata[i,'p_Val'] <- a$p.value
    adata[i,'alternative'] <- 'two.sided'
  }
  else if(adata[i,'SD_Prop'] < adata[1,'SD_Prop']){
    a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                   c(adata[i,'Count'],adata[1,'Count']),alternative = 'less')
    adata[i,'p_Val'] <- a$p.value
    adata[i,'alternative'] <- 'less'
  }
  else{
    a <- prop.test(c(adata[i,'SD_Count'],adata[1,'SD_Count']),
                   c(adata[i,'Count'],adata[1,'Count']),alternative = 'greater')
    adata[i,'p_Val'] <- a$p.value
    adata[i,'alternative'] <- 'greater'
  }
  
 
}



adata$`-log10pVal` <- -log10(adata$p_Val)
adata$Cluster = paste('Cluster',adata$Cluster,sep = '')
adata$Cluster <- factor(adata$Cluster,levels = as.vector(adata$Cluster))






for (i in 1:nrow(adata)) {
    adata$Fold_Change[i] = adata[i,'SD_Prop'] / adata[1,'SD_Prop']
    adata$`log2FC`[i] = log2(adata$Fold_Change[i])
}

adata$label <- ifelse(adata$SD_Prop >= adata[1,'SD_Prop']+0.25 & adata$p_Val < 1e-4,"SD_enrich",
                            ifelse(adata$SD_Prop <= adata[1,'SD_Prop']-0.25 & adata$p_Val < 1e-4,"PR_enrich","None_enrich"))


mad(adata$Count[-1],  constant = 1, na.rm = FALSE, low = FALSE, high = FALSE)
median(adata$Count[-1])


# volcano plot
g1 <- ggplot(adata[-1,], aes(log2FC,`-log10pVal`,fill=`-log10pVal`)) +
  geom_point(size=2,shape=21) +
  scale_fill_distiller(palette = "Reds",direction = 1) +
  theme_classic()+
  # scale_color_manual(values = c("Down-regulated"='red',"NS."='darkgrey',"Up-regulated"='blue'))+
  geom_hline(yintercept=-log10(1e-20),linetype=4, size=0.2)+
  geom_vline(xintercept=c(-1,1),linetype=4, size=0.2)+
  labs(title = 'CD8-CXCL13')+
  theme(plot.title=element_text(hjust=0.5))
g1

label = adata %>% 
  filter(log2FC > 1) %>% 
  filter(p_Val < 1e-20)

g1+geom_text_repel(data = label,
                   aes(log2FC, `-log10pVal`,label=Cluster),
                   size = 2)  




g3 <- ggplot(adata[-1,], aes(SD_Prop,`-log10pVal`,fill=label))+
  geom_point(size=2,shape=21,alpha=.6)+
  theme_classic()+
  scale_fill_manual(values = c("PR_enrich"='red',"None_enrich"='darkgrey',"SD_enrich"='blue'))+
  geom_hline(yintercept=-log10(1e-4),linetype=4, size=0.2)+
  geom_vline(xintercept=c(max(adata[1,'SD_Prop']-0.25,0),min(adata[1,'SD_Prop']+0.25,1)),
             linetype=4, size=0.2)+
  scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
  geom_vline(xintercept=adata[1,'SD_Prop'],linetype=2) +
  geom_rect(aes(xmin=max(adata[1,'SD_Prop']-0.25,0), xmax=min(adata[1,'SD_Prop']+0.25,1), 
                ymin=-Inf, ymax=Inf),fill='grey',alpha = .01)+
  labs(title = 'CD8-CXCL13')+
  theme(plot.title=element_text(hjust=0.5))

label = adata %>% 
  filter(SD_Prop > adata[1,'SD_Prop']+0.25) %>% 
  filter(p_Val < 1e-4)

g4 <- g3+geom_text_repel(data = label,
                         aes(SD_Prop, `-log10pVal`,label=Cluster),
                         size = 3)
g4

# plot
g <- ggplot(adata[-1,],aes(x=Cluster,y=SD_Prop,fill=Count)) +
        geom_segment(aes(x=Cluster,xend=Cluster,y=adata[1,'SD_Prop'],yend=SD_Prop),color='grey',linetype='dotted') +
        geom_point(aes(size=Count),shape=21,color='black') +
        # scale_color_manual(values = c("grey","black")) +
        scale_fill_distiller(palette = "Reds",direction = 1) +
        theme_classic() +
        theme(axis.text.x = element_text (angle = 90, hjust = 1),plot.title = element_text (hjust = 0.5)) +
        ggtitle('CD8-CXCL13') +
        geom_hline(yintercept=adata[1,'SD_Prop'],linetype=2, size=0.2,show)
g


write.csv(adata,'./result/cluster_count.csv',quote = F,row.names = F)


















  df <- data.frame(
  group = c('PR','SD'),
  value = c(4806,598)
)

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pie

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/100)), size=5)

pie + scale_fill_manual(values=c('#D1352B','#4A7CB3')) + blank_theme +
  theme(axis.text.x=element_blank())
