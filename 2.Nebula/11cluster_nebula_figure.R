library(nebula)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(showtext)

setwd('/Users/liupeiwen/BC/New_analysis/2.Nebula/')

# file_list <- list.files('result/',pattern = "mb.csv$")

cluster <- c('CD4_Tcm-LMNA','CD4_Tn-LEF1','CD4-CXCL13',
             'CD4_Treg-FOXP3','CD8-CXCL13','CD8_Tem-GZMK',
             'CD8_Trm-ZNF683','Tact-XIST','macro-MMP9','Bmem-CD27',
             'pB-IGHG1')
df <- data.frame()
top10_all <- data.frame()
for (i in cluster) {
  data <- read.csv(paste0('result/t_',i,'_nebula_result_mb.csv'))
  data$cluster <- i
  data <- data[data$expression %in% c('Up-regulated','Down-regulated'),]
  df <- rbind(df,data)
  top10 <- data %>% top_n(10,abs(logFC))
  top10_all <- rbind(top10_all,top10)
}

df$cluster <- factor(df$cluster,levels = c('CD4_Tcm-LMNA','CD4_Tn-LEF1','CD4-CXCL13',
                                           'CD4_Treg-FOXP3','CD8-CXCL13','CD8_Tem-GZMK',
                                           'CD8_Trm-ZNF683','Tact-XIST','macro-MMP9','Bmem-CD27',
                                           'pB-IGHG1'))
top10_all$cluster <- factor(top10_all$cluster,levels = c('CD4_Tcm-LMNA','CD4_Tn-LEF1','CD4-CXCL13',
                                                         'CD4_Treg-FOXP3','CD8-CXCL13','CD8_Tem-GZMK',
                                                         'CD8_Trm-ZNF683','Tact-XIST','macro-MMP9','Bmem-CD27',
                                                         'pB-IGHG1'))




#叠加每个Cluster Top10基因散点(将散点适当放大强调）：
p <- ggplot()+
  geom_jitter(data = df,
              aes(x = cluster, y = logFC, color = expression),
              size = 0.4,
              width =0.4)

p



#添加X轴的cluster色块标签：
dfcol<-data.frame(x=c(1:11),
                  y=0,
                  label=c('CD4_Tcm-LMNA','CD4_Tn-LEF1','CD4-CXCL13',
                          'CD4_Treg-FOXP3','CD8-CXCL13','CD8_Tem-GZMK',
                          'CD8_Trm-ZNF683','Tact-XIST','macro-MMP9','Bmem-CD27',
                          'pB-IGHG1'))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F","#7E61487F","#82B0D2","#FFBE7A")
p2 <- p + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1,
                     color = "black",
                     fill = mycol,
                     alpha = 1,
                     show.legend = F)
p2


#给每个Cluster差异表达前Top10基因加上标签：
p3 <- p2+
  geom_text_repel(
    data=top10_all,
    aes(x=cluster,y=logFC,label=gene),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"),
    max.overlaps = 10
  )
p3

p4 <- p3 +
  scale_color_manual(values = c("Down-regulated"='red',"Up-regulated"='blue'),
                     labels=c('Up-regulated in non-NR-E','Up-regulated in NR-E'))

p4


#添加cluster数字：
p5 <- p4+
  labs(x="Cluster")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3.5,
            color ="white")
p5

#自定义主题美化：
p6 <- p5+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 10)
  )
p6


