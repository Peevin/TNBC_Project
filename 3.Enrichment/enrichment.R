library(ggplot2)
library(ggh4x)
library(tidyverse)
library(forcats)

setwd('/Users/liupeiwen/BC/New_analysis/3.Enrichment/')

adata <- read.table('./result/CD8-CXCL13.txt',sep = '\t',header = 1)

adata <- adata[c('geneSet','enrichmentRatio', 'FDR')]
adata$enrichmentRatio <- as.numeric(adata$enrichmentRatio)
adata$FDR <- as.numeric(adata$FDR)

ann <- read.csv('./data/KEGG_ann.csv')

adata <- left_join(adata,ann,by=("geneSet"="geneSet"))

adata$geneSet <- substring(adata$geneSet, 6,nchar(adata$geneSet))



adata <- adata[order(adata$enrichmentRatio,decreasing = T),]
adata$geneSet  <- factor(adata$geneSet ,levels = as.vector(rev(adata$geneSet)))


col <- c("Genetic Information Processing"='white',
         "Glycan biosynthesis and metabolism"="#e30039",
         "Energy metabolism"="#f89588",
         "Carbohydrate metabolism"="#f8cb7f",
         "Human Diseases"="white",
         "Nucleotide metabolism"="#76da91",
         "Amino acid metabolism"="#ff6600",
         "Lipid metabolism"="#0e2c82",
         "Cellular Processes" ="white",
         "Xenobiotics biodegradation and metabolism"="white",
         "Organismal Systems" ="white"
         
         )


B <- adata




Tumor.cellularity <- B$geneSet %>% as.data.frame() %>%
  mutate(group=B$annotate) %>%
  mutate(p="")%>%
  ggplot(aes(p,.,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  scale_fill_manual(values = col)+
  theme_void()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none'
        )+
  labs(fill = "Tumor.cellularity")



p1 <- ggplot(adata,aes(x=geneSet,y=enrichmentRatio, fill=FDR)) +
  geom_col(linewidth=0.7) +

  # scale_color_manual(values = col) +
  scale_fill_distiller(palette = "Blues") +
  coord_flip() +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=20)
        # legend.position = 'none'
  )


p1

bottom <- ggplotGrob(Tumor.cellularity)
# p1+annotation_custom(bottom,xmin=-1,xmax=141.5,ymin=-0.03,ymax=0.01)
p1+annotation_custom(bottom,ymin=-0.2,ymax=0)

c <- B$geneSet %>% as.data.frame() %>%
  mutate(group=B$annotate) %>%
  mutate(p="")%>%
  ggplot(aes(p,.,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  scale_fill_manual(values = col)+
  theme_void()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank()
        # legend.position = 'none'
  )+
  labs(fill = "Tumor.cellularity")

c
