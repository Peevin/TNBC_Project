library(ggplot2)
library(nichenetr)

setwd('/Users/liupeiwen/BC/New_analysis/CD8-CXCL13/5.NicheNet/')

Receptor_target <- readRDS('OmniNetworks_NNformat/receptor_TF_matrixWithweights.rds')

TF <- read.table('OmniNetworks_NNformat/allTFs_hg38.txt',col.names = 'TF')


Receptor_TF <- Receptor_target[rownames(Receptor_target) %in% TF$TF,]

up_reg <- c('JUND','RARA','E2F3','GABPA','ZFX','FOXJ2',
            'MTF2','SMAD3','TBX21','MAFF')

Receptor_TF_sub <- as.data.frame(Receptor_TF[rownames(Receptor_TF) %in% up_reg,])

nrow(Receptor_TF_sub)*ncol(Receptor_TF_sub)
vioplot_data <- data.frame()
vioplot_data1 <- data.frame()
  
for (i in rownames(Receptor_TF_sub)) {
  for (j in colnames(Receptor_TF_sub)) {
    vectr <- c(i,j)
    vioplot_data <- rbind(vioplot_data,vectr) # fill data
    
  }
}
for (i in rownames(Receptor_TF_sub)) {
  for (j in colnames(Receptor_TF_sub)) {
    vectr <- Receptor_TF_sub[i,j]
    vioplot_data1 <- rbind(vioplot_data1,vectr) # fill data
    
  }
}



colnames(vioplot_data) <- c('TF','Receptor')
colnames(vioplot_data1) <- 'regulatory_potential'

vioplot_data$regulatory_potential <- vioplot_data1$regulatory_potential

nrow(vioplot_data[vioplot_data$regulatory_potential > 0.013, ])
length(unique(vioplot_data[vioplot_data$regulatory_potential > 0.013, 'TF']))
length(unique(vioplot_data[vioplot_data$regulatory_potential > 0.013, 'Receptor'])
)


ggplot(vioplot_data,aes(x=TF,y=regulatory_potential,fill=TF)) +
  geom_violin()

data_threshold <- vioplot_data[vioplot_data$regulatory_potential > 0.013, ]

my_receptor_TF <- Receptor_target[unique(data_threshold$TF),unique(data_threshold$Receptor)]
my_receptor_TF <- t(my_receptor_TF)
my_receptor_TF <- as.matrix(my_receptor_TF)

p_receptor_TF_network <- my_receptor_TF %>% 
  make_heatmap_ggplot("Receptors","TF", 
                      color = "blue",legend_position = "top", x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradient2() +  
  # ) +  
  theme(axis.text.x = element_text(face = "italic"))
p_receptor_TF_network

