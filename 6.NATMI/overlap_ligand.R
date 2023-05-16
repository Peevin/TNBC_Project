library(readxl)
library(tidyverse)

setwd('/Users/liupeiwen/BC/New_analysis/6.NATMI/')

files = list.files("data/lig_rec/", pattern = "xlsx", full.names = TRUE)
df = map_dfr(files, read_xlsx) 

lig_act <- read.csv('data/ligand_activity.csv')
lig_act <- lig_act %>% 
  filter(pearson>0)

df <- df %>% 
  filter(`target cell`=='t_CD8-CXCL13') %>% 
  filter(ligand %in% lig_act$test_ligand)

myeloid <- c("t_macro-CCL2",
             "t_macro-CD24","t_macro-CFD" ,"t_macro-CX3CR1","t_macro-FOLR2",
             "t_macro-IFI27","t_macro-IGFBP7","t_macro-IL1B","t_macro-IL1RN", 
             "t_macro-MGP","t_macro-MKI67","t_macro-MMP9","t_macro-SLC40A1",
             "t_macro-SPP1","t_macro-TUBA1B")
df <- df %>% 
  filter(`sending cell` %in% myeloid)

df$lig_rec <- paste(df$ligand,df$receptor,sep = '-')


ls <- list()

for (i in unique(df$`sending cell`)) {
  sub <- df[df$`sending cell`==i,]
  ls[[i]] <- sub$lig_rec
}

Reduce(intersect,ls)


