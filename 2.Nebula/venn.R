library(ggvenn)
library(ggplot2)

setwd('/Users/liupeiwen/BC/New_analysis/')

TNBC <- read.csv('2.Nebula/result/t_CD8-CXCL13_nebula_result_mb.csv')

SCC <- read.csv('../SCC/2.NEBULA/result/CD8_ex_nebula_result_mb.csv')

a <- list(TNBC=TNBC[TNBC$expression=='Up-regulated',]$gene,
          SCC=SCC[SCC$expression=='Up-regulated',]$gene)

ggvenn(a,show_percentage = F)
