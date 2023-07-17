library(this.path)
library(tidytree)
library(ggtree)
library(openxlsx)
library(treeio)
library(glue)
this.dir()

setwd(this.dir())

taxa_tbl <- grep('_taxa_table.xlsx$', list.files(), value = TRUE)
if (length(taxa_tbl) == 2) {
    df <- read.xlsx(taxa_tbl[2])
} else {
    df <- read.xlsx(taxa_tbl[1])
}

tre <- read.iqtree('06_iqtree/iqtree_ml.treefile')

ggtree(tre, ) %<+% df + 
    geom_tiplab(aes(subset = (New == 1 & Type == 1), 
                    label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")), 
                color = 'darkblue', parse = T, ) + 
    geom_tiplab(aes(subset = (New == 0 & Type == 1), 
                    label = glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection})~bold({Number1})")), 
                color = 'black', parse = T, ) + 
    geom_tiplab(aes(subset = (New == 0 & Type == 0), 
                    label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection,'~', Number1)), 
                color = 'black', parse = T, ) + 
    geom_tiplab(aes(subset = (New == 1 & Type == 0), 
                    label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection,'~', Number1)), 
                color = 'darkblue', parse = T, ) + 
    
    geom_text2(aes(label= UFboot, subset = (UFboot >=95)), vjust = 0, hjust = 1, size = 2) +
    xlim(NA,0.3) + geom_treescale()
ggtree(tre, ) + geom_tiplab() + xlim(NA,0.3)

