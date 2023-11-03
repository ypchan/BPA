library(this.path)
library(tidytree)
library(ggtree)
library(openxlsx)
library(treeio)
library(glue)
library(midrootiqtree)

# set the path of this script as the working directory
setwd(this.dir())

# read metadate table into R, the table must be end with '_taxa_table.xlsx'
# otherwise you should specify the table name
taxa_tbl <- grep('_taxa_table.xlsx$', list.files(), value = TRUE)
if (length(taxa_tbl) == 2) {
    df <- read.xlsx(taxa_tbl[2])
} else {
    df <- read.xlsx(taxa_tbl[1])
}

# convert ternary tree to binary tree and read the tree as treedata object
outgroup_label <- 'Melanconiella_hyperopta_CBS_131696'
tree_string <- readLines('05_iqtree/iqtree_ml.treefile')
binary_tree_string <- mid_root_iqtree(tree_string, outgroup_label)
tree_connection <- textConnection(binary_tree_string, open = "r")
tre <- read.iqtree(tree_connection)
close(tree_connection)

# we can adjust the width of the tree by modfying xlim(NA, width)
ggtree(tre, size=0.25) %<+% df + 
    geom_tiplab(aes(subset = (New == 1 & Type == 1), 
                    label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
                color = 'red', parse = T, size = 3) + 
    geom_tiplab(aes(subset = (New == 0 & Type == 1), 
                    label = glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
                color = 'black', parse = T, size = 3) + 
    geom_tiplab(aes(subset = (New == 0 & Type == 0), 
                    label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
                color = 'black', parse = T, size = 3) + 
    geom_tiplab(aes(subset = (New == 1 & Type == 0), 
                    label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
                color = 'red', parse = T, size = 3) + 
    
    geom_text2(aes(label= UFboot, subset = (UFboot >=95)), vjust = 0, hjust = 1, size = 2) +
    xlim(NA,0.65) + geom_treescale() 

# check topology if some branches are too short to recognize.
ggtree(tre, ) + geom_tiplab() + xlim(NA,0.3)

