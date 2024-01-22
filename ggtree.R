library(this.path)
library(tidytree)
library(ggtree)
library(openxlsx)
library(treeio)
library(glue)
if ( !require(midrootiqtree) ) {
    devtools::install_github('ypchan/midrootiqtree')
} else {
    library(midrootiqtree)
    }



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

##-----------------------------IQTREE----------------------------------------------------------------------
# convert ternary tree to binary tree and read the tree as treedata object
outgroup_label <- 'Stemphylium_vesicarium_CBS_191.86'
tree_string <- readLines('05_iqtree/iqtree_ml.treefile')
binary_tree_string <- mid_root_iqtree(tree_string, outgroup_label)
cat(binary_tree_string, file = 'iqtree_ml.binary.treefile')
tree_connection <- textConnection(binary_tree_string, open = "r")
tre <- read.iqtree(tree_connection)
close(tree_connection)

# we can adjust the width of the tree by modfying xlim(NA, width)
ggtree(tre, size=0.25) %<+% df + 
  geom_tiplab(aes(subset = (New == 1 & Type == 1), 
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'red', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 1), 
                  label = glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 1 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'red', parse = T, size = 2.5) + 
  
  geom_text2(aes(label= UFboot, subset = (UFboot >=95)), vjust = 0, hjust = 1.2, size = 2) +
  xlim(NA,0.65) + geom_treescale() 

##-------------------MPboot------------------------------------------------------------------------
mpboot_tre <- read.iqtree('07_mpboot/concatenated.mpboot.treefile')


outgroup_label <- 'Stemphylium_vesicarium_CBS_191.86'
tree_string <- readLines('07_mpboot/concatenated.mpboot.treefile')
binary_tree_string <- mid_root_iqtree(tree_string, outgroup_label)
#cat(binary_tree_string, file = 'iqtree_ml.binary.treefile')
tree_connection <- textConnection(binary_tree_string, open = "r")
mpboot_tre <- read.iqtree(tree_connection)
close(tree_connection)

# check topology if some branches are too short to recognize.
ggtree(mpboot_tre,  size=0.25) + geom_tiplab()
ggtree(mpboot_tre, size=0.25) %<+% df + 
  geom_tiplab(aes(subset = (New == 1 & Type == 1), 
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'red', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 1), 
                  label = glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 1 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'red', parse = T, size = 2.5) + 
  
  geom_text2(aes(label= UFboot, subset = (UFboot >=95)), vjust = 0, hjust = 1, size = 2) +
  xlim(NA,50) + geom_treescale() # please remember to adjust the xlim


## ---------------------------------mrbayes----------------------------------------------------------------
mrbayes_tre <- read.mrbayes('06_mrbayes/run_mrbayes.nexus.con.tre')
ggtree(mrbayes_tre, size=0.25) %<+% df + 
  geom_tiplab(aes(subset = (New == 1 & Type == 1), 
                  label=glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'red', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 1), 
                  label = glue("bolditalic({Genus})~bolditalic({Epithet})~bold({Collection1})~bold({Number1})")), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 0 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'black', parse = T, size = 2.5) + 
  geom_tiplab(aes(subset = (New == 1 & Type == 0), 
                  label = paste0('italic(', Genus, ')~italic(', Epithet, ')', '~', Collection1,'~', Number1)), 
              color = 'red', parse = T, size = 2.5) + 
  
  geom_text2(aes(label= prob, subset = (prob >=0.75 & !isTip)), vjust = 0, hjust = 1, size = 2) +
  xlim(NA,0.65) + geom_treescale() 
