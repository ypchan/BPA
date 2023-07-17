
library(this.path)
library(openxlsx)
library(ape)
library(tidyverse)

# set working directory
setwd(this.dir())

# get taxa table by capturing the suffix
taxa_tbl <- grep('_taxa_table.xlsx$', list.files('../'), value = TRUE)
if (length(taxa_tbl) == 2) {
    df <- read.xlsx(paste0('../',taxa_tbl[2]))
} else {
    df <- read.xlsx(paste0('../',taxa_tbl[1]))
}
head(df)

df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
head(df)

# prefix barcode FASTA files
species <- 'Dinemasporium'

# Markers will be downloded, marker must match the colunames.
maker_lst <- c('ITS',  'LSU')

for (marker in maker_lst) {
    outfilename <- paste(species,marker, sep='_')
    outfilename <- paste(outfilename,'fasta', sep='.')
    cat(outfilename,'\n')
    
    df_marker <- df %>% drop_na(all_of(marker))
    marker_obj <- read.GenBank(unlist(df_marker[marker]))
    names(marker_obj) <- df_marker$longLabel
    write.FASTA(marker_obj, outfilename)
}

