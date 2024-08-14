library(this.path)
library(Biostrings)

# check where are you in
this.dir()
# set working directory

setwd('')

# input trimmed multiple sequence alignments
msa1 <- readDNAMultipleAlignment('Lachnum/03_triaml/ITS.mafft.trimal.fna')
msa2 <- readDNAMultipleAlignment('Lachnum/03_triaml/ITS.mafft.trimal.fna')
msa3 <- readDNAMultipleAlignment('Lachnum/03_triaml/ITS.mafft.trimal.fna')

# add place holders
msa_lst <- list(msa1, msa2, msa3)
seq_identifiers <- c()

for ( i in seq_along(msa_lst) ) {
  msa <- msa_lst[[i]]
  rownames(msa) <- gsub("_R_", "", rownames(msa))
  seq_identifiers <- c(seq_identifiers, rownames(msa))
  msa_lst[[i]] <- msa
}
unique_seq_identifiers <- unique(seq_identifiers)

unique_seq_identifiers
msa_lst[[1]]

for ( i in seq_along(msa_lst) ) {
  msa <- msa_lst[[i]]
  
  for ( identifier in unique_seq_identifiers ) {
    if ( !(identifier %in% rownames(msa)) ) {
      placeholder <- DNAStringSet(strrep("-", nchar(msa)))
      names(placeholder) <- identifier
      msa <- replaceAt(msa, at = identifier, value = placeholder)
    }
  }
  msa_lst[[i]] <- msa
}

for ( identifier in unique_seq_identifiers ) {
  concatenated_seq <- ""
  for ( msa in msa_lst ) {
    concatenated_seq <- paste0(concatenated_seq, msa@unmasked[identifier])
  }
  cat('>', identifier, '\n', sep = "", file = 'cconcatenated.fna', append = TRUE)
  cat(concatenated_seq, '\n', file = 'cconcatenated.fna', append = TRUE)
}
