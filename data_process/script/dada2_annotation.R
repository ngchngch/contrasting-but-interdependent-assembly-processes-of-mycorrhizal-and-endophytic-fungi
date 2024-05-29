library(dada2)
#?assignTaxonomy
library(stringr)

library(seqinr)

setwd('XXXXXXXXXXXXXXXXXXXXXX')
source("function/01_1_function.R")

ref <-'sh_general_release_dynamic_all_29.11.2022.fasta'

  #'/media/allmember/HDD12TBA/All/Noguchi/function/referenceDB/UNITE_public_16.10.2022.fasta'
fg1_seq <- read.fasta('05_Merge_sequence_output/OTUseq_0.97.fasta')
tab <- readRDS('05_Merge_sequence_0utput/Fungi/seqOTUtab.rds')
svec1 <- sapply(fg1_seq,paste,collapse="")

tab1 <- tab
colnames(tab1) <- svec1[colnames(tab)]
tax1 <- assignTaxonomy((svec1),ref,multithread=TRUE)


save.dir <- "06_Taxonomy_annotation"
dir.create(save.dir)
###

## ===================================================== ##
# -- Compile table

tax_compile <- function(taxa.print,name_vec){
  if( max( sapply(strsplit(taxa.print[,"Kingdom"], "k__"), length), na.rm=TRUE ) > 1){
    taxa.print[,"Kingdom"] <- sapply(strsplit(taxa.print[,"Kingdom"], "k__"),  "[" ,2)
  }
  
  taxa.print[,"Phylum"] <- gsub("p__", "", taxa.print[,"Phylum"])
  taxa.print[,"Class"] <- gsub("c__", "", taxa.print[,"Class"])
  taxa.print[,"Order"] <- gsub("o__", "", taxa.print[,"Order"])
  taxa.print[,"Family"] <- gsub("f__", "", taxa.print[,"Family"])
  taxa.print[,"Genus"] <- gsub("g__", "", taxa.print[,"Genus"])
  taxa.print[,"Species"] <- gsub("s__", "", taxa.print[,"Species"])
  
  taxa.print[is.na(taxa.print)] <- "Unidentified"
  taxa.print <- gsub("unidentified", "Unidentified", taxa.print)
  
  
  rownames(taxa.print) <- names(name_vec)
  return(taxa.print)
} 


#####
tax12 <- tax_compile(tax1,svec1)
saveRDS(tax12,sprintf("%s/fungi_annnotaion_results_dada2.rds",save.dir))
