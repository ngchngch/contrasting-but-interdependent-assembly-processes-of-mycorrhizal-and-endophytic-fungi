
##
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- "parallel";library(package = lib, character.only=TRUE);packageVersion(lib)
########
#dir <- '/media/allmember/HDD12TBA/All/Noguchi/Sugadaira'
##


#
setwd("/Volumes/Transcend/data/sugadaira_bacteria_2023/analysis")

save.dir <- "0_data_proccese/231106_OTU_convert_97"

dir.create(save.dir)

##read original functiions
source('/Volumes/Transcend/function/01_1_function.R')

#set param
vsearchpath="/Users/noguchimikihito/vsearch-2.21.1-macos-x86_64 2/bin/vsearch"

minident <- 0.97

raw_seq <- '/Volumes/Transcend/data/sugadaira_bacteria_2023/sequence_data/Linux222/merge/230920_merge_and_decontam/nonchim_seq.fasta'
###
#$Fungi
seqtab2 <- readRDS('/Volumes/Transcend/data/sugadaira_bacteria_2023/sequence_data/Linux222/merge/230920_merge_and_decontam/seqtab_rmChimera.rds')

for(l in 1:length(minident)){
  ident_th <- minident[l] 
  
  dir.create(sprintf("%s/%s",save.dir,ident_th))
  
  system2(command = vsearchpath, 
          args = c("--cluster_fast $input",raw_seq,
                   "--id",ident_th ,
                   "--mothur_shared_out", sprintf("%s/%s/ASV_OTU_corestab_%s.txt", save.dir,ident_th, ident_th),
                   "--centroids", sprintf("%s/%s/OTUseq_%s.fasta", save.dir,ident_th, ident_th),
                   "--msaout",  sprintf("%s/%s/seqAlign_%s.txt", save.dir,ident_th, ident_th) ) )
  
  otu <- read.table(sprintf("%s/%s/ASV_OTU_corestab_%s.txt", save.dir,ident_th, ident_th),
                    header=TRUE, row.names=2)[,-c(1:2)]
  
  ## ========================================== ##
  
  otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab2),
                   dimnames=list(rownames(seqtab2), colnames(otu)))
  
  for(i in 1:ncol(otu)){
    
    if( sum(otu[,i])>1 ){
      memberSeq <- rownames(otu)[which(otu[,i]>0)]
      otutab[,i] <- rowSums(seqtab2[,which(colnames(seqtab2) %in% memberSeq) ])
    }else{
      centroidSeq <- colnames(otu)[i]
      otutab[,i] <- seqtab2[, which(colnames(seqtab2) == centroidSeq) ]
    }
    
  }
  
  saveRDS(otutab, sprintf("%s/%s/seqOTUtab_th%s.rds", save.dir,ident_th,ident_th))
  
}
  
