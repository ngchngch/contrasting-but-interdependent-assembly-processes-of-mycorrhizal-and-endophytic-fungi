
##
library(phyloseq)
library(decontam)
library(ggplot2)
library(seqinr)
library(dada2)
########
#dir <- '/media/allmember/HDD12TBA/All/Noguchi/Sugadaira'
setwd(dir)

##

source("function/01_1_function.R")
##
save.dir <- "05_Merge_sequence_output"

dir.create(save.dir)

#set param
vsearchpath="XXXXXX=vsearch path=XXXXXX"
minident <- c(0.97,0.93)

###read sequence files
#$Fungi

seqfiles <- list.files("04_Denoising/Fungi",full.names = T)

seqlist <- lapply(seqfiles,readRDS)


#read sample infomation

info <- read.csv("Raw_data/sample_info.csv",row.names = 1)



###############################################################33

#named each seq files
names(seqlist) <- sapply(seqfiles,function(x){
  nam_vec <- strsplit(strsplit(x,"/")[[1]][4],"_")[[1]]
  paste(nam_vec[1],nam_vec[2],sep="_")
})

#merge same libray runs
runs <- unique(substr(names(seqlist),1,4))

seqlist2 <- list(NULL)

for(i in 1:length(runs)){#i <- 1
  rnam <- names(seqlist)[grep(runs[i],names(seqlist))]

  slist <- list(NULL)
  for(j in 1:length(rnam)){
    slist[[j]]<- seqlist[[rnam[j]]]
  }
  if(length(slist)>1){
    
    seqlist2[[i]] <- mergeSequenceTables(tables=slist,repeats = "sum")
  }else{
    seqlist2[[i]] <- slist[[1]]
  }
    #rownames(seqlist2[[i]]) <- paste(rownames(seqlist2[[i]]),runs[i],sep=":")
}


#detect & delete contami asv reads

decont_fg  <- lapply(seqlist2,function(x){# x <- seplist[[1]]
    Decontam_df(x,info,th=0.1)
})

#merge sequence table
#names(decont_fg)
length(decont_fg)
fg_tab1 <- mergeSequenceTables(tables=decont_fg[1:4],repeats = "sum")
fg_tab2 <- decont_fg[[5]]
rownames(fg_tab2) <- paste(rownames(fg_tab2),5,sep=":")
fg_tab3 <- decont_fg[[6]]
rownames(fg_tab3) <- paste(rownames(fg_tab3),6,sep=":")
fg_tab4 <- decont_fg[[7]]
rownames(fg_tab4) <- paste(rownames(fg_tab4),7,sep=":")

fg_tab <- mergeSequenceTables(fg_tab1,fg_tab2,fg_tab3,fg_tab4)

#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# otu clustering

  exp="Fungi"
  dir.create(sprintf("%s/%s",save.dir,exp))
  
  mseqtab_f <- removeBimeraDenovo(fg_tab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  
  seq.mat <- cbind(colnames(mseqtab_f),sprintf('X_%s', formatC(1:ncol(mseqtab_f), width = nchar(ncol(mseqtab_f)), flag = "0"))) 
  write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/%s/merge_seq.fasta", save.dir,exp) )
  
  seqtab2 <- mseqtab_f
  colnames(seqtab2) <- seq.mat[,2]
  saveRDS(seqtab2, sprintf('%s/%s/mASVseqtab.rds', save.dir,exp))
  write.csv(cbind(sample=rownames(seqtab2), seqtab2), sprintf('%s/%s/mASVseqtab.csv', save.dir,exp), row.names=FALSE)
  
  for(id in 1:length(minident)){
    mini <- minident[id]
    dir.create(sprintf("%s/%s/%s",save.dir,exp,mini))
    
    system2(command = vsearchpath, 
            args = c("--cluster_fast $input", sprintf("%s/%s/%s/merge_seq.fasta", save.dir,exp,mini),
                     "--id", mini,
                     "--mothur_shared_out", sprintf("%s/%s/%s/ASV_OTU_corestab_%s.txt", save.dir,exp, mini, mini),
                     "--centroids", sprintf("%s/%s/%s/OTUseq_%s.fasta", save.dir,exp, mini, mini),
                     "--msaout",  sprintf("%s/%s/%s/seqAlign_%s.txt", save.dir,exp, mini, mini) ) )
    
    otu <- read.table(sprintf("%s/%s/%s/ASV_OTU_corestab_%s.txt", save.dir,exp, mini, mini),
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
    
    saveRDS(otutab, sprintf("%s/%s/%s/seqOTUtab.rds", save.dir,exp, mini))
    write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/%s/%s/seqOTUtab.csv',save.dir,exp, mini), row.names=FALSE)
    
  }

  ###check
  nam <- rowSums(otutab)[which(rowSums(otutab)>5000)]
  rowSums(otutab)[grep("NS",rownames(otutab))]
  length(nam[grep("NS",names(nam))])
  ################plant############################################################
  
  seqfiles <- list.files("04_Denoising/Plant",full.names = T)
  
  seqlist <- lapply(seqfiles,readRDS)
  
  
  
  ###############################################################33
  
  #named each seq files
  names(seqlist) <- sapply(seqfiles,function(x){
    nam_vec <- strsplit(strsplit(x,"/")[[1]][4],"_")[[1]]
    paste(nam_vec[1],nam_vec[2],sep="_")
  })
  
  #merge same libray runs
  runs <- unique(substr(names(seqlist),1,4))
  
  seqlist2 <- list(NULL)
  
  for(i in 1:length(runs)){#i <- 1
    rnam <- names(seqlist)[grep(runs[i],names(seqlist))]
    
    slist <- list(NULL)
    for(j in 1:length(rnam)){
      slist[[j]]<- seqlist[[rnam[j]]]
    }
    if(length(slist)>1){
      
      seqlist2[[i]] <- mergeSequenceTables(tables=slist,repeats = "sum")
    }else{
      seqlist2[[i]] <- slist[[1]]
    }
    #rownames(seqlist2[[i]]) <- paste(rownames(seqlist2[[i]]),runs[i],sep=":")
  }
  
  
  #detect & delete contami asv reads
  
  decont_pl  <- lapply(seqlist2,function(x){# x <- seplist[[1]]
    Decontam_df(x,info,th=0.1)
  })
  
  #merge sequence table
  #names(decont_fg)
  length(decont_pl)
  pl_tab1 <- mergeSequenceTables(tables=decont_pl[1:4],repeats = "sum")
  pl_tab3 <- decont_pl[[5]]
  rownames(pl_tab3) <- paste(rownames(pl_tab3),6,sep=":")
  pl_tab4 <- decont_pl[[6]]
  rownames(pl_tab4) <- paste(rownames(pl_tab4),7,sep=":")
  
  pl_tab <- mergeSequenceTables(pl_tab1,pl_tab3,pl_tab4)
  
  #seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  # otu clustering
  
  exp="Plant"
  dir.create(sprintf("%s/%s",save.dir,exp))
  
  mseqtab_f <- removeBimeraDenovo(pl_tab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  
  seq.mat <- cbind(colnames(mseqtab_f),sprintf('X_%s', formatC(1:ncol(mseqtab_f), width = nchar(ncol(mseqtab_f)), flag = "0"))) 
  write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/%s/merge_seq.fasta", save.dir,exp) )
  
  seqtab2 <- mseqtab_f
  colnames(seqtab2) <- seq.mat[,2]
  saveRDS(seqtab2, sprintf('%s/%s/mASVseqtab.rds', save.dir,exp))
  write.csv(cbind(sample=rownames(seqtab2), seqtab2), sprintf('%s/%s/mASVseqtab.csv', save.dir,exp), row.names=FALSE)
  
  system2(command = vsearchpath, 
          args = c("--cluster_fast $input", sprintf("%s/%s/merge_seq.fasta", save.dir,exp),
                   "--id", minident,
                   "--mothur_shared_out", sprintf("%s/%s/ASV_OTU_corestab_%s.txt", save.dir,exp, minident[1]),
                   "--centroids", sprintf("%s/%s/OTUseq_%s.fasta", save.dir,exp, minident[1]),
                   "--msaout",  sprintf("%s/%s/seqAlign_%s.txt", save.dir,exp, minident[1]) ) )
  
  otu <- read.table(sprintf("%s/%s/ASV_OTU_corestab_%s.txt", save.dir,exp, minident[1]),
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
  
  saveRDS(otutab, sprintf("%s/%s/seqOTUtab.rds", save.dir,exp))
  write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/%s/seqOTUtab.csv',save.dir,exp), row.names=FALSE)
  
  
  
