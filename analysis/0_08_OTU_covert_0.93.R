
##
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- "parallel";library(package = lib, character.only=TRUE);packageVersion(lib)
########
#dir <- '/media/allmember/HDD12TBA/All/Noguchi/Sugadaira'
##

##
save.dir <-sprintf("%s/08_OTUtable_convert_th0.93",current_dir)

dir.create(save.dir)

#set param
vsearchpath="/Users/noguchimikihito/vsearch-2.21.1-macos-x86_64 2/bin/vsearch"

minident <- 0.93

raw_seq <- 'bioinfo_output/sample_process/fungi/OTUseq_0.97.fasta'
###read sequence files

info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_ro <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_rootF.rds")
df_root <- df_ro[[1]]

cov_th <- df_ro[[2]]

df_so <- readRDS('0_data_processing/04_rarefaction/covrarefy_sqtb_soilF.rds')
df_soil <- df_so[[1]]

cov_th2 <- df_so[[2]]


taxa <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")


taxa_f <- cbind(taxa,Genus2=NA)
taxa_f[,"Genus2"] <- apply(taxa,1,function(x){#x <- taxa_f[1,]
  if(length(grep("Unidentified",x["Genus"]))==0){sprintf("*%s*",x["Genus"])}else{
    x["Genus"]
  }})

info2 <- rowSelect(info,"target",sp_Names = "fungi")
info3 <- rowSelect(info2,"plant",sp_Names = c("Unidentified","-"),invert = T)
###3
df_root2 <- df_root[which(rownames(df_root) %in% rownames(info3)),]
df_root2 <- df_root2[,colSums(df_root2)>0]

#$Fungi
seqtab2 <- readRDS('bioinfo_output/sample_process/fungi/seqOTUtab.rds')

  ident_th <- minident 
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

  
####==covarage_based rarefaction
  
  otab2 <- otutab[which(rownames(otutab) %in% rownames(df_root2)),which(colnames(otutab) %in% colnames(df_root2))]


  cov_tab <- covrarefy(otab2[which(rownames(otab2) %in% rownames(df_root2)),
                             which(colnames(otab2) %in% colnames(df_root2))],
                       readth = 1000,
            ncore=8)
 
  saveRDS(cov_tab,sprintf("%s/covrared_tabs.rds", save.dir))

  
  cov_tab_s <- covrarefy(otutab[which(rownames(otutab) %in% rownames(df_soil)),
                               which(colnames(otutab) %in% colnames(df_soil))],
                         readth = 4000,
                         ncore=8)
  
  
  saveRDS(cov_tab_s,sprintf("%s/covrared_tabs_soil.rds", save.dir))
  