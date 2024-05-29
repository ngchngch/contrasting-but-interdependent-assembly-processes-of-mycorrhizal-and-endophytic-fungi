library(dada2)
#?assignTaxonomy
library(stringr)

library(seqinr)

setwd('/media/allmember/HDD12TBA/All/Noguchi/Sugadaira')
source("/media/allmember/HDD12TBA/All/Noguchi/function/01_1_function.R")

ref <-'/media/allmember/HDD12TBA/All/Noguchi/function/referenceDB/sh_general_release_dynamic_all_29.11.2022.fasta'

  #'/media/allmember/HDD12TBA/All/Noguchi/function/referenceDB/UNITE_public_16.10.2022.fasta'
fg1_seq <- read.fasta('/media/allmember/HDD12TBA/All/Noguchi/Sugadaira/analysis/230414_sample_process/fungi/OTUseq_0.97.fasta')
tab <- readRDS('/media/allmember/HDD12TBA/All/Noguchi/Sugadaira/analysis/230414_sample_process/fungi/seqOTUtab.rds')
svec1 <- sapply(fg1_seq,paste,collapse="")

tab1 <- tab
colnames(tab1) <- svec1[colnames(tab)]
tax1 <- assignTaxonomy((svec1),ref,multithread=TRUE)


save.dir <- "analysis/230414_dada2annotation"
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

#####
tab2 <- tab[,rownames(tax12)]
colnames(tab)[which(colSums(tab[,which(tax12[,1] == "Unidentified")])>10000)]
mean(rowSums(tab[,which(tax12[,1] == "Fungi")]))

df1 <- tab[,which(tax12[,1] == "Fungi")]

######

df1_th <- th_str(df1,5000,1000,str_detect(rownames(df1),pattern = "NS"))

uniq_nam <- unique(matrix(unlist(strsplit(rownames(df1_th),":")),
       ncol=2,byrow=T)[,1])

length(uniq_nam[grep("NS",uniq_nam)])
length(uniq_nam[grep("NF",uniq_nam)])
###############3

sinfo <- read.csv('/media/allmember/Transcend/data/Sugadaira_sequence/All_Merged/analysis/230308_merge_sample_data/comp_sample_info.csv',row.names=1)
  #readRDS('/media/allmember/HDD12TBA/All/Noguchi/221223_SAI_randamization/comp_sample_info.rds')

sinfo2 <- sinfo[which(!sinfo[,"plant"] %in% c("contami")),]


sinfo3 <- sinfo2[which(!sinfo2[,"sampling"] == "seedling"),]


usinfo <- sinfo[which(sinfo[,"plant"] == "Unidentified"),]

noplant <- rownames(usinfo)[grep("NF",rownames(usinfo))]

sample <- intersect(th123[,1],rownames(sinfo3))
length(sample)

nrow(th123[th123[,2] < 1500,])
##3
nosoil <- setdiff(rownames(sinfo3)[grep("NS",rownames(sinfo3))],th123[grep("NS",rownames(th123)),1])

table(sinfo3[grep("NF",rownames(sinfo3)),"plant"])

noroot <- setdiff(rownames(sinfo3)[grep("NF",rownames(sinfo3))],sample[grep("NF",sample)])

rowSums(df1[nosoil,])
###
write.csv(rbind(cbind(id=nosoil,tag="Soil"),
                       cbind(id=noroot,tag="Root"),
                       cbind(id=noplant,tag="Plant")),sprintf("%s/reseqsamples.csv",save.dir))

reseq <-rbind(cbind(id=nosoil,tag="Soil"),
           cbind(id=noroot,tag="Root"),
           cbind(id=noplant,tag="Plant"))


spos <- read.csv('/media/allmember/Transcend/data/Sugadaira_sequence/All_Merged/analysis/230321_reseq/sample_position.csv',row.names = 1)

#reseq2 <- read.csv('/media/allmember/Transcend/data/Sugadaira_sequence/All_Merged/analysis/230321_reseq/reseq_samples.csv',row.names = 2)

intersect(reseq[,1],rownames(reseq2))
cbind(reseq,spos[reseq[,1],])
write.csv(cbind(reseq,spos[reseq[,1],]), sprintf("%s/reseqsample_posiution_old.csv",save.dir))
