lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "0_data_processing"

dir.create(current_dir)

save.dir <- sprintf("%s/01_extract_pseq",current_dir)

dir.create(save.dir)

###############33
#ITS2
df_plant <- readRDS("data_process/05_Merge_sequence_output/Plant/seqOTUtab.rds") 

tax_p <- read.table("data_process/06_Taxonomy_annotaion/OTUseq_0.97.5nn.tsv")
tax_p2 <- tax_p[which(tax_p[,"Kingdom"] == "Viridiplantae"),]

info <- read.csv("Raw_data/metadata/sample_info.csv",
                 row.names = 1,header=T)
info2 <- info[which(info[,"target"] == "plant"),]

seq <- read.fasta("data_process/05_Merge_sequence_output/Plant/OTUseq_Plant_0.97.fasta")


ref_pl <- df_plant[which(rownames(df_plant) %in% rownames(info2[which(info2[,"sampling"]=="reference"),])),]



df_pl <- df_plant[,rownames(tax_p2)]
df_pl2 <- df_pl[,colSums(df_pl)>0]

nms <- c()
for(i in 1:ncol(df_pl2)){# i <- 1
  
  if(all(df_pl2[rownames(ref_pl),i]==0)){
    nms[i] <- paste(tax_p2[colnames(df_pl2)[i],"Genus"],colnames(df_pl2)[i],sep="_")
  }else{
    ref_nam <- info2[rownames(ref_pl)[df_pl2[rownames(ref_pl),i]>0],
                   "sample_type"]
    
    nms[i] <- paste("ref",
                    str_sub(as.vector(ref_nam),start=1,end=3),tax_p2[colnames(df_pl2)[i],
                                                                        "Genus"],
                    colnames(df_pl2)[i],sep="_")
  }
  
}


true_seq <- seq[colnames(df_pl2)]

write.fasta(true_seq,nms,sprintf("%s/plant_seq_for_tree.fasta",save.dir))

###########3
