##read packages
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'parallel';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'dplyr';library(package = lib, character.only=TRUE);packageVersion(lib)


##set directory
current_dir <- "0_data_processing"

dir.create(current_dir)

save.dir <- sprintf("%s/02_root_plant_annotation",current_dir)

dir.create(save.dir)

##read table
info <- read.csv("Raw_data/metadata/sample_info.csv")
rownames(info) <- info[,"ID"]

#read_sampleanamelist
snam_l <- read.csv("Raw_data/metadata/sample_list.csv",
                      row.names=1,header=T)
snam_list <- snam_l[grep("NF",snam_l[,1]),]
rownames(snam_list) <- snam_list[,1]

##add Manual BLAST search info to taxa_fp1

taxa_fp2 <- read.csv(sprintf("%s/01_extract_pseq/Manual/Plant_seq_in_ITS1_BLAST.csv",current_dir),
                     row.names=1,header=T)


df_f <- readRDS("data_process/05_Merge_sequence_output/Fungi/0.97/seqOTUtab.rds")
####

P_annotaion <- function(df_f,taxa_fp2){
  unlist(apply(df_f[,rownames(taxa_fp2)],1,function(x){#x <- df_f[1,rownames(taxa_fp)]
    
    if(length(unique(taxa_fp2[names(x[x>0]),"Genus_5nn_BLAST"]))==1){
      return(unique(taxa_fp2[names(x[x>0]),"Genus_5nn_BLAST"]))
    }else{
      if(length(unique(taxa_fp2[names(x[x>0]),"Genus_5nn_BLAST"]))==0){
        return("Unidentified")
      }else{
        
        return(paste("contami",paste(unique(taxa_fp2[names(x[x>0]),"Genus_5nn_BLAST"]),collapse="_"),sep="_"))
      }
    }
    
  }))
  
} 
#

plant_ann_f <- P_annotaion(df_f,taxa_fp2)


plant_ann_f2 <- plant_ann_f[rownames(snam_list)]
plant_ann_f2[grep("contami",plant_ann_f2)] <- "contami"

names(plant_ann_f2) <- snam_list[names(plant_ann_f2),2]


#setdiff(names(plant_ann_f2),rownames(snam_list))
#==plant===#
##ITS_3pl62F1 
df_p <- readRDS("data_process/05_Merge_sequence_output/Plant/seqOTUtab.rds")

#Plant annotation result with Manual BLAST search
taxa_p3 <- read.csv("0_data_process/01_extract_pseq/Manual/Plant_annotaion_Genus.csv",row.names=1,header=T)

df_p2 <- df_p[,which(colnames(df_p) %in% rownames(taxa_p3))]


#separate sample by run
run1 <- rownames(df_p2[grep(":",rownames(df_p2),invert=T),])

run3 <- sapply(strsplit(rownames(df_p2[grep(":6",rownames(df_p2)),]),":"),
               function(x){x[1]})
run4 <- sapply(strsplit(rownames(df_p2[grep(":7",rownames(df_p2)),]),":"),
               function(x){x[1]})


txtab_p <- Taxa.mat(df_p2[run1,],
                    taxa_p3,"Genus_5nn_BLAST")

df_p3 <- df_p[run1,]

plant_ann <- P_annotaion(df_p3,taxa_p3)
  


plant_ann2 <- plant_ann
plant_ann2[grep("contami",plant_ann2)] <- "contami"


unident_run1 <- names(plant_ann2[grep("Unidentified",plant_ann2)])

#####run3

df_p4 <- df_p[paste(intersect(run3,unident_run1),":6",sep=""),]

plant_ann_r3 <- P_annotaion(df_p4,taxa_p3)

plant_ann2_r3 <- plant_ann_r3
plant_ann2_r3[grep("contami",plant_ann2_r3)] <- "contami"

table(plant_ann2_r3)

plant_ann2_2 <- plant_ann2
new_ident_sample <- gsub(":6","",names(plant_ann2_r3)[grep("Unidentified",plant_ann2_r3,invert = T)])
plant_ann2_2[new_ident_sample] <- plant_ann2_r3[paste(new_ident_sample,6,sep=":")]

unident_run13 <- names(plant_ann2_2[grep("Unidentified",plant_ann2_2)])

#####run4
# noadditional sequenceresult in run4

#===Merge 2 region result=====#
pann_f <- plant_ann_f2
names(pann_f) <- gsub("NF","NP",names(plant_ann_f2))


m_rslt <-unique(merge(cbind(ID=names(plant_ann2_2),plant_ann2_2),
                      cbind(ID=names(pann_f),pann_f),by="ID",all=T))
m_rslt[is.na(m_rslt[,3]),3] <- "Unidentified"
m_rslt[is.na(m_rslt[,2]),2] <- "Unidentified"

rownames(m_rslt) <- m_rslt[,1]
merge_pann <- apply(m_rslt,1,function(x){
  if(length(grep("Unidentified",x[2]))==0){
    return(x[2])
  }else{
    if(length(grep("Unidentified",x[3]))==0){
      return(x[3])
    }else{
      return(x[2])
    }
  }
})
table(merge_pann)

p_ann_tab <- cbind(ID=gsub("NP","",names(merge_pann)),merge_pann)


info2 <- cbind(info[which(info[,"sample_type"] == "root" &info[,"sampling"] == "root"),],plant="Unidentified")

for(i in 1:nrow(p_ann_tab)){
  info2[grep(p_ann_tab[i,1],rownames(info2)),"plant"] <- p_ann_tab[i,2]
}

info3 <- rbind(info2,cbind(info[setdiff(rownames(info),rownames(info2)),],plant="-"))

sum(table(info2[which(info2[,"target"] == "fungi"),"plant"]))


saveRDS(info3,sprintf("%s/sample_info_with_plant.rds",save.dir))
write.csv(info3,sprintf("%s/sample_info_with_plant.csv",save.dir))

