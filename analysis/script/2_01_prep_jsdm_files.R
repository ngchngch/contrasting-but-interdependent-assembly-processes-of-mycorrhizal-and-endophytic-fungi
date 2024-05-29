
lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- "sjSDM";library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- "parallel";library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

####

##set directory
current_dir <- "2_jSDM"

dir.create(current_dir)

save.dir <- sprintf("%s/01_prep_jsdm_files",current_dir)

dir.create(save.dir)

#==========read data====#

info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_root <- readRDS("0_data_processing/04_rarefaction/OTU97/covrarefy_sqtb_rootF.rds")[[1]]
df_soil <- readRDS("0_data_processing/04_rarefaction/OTU97/covrarefy_sqtb_soilF.rds")[[1]]
taxa <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("data_process/05_Merge_sequence_output/Fungi/0.97/OTUseq_0.97.fasta")

taxa_f <- cbind(taxa,Genus2=NA)
taxa_f[,"Genus2"] <- apply(taxa,1,function(x){#x <- taxa_f[1,]
  if(length(grep("Unidentified",x["Genus"]))==0){sprintf("*%s*",x["Genus"])}else{
    x["Genus"]
  }})


No_soil <- setdiff(rownames(info)[grep("NS",rownames(info))],rownames(df_soil))

info2 <- info[which(!info[,"site"] %in% info[No_soil,"site"]),]
df <- df_root[which(rownames(df_root) %in% rownames(info2)),
              which(colnames(df_root) %in% rownames(rowSelect(taxa,"Kingdom","Fungi")))]

info3 <- info2[rownames(df),]

#####3

#sample_data(ps_sl)$plant
plant_th <- 30
pnumber <- table(info3[grep("NF",rownames(info2)),"plant"])
pnam <- setdiff(names(pnumber[pnumber>=plant_th]),"Unidentified")

info35 <- info3[which(info3$plant %in% pnam),]

#soil comm
sinfo <- info
rownames(sinfo) <- sinfo[,"ID"]
df_s <- df_soil/rowSums(df_soil)

dist_s <- vegdist(df_s,method="bray")

pco <- cmdscale(dist_s, k=10, eig=TRUE)
k <- max(which(cumsum(pco$eig/sum(pco$eig))<0.9))

res.pco <- cmdscale(dist_s, k=k+1, eig=TRUE)
score <- apply(res.pco$points,2,scale)
colnames(score) <- sprintf("Soil_PCo%s",1:(k+1))

#rownames(df_s)
res <- cbind(site=sinfo[rownames(res.pco$points),"site"],as.data.frame(score))

info4 <- merge(info35,res,by="site")


rownames(info4) <- info4[,"ID"]

df_sl <- df[which(rownames(df) %in% rownames(info4)),]

table(info4[,"plant"])

th <- 30
mat <- as.matrix(df_sl)
sum(colSums(mat>0)>th)
Occ<- mat[,colSums(mat>0)>th]

Occ[Occ>0] <- 1

Env <- info4[rownames(Occ),c(6,(ncol(info35)+1):ncol(info4))]

#

SP <- cbind(Y=info3[rownames(Occ),"y_m"]-min(info3[rownames(Occ),"y_m"]),
            X=info3[rownames(Occ),"x_m"]-min(info3[rownames(Occ),"x_m"])) # spatial coordinates (no effect on species occurences)

saveRDS(list(Occ=Occ,Env=Env,SP=SP),sprintf("%s/jsdm_files_sPC90per.rds",save.dir))

nam_otu97 <- colnames(Occ)
####
write.fasta(seq[colnames(Occ)],names = colnames(Occ),
            sprintf("%s/OTU_in_jSDM.fasta",save.dir))

write.csv(taxa[colnames(Occ),],sprintf("%s/OTU_in_jSDM.csv",save.dir))

##############################
#with soil pH

ph <- read.csv("metadata/sugadaira_soil_pH.csv",
               row.names=1,header=T)


#

info_ph <- info3[which(info3$site2 %in% rownames(ph)),]

info_ph2 <- cbind(info_ph,ph=ph[info_ph$site2,])
#####3

#sample_data(ps_sl)$plant
plant_th <- 15
pnumber <- table(info_ph2[grep("NF",rownames(info_ph2)),"plant"])
pnam <- setdiff(names(pnumber[pnumber>plant_th]),"Unidentified")

info35 <- info_ph2[which(info_ph2$plant %in% pnam),]

#soil comm
sinfo1 <- info[which(info$site2 %in% unique(info35$site2)),]

sinfo <- sinfo1[which(sinfo1$sample_type =="soil"),]

rownames(sinfo) <- sinfo[,"ID"]

df_s <- df_soil/rowSums(df_soil)

dist_s <- vegdist(df_s[which(rownames(df_s) %in% rownames(sinfo)),],method="bray")

pco <- cmdscale(dist_s, k=10, eig=TRUE)
k <- max(which(cumsum(pco$eig/sum(pco$eig))<0.9))

res.pco <- cmdscale(dist_s, k=k+1, eig=TRUE)
score <- apply(res.pco$points,2,scale)
colnames(score) <- sprintf("Soil_PCo%s",1:(k+1))

res <- cbind(site=sinfo[rownames(res.pco$points),"site"],as.data.frame(score))

info4 <- merge(info35,res,by="site")
rownames(info4) <- info4[,"ID"]


df_sl <- df[which(rownames(df) %in% rownames(info4)),]

th <- 15
mat <- as.matrix(df_sl)
sum(colSums(mat>0)>th)
Occ<- mat[,colSums(mat>0)>th]

Occ[Occ>0] <- 1

Env <- info4[rownames(Occ),c(6,12:ncol(info4))]

####
SP <- cbind(Y=info3[rownames(Occ),"y_m"]-min(info3[rownames(Occ),"y_m"]),
            X=info3[rownames(Occ),"x_m"]-min(info3[rownames(Occ),"x_m"])) # spatial coordinates (no effect on species occurences)

saveRDS(list(Occ=Occ,Env=Env,SP=SP),sprintf("%s/jsdm_files_ph.rds",save.dir))


##################################3
#OTU93
df_root93 <- readRDS("0_data_processing/04_rarefaction/OTU93/covrarefy_sqtb_rootF.rds")[[1]]
ncol(df_root93)

df_soil93 <-  readRDS("0_data_processing/04_rarefaction/OTU93/covrarefy_sqtb_soilF.rds")[[1]]

No_soil <- setdiff(rownames(info)[grep("NS",rownames(info))],rownames(df_soil93))

info2 <- info[which(!info[,"site"] %in% info[No_soil,"site"]),]

df <- df_root93[which(rownames(df_root93) %in% rownames(info2)),
              which(colnames(df_root93) %in% rownames(rowSelect(taxa,"Kingdom","Fungi")))]

info3 <- info2[rownames(df),]

#sample_data(ps_sl)$plant
plant_th <- 30
pnumber <- table(info3[grep("NF",rownames(info2)),"plant"])
pnam <- setdiff(names(pnumber[pnumber>=plant_th]),"Unidentified")

info35 <- info3[which(info3$plant %in% pnam),]

#soil comm
sinfo <- info
rownames(sinfo) <- sinfo[,"ID"]
df_s <- df_soil93/rowSums(df_soil93)

dist_s <- vegdist(df_s,method="bray")

pco <- cmdscale(dist_s, k=10, eig=TRUE)
k <- max(which(cumsum(pco$eig/sum(pco$eig))<0.9))

res.pco <- cmdscale(dist_s, k=k+1, eig=TRUE)
score <- apply(res.pco$points,2,scale)
colnames(score) <- sprintf("Soil_PCo%s",1:(k+1))

#rownames(df_s)
res <- cbind(site=sinfo[rownames(res.pco$points),"site"],as.data.frame(score))

info4 <- merge(info35,res,by="site")


rownames(info4) <- info4[,"ID"]

df_sl <- df[which(rownames(df) %in% rownames(info4)),]

table(info4[,"plant"])

th <- 30
mat <- as.matrix(df_sl)
sum(colSums(mat>0)>th)
Occ<- mat[,colSums(mat>0)>th]

Occ[Occ>0] <- 1

Env <- info4[rownames(Occ),c(6,(ncol(info35)+1):ncol(info4))]

#
SP <- cbind(Y=info3[rownames(Occ),"y_m"]-min(info3[rownames(Occ),"y_m"]),
            X=info3[rownames(Occ),"x_m"]-min(info3[rownames(Occ),"x_m"])) # spatial coordinates (no effect on species occurences)

saveRDS(list(Occ=Occ,Env=Env,SP=SP),sprintf("%s/jsdm_files_sPC90per_OTU93.rds",save.dir))
