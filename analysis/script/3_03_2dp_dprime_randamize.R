lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'bipartite';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'parallel';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "3_additional_analysis"

dir.create(current_dir)
save.dir <- sprintf("%s/03_dprime_2dp_randamize",current_dir)

dir.create(save.dir)

#==========read data====#

info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_rt <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_rootF.rds")[[1]]
taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

jsdm_files <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per.rds")

##
df_root <- df_rt[rownames(jsdm_files$Occ),colnames(jsdm_files$Occ)]


info$plant2 <- ifelse(info$plant!="Unidentified",sprintf("*%s*",info$plant),"Unidentified")
info_sl <- info4

Oc <- as.matrix(df_root)
Oc[Oc>0] <- 1

Oc2 <- Oc

or_res <- Taxa.mat(t(Oc2),info[rownames(Oc2),],"plant2")

#randamization
ncore <- 88
rand <- 10000

mat <- Oc2

#in site randamize

rand_mat <- list(NULL)
for(i in 1:rand){#i <- 1
  rand_mat[[i]] <- blockSample(mat,info_sl[rownames(mat),"site2"],rownames(mat))
}

rand_2dp <- mclapply(rand_mat,function(x){
  Taxa.mat(t(x),info_sl[rownames(x),],"plant2")[,colnames(or_res)]
},mc.cores=ncore)


saveRDS(rand_2dp,sprintf("%s/2dp_plant_randamize_insite.rds",save.dir))

z_val <- matrix(NA,nrow=nrow(or_res),ncol=ncol(or_res))
p_val <- matrix(NA,nrow=nrow(or_res),ncol=ncol(or_res))
for(i in 1:ncol(or_res)){
  for(j in 1:nrow(or_res)){
    rres <- sapply(rand_2dp,function(x){x[j,i]})
    z_val[j,i] <- (or_res[j,i]-mean(rres))/sd(rres)
    p_val[j,i] <- ifelse( z_val[j,i]>0,sum(or_res[j,i]<=rres)/rand,
                          sum(or_res[j,i]>=rres)/rand)
}
}

colnames(z_val) <- colnames(or_res)
rownames(z_val) <- rownames(or_res)
saveRDS(z_val,sprintf("%s/2dp_plant_zmat_plant.rds",save.dir))

rand_res <- cbind(gather(as.data.frame(cbind(OTU=rownames(or_res),z_val)),plant,zval,-1),
                  raw_p=gather(as.data.frame(p_val),plant,raw_p)[,2])


rand_res$raw_p <- rand_res$raw_p
rand_res$p_BH <- p.adjust(rand_res$raw_p,method = "BH")
rand_res$signif <- sapply(ifelse(is.na(rand_res$p_BH),1,rand_res$p_BH),function(x){sig(x,bin=T)})

saveRDS(rand_res,sprintf("%s/Z_val_2dp_plant.rds",save.dir))
#rand_res <- readRDS("3_additional_analysis/03_dprime_2dp_randamize/result_inlinux/result/Z_val_2dp_plant.rds")
write.csv(rand_res,sprintf("%s/Z_val_2dp_plant.csv",save.dir))
###################################################
#dprime
or_dpr_p <- dfun(t(or_res))$dprime
or_dpr_f <- dfun(or_res)$dprime

rand_dpr <- mclapply(rand_mat,function(x){
  tab <- Taxa.mat(t(x),info_sl[rownames(x),],"plant2")[,colnames(or_res)]
  dpl <- dfun(t(tab))$dprime
  dfg <- dfun(tab)$dprime
  return(list(Plant=dpl,Fungi=dfg))
},mc.cores=ncore)

saveRDS(rand_dpr,sprintf("%s/dprime_randamize_insite.rds",save.dir))

#rand_dpr <- readRDS("3_additional_analysis/03_dprime_2dp_randamize/result_inlinux/result/dprime_randamize_insite.rds")
#Fungi
z_dpr_f <- c()
p_dpr_f <- c()

  for(j in 1:length(or_dpr_f)){#j <-4
    rres <- sapply(rand_dpr,function(x){
      r <- x[["Fungi"]][names(or_dpr_f)[j]]
      return(r)})
    z_dpr_f[j] <- (or_dpr_f[j]-mean(rres))/sd(rres)
    #one-tailed test
    p_dpr_f[j] <- sum(or_dpr_f[j]<=rres)/rand
  }

names(z_dpr_f) <- names(or_dpr_f)
pf_adj <- p.adjust(p_dpr_f,method = "BH")

zres_dprf <- data.frame(name=names(or_dpr_f),
                        original=or_dpr_f,
                        zval=z_dpr_f,
                        p_BH=pf_adj,
                        signif=sapply(pf_adj,sig))

saveRDS(zres_dprf,sprintf("%s/dprime_zval_fungi.rds",save.dir))
write.csv(zres_dprf,sprintf("%s/dprime_zval_fungi.csv",save.dir))
#Plant
z_dpr_p <- c()
p_dpr_p <- c()

for(j in 1:length(or_dpr_p)){
  rres <- sapply(rand_dpr,function(x){
    r <- x[["Plant"]][names(or_dpr_p)]
    r[i]})
  z_dpr_p[i] <- (or_dpr_p[i]-mean(rres))/sd(rres)
  #one-tailed test
  p_dpr_p[i] <- sum(or_dpr_p[i]<=rres)/rand
}

names(z_dpr_p) <- names(or_dpr_p)
pp_adj <- p.adjust(p_dpr_p,method = "BH")

zres_dprp <- data.frame(name=names(or_dpr_p),
                        original=or_dpr_p,
                        zval=z_dpr_p,
                        p_BH=pp_adj,
                        signif=sapply(pp_adj,sig))

saveRDS(zres_dprp,sprintf("%s/dprime_zval_plant.rds",save.dir))
write.csv(zres_dprp,sprintf("%s/dprime_zval_plant.csv",save.dir))