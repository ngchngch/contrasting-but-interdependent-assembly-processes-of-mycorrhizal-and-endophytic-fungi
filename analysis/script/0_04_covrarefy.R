##read packages
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'parallel';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'dplyr';library(package = lib, character.only=TRUE);packageVersion(lib)


####

##set directory
current_dir <- "0_data_processing"

dir.create(current_dir)

save.dir <- sprintf("%s/04_rarefaction",current_dir)

dir.create(save.dir)

#==========read data====#

rtl <- list.files(sprintf("%s/03_make_merge_seqtab",current_dir),pattern = "No_rarefy_sqtb_rootF",
                  recursive = T,full.names = T)

root_list <- lapply(rtl,readRDS)                  

##change here
names(root_list) <- sapply(strsplit(rtl,"/"),function(x){x[3]})

sll <- list.files(sprintf("%s/03_make_merge_seqtab",current_dir),pattern = "No_rarefy_sqtb_soilF",
                  recursive = T,full.names = T)

soil_list <- lapply(sll,readRDS)

##change here
names(soil_list) <- sapply(strsplit(sll,"/"),function(x){x[3]})

#==========rarefaction=================#


covrarefy <- function(df,readth=0,ncore){
  
  OTU_table <- df[rowSums(df)>readth,]
  s_comp <- list(NULL)
  for(i in 1:nrow(OTU_table)){
    s_comp[[i]]<-OTU_table[i,]
  }
  rareslopelist<-mclapply(s_comp,function(x){
    rareslope(x,1:(sum(x)-1))
  },mc.cores=ncore)
  
  getmincov<-unlist(mclapply(rareslopelist,function(x){
    x[length(x)]
  },mc.cores=ncore))
  
  #histogram(getmincov)
  cov_th <- max(getmincov)
  #histogram(getmincov[getmincov<=cov_th])
  
  #pass_sample <- rownames(OTU_table)[getmincov<=cov_th]
  
  #length(pass_sample)/nrow(OTU_table)
  
  
  #指定したカバレッジに到達した（＝傾きが指定値を下回る）瞬間のリード数をサンプルごとに採ってくる
  cvrfun<-function(x){min(which(x<=max(getmincov[getmincov<=cov_th])))+1} #関数を設定。上記で1を引いた分を足し戻す
  cvrrare<-unlist(mclapply(rareslopelist,cvrfun,mc.cores = ncore))　#lapply+unlistでベクトル形式にして一括で値を取得
  
  cvrrare2 <- cvrrare[getmincov<=cov_th]
  set.seed(123) #再現性をとるためにランダム変数を固定（数字は何でもいい）
  OTU_covrared<-rrarefy(OTU_table[getmincov<=cov_th,],cvrrare2) #得られたリード数に沿って各サンプルからリサンプリング
  
  
  OTU_covrared2 <- OTU_covrared[,colSums(OTU_covrared)>0]
  return(list(table=OTU_covrared2,cov_th=(1-cov_th)))
}

################

for(i in 1:length(rtl)){#i <- 1
  ident_th <- names(root_list)[i]
  
  new_dir <- sprintf("%s/%s",save.dir,ident_th)
  dir.create(new_dir)
  
  root <- root_list[[ident_th]]
  soil <- soil_list[[ident_th]]
  
  #===root=========#
  rdf_root <- covrarefy(root,ncore=8)
  
  saveRDS(rdf_root,sprintf("%s/covrarefy_sqtb_rootF.rds",new_dir))
  #============soil=====================#
  rdf_soil <- covrarefy(soil,ncore=8)
  
  saveRDS(rdf_soil,sprintf("%s/covrarefy_sqtb_soilF.rds",new_dir))
  
}
