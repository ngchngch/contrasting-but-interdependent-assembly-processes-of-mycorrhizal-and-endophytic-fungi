lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

save.dir <- sprintf("%s/03_make_merge_seqtab",current_dir)

dir.create(save.dir)

###############33
#OTU97

tab_list <- list.files("bioinfo_output/00_OTUtable_convert",
                       pattern = "seqOTUtab_th",
                       recursive = T,full.names = T)

tab <- lapply(tab_list,readRDS)
  
names(tab) <- sapply(strsplit(tab_list,"_"),function(x){
  paste("OTU",str_sub(x[6],5,6),sep="")
})

tax <- readRDS("bioinfo_output/230414_dada2annotation/fungi_annnotaion_results_dada2.rds")

tax_f <- tax[which(tax[,"Kingdom"] == "Fungi"),]

tab_f_list <-lapply(tab,function(x){x[,which(colnames(x) %in% rownames(tax_f))]}) 


for(i in 1:length(tab_list)){
  new_dir <- sprintf("%s/%s",save.dir,names(tab_f_list)[i])
  dir.create(new_dir)
  
  tab_f <- tab_f_list[[i]]
  #select samples by read number
  enh_sam <- th_str(tab_f,5000,1000,str_detect(rownames(tab_f),pattern="NS"))
  
  #separate sample by run
  run1 <- enh_sam[grep(":",enh_sam[,1],invert=T),1]
  run2 <- sapply(strsplit(enh_sam[grep(":5",enh_sam[,1]),1],":"),
                 function(x){x[1]})
  run3 <- sapply(strsplit(enh_sam[grep(":6",enh_sam[,1]),1],":"),
                 function(x){x[1]})
  run4 <- sapply(strsplit(enh_sam[grep(":7",enh_sam[,1]),1],":"),
                 function(x){x[1]})
  
  #merge runs
  run2_sample <- setdiff(run2,run1)
  run3_sample <- setdiff(run3,unique(c(run1,run2)))
  run4_sample <- setdiff(run4,unique(c(run1,run2,run3)))
  
  tab_sl <- tab_f[c(run1,
                    paste(run2_sample,5,sep=":"),
                    paste(run3_sample,6,sep=":"),
                    paste(run4_sample,7,sep=":")),]
  
  nrow(tab_sl)
  
  write.csv(cbind(sample_name_withRun=c(run1[grep("NS|NF",run1)],
                                        paste(run2_sample,5,sep=":"),
                                        paste(run3_sample,6,sep=":"),
                                        paste(run4_sample,7,sep=":")),
                  tru_name=c(run1[grep("NS|NF",run1)],run2_sample,run3_sample,run4_sample)),
            sprintf("%s/sample_list.csv",new_dir))
  
  if(all(table(sapply(strsplit(rownames(tab_sl),":"),
                      function(x){x[1]}))==1)){
    rownames(tab_sl) <- sapply(strsplit(rownames(tab_sl),":"),
                               function(x){x[1]})
  }
  
  
  #sample(rownames(tab_sl))
  ##########separate root/soil sample
  
  soil_tab <- tab_sl[grep("NS",rownames(tab_sl)),]
  
  root_tab <- tab_sl[grep("NF",rownames(tab_sl)),]
  
  
  saveRDS(soil_tab,sprintf("%s/No_rarefy_sqtb_soilF.rds",new_dir))
  
  saveRDS(root_tab,sprintf("%s/No_rarefy_sqtb_rootF.rds",new_dir))
  
  ######
  
}
