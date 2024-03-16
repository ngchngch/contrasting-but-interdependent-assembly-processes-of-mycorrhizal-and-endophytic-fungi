
##read packages
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggtext';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'cowplot';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'igraph';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggnetwork';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'graphlayouts';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'SpiecEasi';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

save.dir <- sprintf("%s/03_slr_latent_jsdm_correspondence",current_dir)

dir.create(save.dir)

source("Script/conet_functions_231025.R")


colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####

taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("bioinfo_output/00_OTUtable_convert/0.97/OTUseq_0.97.fasta")

######
all_lat <- list.files(path="4_co-occurrence_network_slr/01_conet_estimate",
                      pattern = "conet_all_slr",recursive = T,
                      full.names = T)

netlist <- lapply(all_lat,readRDS)

names(netlist) <- sapply(strsplit(all_lat,"/"),function(x)x[3])

#####
##jsdm res

jsdm_otu <- list.files(path=jsdm_data_path,
                       pattern = "jsdm_sdist_result_iter100_Msamp5000",recursive = T,
                       full.names = T)

ll_lt <- lapply(lls,readRDS)

otu_nam <- lapply(jsdm_otu,function(x){
  a <- readRDS(x)
  return(a$species)
})

cov <- list.files(path=jsdm_data_path,
                  pattern = "OTU_cov_",recursive = T,
                  full.names = T)
covl <- lapply(cov,readRDS)

names(covl) <- sapply(strsplit(cov,"/"),function(x)x[3])



for(i in 1:length(netlist)){#i <-2
  ident_th <- names(netlist)[i]
  
  cur_dir <- sprintf("%s/%s",save.dir,ident_th)
  dir.create(cur_dir)
  #################
  #BIC with latent models
  bic <-R_BIC(netlist[[i]],netlist[[i]]$rank0$est$data)
  
  
  opt_mod <- which.min(bic)
  
  ###############
  #jsdm correlation
  cov2 <- cov2cor(covl[[i]])
  
  diag(cov2) <- NA
  dimnames(cov2) <- list(otu_nam[[i]],otu_nam[[i]])
  
  
  
  #slr_opt
  slr_net <- netlist[[i]]
  nlay <- slr_net[[names(opt_mod)]]
  
  
  L <- nlay$est$resid[[getOptInd(nlay)]]
  
  dimnames(L) <- list(colnames(nlay$est$data),colnames(nlay$est$data))
  L2 <- L[rownames(cov2),colnames(cov2)]
  
  diag(L2) <- NA
  
  
  df_cor <- as.data.frame(matrix(NA,nrow=nrow(cov2),ncol=3))
  colnames(df_cor) <- c("OTU","jsdm_low","P_jlow")
  
  for(k in 1:nrow(cov2)){#k <- 2
    cor_jlow <- cor.test(cov2[k,],L2[k,],method="spearman")
    
    df_cor[k,"OTU"] <- rownames(cov2)[k]
    df_cor[k,"jsdm_low"] <-cor_jlow$estimate
    df_cor[k,"P_jlow"] <-cor_jlow$p.value
    
    
  }
  
  
  df_cor$P_jlow_adj <- p.adjust(df_cor$P_jlow,method = "BH")
  
  
  df_cor2 <- df_cor
  df_cor2$Guild <- taxa_f[df_cor2$OTU,"Guild"]
  
  df_cor2$Guild2 <-factor(df_cor2$Guild,levels=c(setdiff(unique(df_cor2$Guild),
                                                         c("Unassigned","Unidentified")),
                                                 "Unassigned","Unidentified"))
  
  
  write.csv(df_cor2,sprintf("%s/jSDM_SLR_cor_correspondence.csv",cur_dir))
  
  ####
  
  mj_gld <- setdiff(names(table(df_cor2$Guild))[table(df_cor2$Guild)>5],c("Unassigned","Unidentified"))
  df_cor3 <- df_cor2[which(df_cor2$Guild %in% mj_gld),]
  df_cor3$Guild3 <- factor(df_cor3$Guild,levels = c("EcMF","AMF","Endophyte","Other_RAF"))
  # 正確法
  
  steel <- multiComp_SD(data=df_cor3,
               groups=df_cor3$Guild3,
               value = "jsdm_low")
  
  write.csv(steel$result_test,sprintf("%s/Steel_Dwass_res.csv",cur_dir))

  anno <- steel$alphabet[levels(df_cor3$Guild3),"alpha"]
    #####plot
  
  g <- ggplot(df_cor3,aes(x=Guild3,y=jsdm_low))+
    geom_boxplot(aes(fill=Guild3),outlier.color = NA,show.legend = F)+
    geom_star(starshape=15,position = position_jitter(width=0.2),
              fill="gray",alpha=0.3,size=2)+
    stat_summary(geom = 'text', label =anno, fun = max, vjust = -1,size=8)+
    theme_light()+
    theme(axis.title = element_text(size=18),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=15,angle=45,hjust=1),
          panel.background = element_rect(colour = "black",linewidth=1.1))+
    scale_starshape_manual(values=c(23,28,15))+
    labs(y="Correlation coefficieant\n[ residual correlation in jSDM vs latent variables in SLR ]",
          x="Functional group")+
    scale_fill_manual(values=colfng_gld[levels(df_cor2$Guild2),2])+
    coord_cartesian(ylim=c(0,1))
  g
  
  
  ggsave(sprintf("%s/jsdm_SLR_correspondence.pdf",cur_dir),h=10,w=6)
  
  
}
