
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

save.dir <- sprintf("%s/03_conet_slr_jsdm_correspondence",current_dir)

dir.create(save.dir)

source("Script/conet_functions_231025.R")


colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####

taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("bioinfo_output/sample_process/fungi/OTUseq_0.97.fasta")

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



for(i in 1:length(netlist)){#i <- 1
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
  
  #optimal inverse covariance -> covariance
  slr.icov <- solve(as.matrix(nlay$est$icov[[getOptInd(nlay)]])) 
  #head(slr.icov)
  #covariance -> correlation
  slr.cor <- cov2cor(slr.icov)
  
  #optimal network
  slr.refit <- as.matrix(getRefit(nlay))
  slr.cor[which(slr.refit==FALSE)] <- 0
  dimnames(slr.cor) <- list(colnames(nlay$est$data),colnames(nlay$est$data))
  
  cor2 <- slr.cor[rownames(cov2),]
  
  diag(cor2) <- NA
  
  L <- nlay$est$resid[[getOptInd(nlay)]]
  
  dimnames(L) <- list(colnames(nlay$est$data),colnames(nlay$est$data))
  L2 <- L[rownames(cov2),]
  
  diag(L2) <- NA
  
  
  df_cor <- as.data.frame(matrix(NA,nrow=nrow(cov2),ncol=5))
  colnames(df_cor) <- c("OTU","jsdm_sp","P_jsp","jsdm_low","P_jlow")
  
  for(k in 1:nrow(cov2)){#k <- 2
    cor_jsp <- cor.test(cov2[k,],cor2[k,],method="spearman")
    cor_jlow <- cor.test(cov2[k,],L2[k,],method="spearman")
    
    df_cor[k,"OTU"] <- rownames(cov2)[k]
    df_cor[k,"jsdm_sp"] <-cor_jsp$estimate
    df_cor[k,"P_jsp"] <-cor_jsp$p.value
    df_cor[k,"jsdm_low"] <-cor_jlow$estimate
    df_cor[k,"P_jlow"] <-cor_jlow$p.value
    
    
  }
  
  
  df_cor$P_jsp_adj <- p.adjust(df_cor$P_jsp,method = "BH")
  df_cor$P_jlow_adj <- p.adjust(df_cor$P_jlow,method = "BH")
  
  ##delete row cntaining NA (no correlaton in SLR)
  df_cor2 <- na.omit(df_cor)
  
  df_cor2$sig_jsp <- sapply(df_cor2$P_jsp_adj,function(x)sig(x,star = F))
  
  df_cor2$signif <- gsub("p","*P* ",df_cor2$sig_jsp)
  
  
  df_cor2$signif <- gsub("<","< ",df_cor2$signif)
  
  df_cor2[which(df_cor2$signif == ""),"signif"] <- "N.S." 
  
  
  df_cor2$Guild <- taxa_f[df_cor2$OTU,"Guild"]
  
  df_cor2$Guild2 <-factor(df_cor2$Guild,levels=c(setdiff(unique(df_cor2$Guild),
                                                         c("Unassigned","Unidentified")),
                                                 "Unassigned","Unidentified"))
  
  
  write.csv(df_cor2,sprintf("%s/jSDM_SLR_cor_correspondence.csv",cur_dir))
  
  ####
  
  mj_gld <- setdiff(names(table(df_cor2$Guild))[table(df_cor2$Guild)>5],c("Unassigned","Unidentified"))
  df_cor3 <- df_cor2[which(df_cor2$Guild %in% mj_gld),]
  
  # 正確法
  tab <- data.frame(mj_gld,
                    c(1:length(mj_gld)))
  rownames(tab) <- tab[,1]
  
  set.seed(0)
  p_exact = NSM3::pSDCFlig(df_cor3$jsdm_sp,
                           df_cor3$Guild2)
  
  p_exact2 = NSM3::pSDCFlig(df_cor3$jsdm_low,
                            df_cor3$Guild2)
  
  
  res <- data.frame(label=p_exact$labels,
                    Sparse=p_exact$p.val,
                    latent=p_exact2$p.val)
  
  write.csv(res,sprintf("%s/Steel_Dwass_res.csv",cur_dir))
  #####plot
  
  g <- ggplot(df_cor2,aes(x=jsdm_low,y=jsdm_sp))+
    geom_star(aes(fill=Guild,starshape=signif),size=4,show.legend = F)+
    theme_light()+
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=15),
          panel.background = element_rect(colour = "black",linewidth=1.1))+
    scale_starshape_manual(values=c(23,28,15))+
    scale_fill_manual(values=colfng_gld[levels(df_cor2$Guild2),2])
  g
  
  
  gdens1 <- ggplot(df_cor3,aes(x=jsdm_low))+
    geom_density(aes(group=Guild),fill="transparent",color="black",
                 linewidth=3,show.legend = F)+
    geom_density(aes(fill=Guild,color=Guild),
                 linewidth=1.5,alpha=0.3,show.legend = F)+
    theme_classic()+
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=15),
          panel.background = element_rect(colour = "black",linewidth=1.1))+
    scale_fill_manual(values=colfng_gld[levels(df_cor3$Guild2),2])+
    scale_color_manual(values=colfng_gld[levels(df_cor3$Guild2),2])
  
  
  gdens1
  
  gdens2 <- ggplot(df_cor3,aes(x=jsdm_sp))+
    geom_density(aes(group=Guild),fill="transparent",color="black",
                 linewidth=3,show.legend = F)+
    geom_density(aes(fill=Guild,color=Guild),
                 linewidth=1.5,alpha=0.3,show.legend = F)+
    theme_classic()+
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=15),
          panel.background = element_rect(colour = "black",linewidth=1.1))+
    scale_fill_manual(values=colfng_gld[levels(df_cor3$Guild2),2])+
    scale_color_manual(values=colfng_gld[levels(df_cor3$Guild2),2])+
    coord_flip()
  
  gdens2
  
  
  
  plot_grid(NULL,NULL,NULL,
            gdens1+theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
            NULL,NULL,
            g+labs(y="Correlation coefficieant\n[ residual correlation in jSDM vs Sparse correlation in SLR ]",
                   x="Correlation coefficieant\n[ residual correlation in jSDM vs latent variables in SLR ]"),
            gdens2+theme(axis.title.y = element_blank(),axis.text.y = element_blank()),NULL,
            nrow=3,align = "hv",
            rel_widths = c(.8,.3,.05),rel_heights = c(.1,.3,.8),
            labels = c(NA,NA,NA,"a",NA,NA,"b","c",NA),label_size = 24,vjust = -0.3)
  
  
  ggsave(sprintf("%s/jsdm_SLR_correspondence.pdf",cur_dir),h=12,w=11.5)
  
  ggsave(plot=g_legend(g+geom_star(aes(fill=Guild,starshape=signif),size=4,show.legend = T)+
                         labs(starshape="Significance of correlation\n[ jSDM vs Sparse correlation ]",
                              fill="Functional group")+
                         theme(legend.title = element_markdown(),
                               legend.text = element_markdown())+
                         guides(fill=guide_legend(override.aes = list(starshape=15)))),
         sprintf("%s/legend_jsdm_SLR_correspondence.pdf",cur_dir),h=5,w=8)
  
  
}
