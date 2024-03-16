
##read packages
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggtext';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggforce';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'igraph';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggnetwork';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'graphlayouts';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'SpiecEasi';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

save.dir <- sprintf("%s/02_conet_analysis_3method",current_dir)

dir.create(save.dir)

source("Script/conet_functions_231025.R")


colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####

taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("bioinfo_output/sample_process/fungi/OTUseq_0.97.fasta")

######
all_lat <- list.files(path="4_co-occurrence_network/01_conet_estimate",
                      pattern = "conet_all_slr",recursive = T,
                      full.names = T)

netlist <- lapply(all_lat,readRDS)

names(netlist) <- sapply(strsplit(all_lat,"/"),function(x)x[11])

#optimal slr networks

slr <- list.files(path="4_co-occurrence_network/01_conet_estimate",
                  pattern = "conet_opt_slr",recursive = T,
                  full.names = T)

slrlist <- lapply(slr,readRDS)

names(slrlist) <- sapply(strsplit(slr,"/"),function(x)x[3])

#####
occtab <- list.files(path=conet_data_path,
                     pattern = "jsdm_occtable",recursive = T,
                     full.names = T)


tablist <- lapply(occtab,readRDS)

namlist <- lapply(tablist,colnames)

#mb network

mb <- list.files(path="4_co-occurrence_network/01_conet_estimate",
                 pattern = "conet_mb",recursive = T,
                 full.names = T)

mblist <- lapply(mb,readRDS)

names(mblist) <- sapply(strsplit(mb,"/"),function(x)x[3])

#glasso network


glasso <- list.files(path="4_co-occurrence_network/01_conet_estimate",
                     pattern = "conet_glasso",recursive = T,
                     full.names = T)

glassolist <- lapply(glasso,readRDS)

names(glassolist) <- sapply(strsplit(glasso,"/"),function(x)x[3])

#raw_read table
rtab <- list.files(path=conet_data_path,
                   pattern = "conet_input",recursive = T,
                   full.names = T)

rtablist <- lapply(rtab,readRDS)

names(rtablist) <- sapply(strsplit(rtab,"/"),function(x)x[4])

##jsdm res
lls <- list.files(path=jsdm_data_path,
                  pattern = "restricted_lls_mat_",recursive = T,
                  full.names = T)

jsdm_otu <- list.files(path=jsdm_data_path,
                       pattern = "jsdm_sdist_result_iter100_Msamp5000",recursive = T,
                       full.names = T)

ll_lt <- lapply(lls,readRDS)

otu_nam <- lapply(jsdm_otu,function(x){
  a <- readRDS(x)
  return(a$species)
})

ll_list <- list(NULL)
for(i in 1:length(ll_lt)){
  x <- ll_lt[[i]]
  y <- otu_nam[[i]]
  ll_list[[i]] <- data.frame(OTU=y,sapply(x,colSums))
}



ll_b2 <- lapply(ll_list,function(x){
  a <- cbind(llr=x$NoBio-x$Full,taxa_f[x$OTU,])
  rownames(a) <- x$OTU
  
  return(a)
})

names(ll_b2) <- sapply(strsplit(lls,"/"),function(x)x[11])

##########

for(i in 1:length(rtablist)){
  show.progress(i,1:length(rtablist))
  #i <-1
  ident_th <- names(rtablist)[i]
  ll <- ll_b2[[i]]
  
  cur_dir <- sprintf("%s/%s",save.dir,ident_th)
  dir.create(cur_dir)
  #################
  #BIC with latent models
  bic <-R_BIC(netlist[[i]],rtablist[[i]])
  
  bic_models <- drow_ModsBIC(bic)
  
  ggsave(plot=bic_models,sprintf("%s/BIC_Conet_models_%s.pdf",cur_dir,ident_th),h=7,w=8)
  
  sink(sprintf("%s/BIC_Conet_models_%s.txt",cur_dir,ident_th))
  print(bic)
  sink()
  ####################
  
  opt_mod <- which.min(bic)
  
  ###############
  #edge weight
  #slr_opt
  ew <- Edge_weight(net=slrlist[[i]],zero=T,method="slr")
  
  
  ew1 <- ew[which(taxa_f[ew$OTU1,"Guild"] != "AMF" | taxa_f[ew$OTU2,"Guild"] != "AMF"),]
  ew_p <- ew1[order(ew1$cor)[1:20],]
  
  ew_n <- ew1[order(ew1$cor,decreasing=T)[1:20],]
  
  write.csv(cbind(ew_p,
                  OTU1=taxa_f[ew_p$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_p$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/SLR_Top20_correlation_pairs_P_%s_%s.csv",
                    cur_dir,ident_th,names(opt_mod)))
  write.csv(cbind(ew_n,
                  OTU1=taxa_f[ew_n$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_n$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/SLR_Top20_correlation_pairs_N_%s_%s.csv",
                    cur_dir,ident_th,names(opt_mod)))
  
  #MB
  ew_mb <- Edge_weight(net=mblist[[i]],zero=T,method="mb")
  
  
  ew_mb1 <- ew_mb[which(taxa_f[ew_mb$OTU1,"Guild"] != "AMF" | taxa_f[ew_mb$OTU2,"Guild"] != "AMF"),]
  ew_mb_p <- ew_mb1[order(ew_mb1$cor)[1:20],]
  
  ew_mb_n <- ew_mb1[order(ew_mb1$cor,decreasing=T)[1:20],]
  
  write.csv(cbind(ew_mb_p,
                  OTU1=taxa_f[ew_mb_p$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_mb_p$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/MB_Top20_correlation_pairs_P_%s.csv",
                    cur_dir,ident_th))
  write.csv(cbind(ew_mb_n,
                  OTU1=taxa_f[ew_mb_n$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_mb_n$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/MB_Top20_correlation_pairs_N_%s.csv",
                    cur_dir,ident_th))
  #glasso
  ew_gl <- Edge_weight(net=glassolist[[i]],zero=T,method="glasso")
  
  
  ew_gl1 <- ew_gl[which(taxa_f[ew_gl$OTU1,"Guild"] != "AMF" | taxa_f[ew_gl$OTU2,"Guild"] != "AMF"),]
  ew_gl_p <- ew_gl1[order(ew_gl1$cor)[1:20],]
  
  ew_gl_n <- ew_gl1[order(ew_gl1$cor,decreasing=T)[1:20],]
  
  write.csv(cbind(ew_gl_p,
                  OTU1=taxa_f[ew_gl_p$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_gl_p$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/glasso_Top20_correlation_pairs_P_%s.csv",
                    cur_dir,ident_th))
  write.csv(cbind(ew_gl_n,
                  OTU1=taxa_f[ew_gl_n$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_gl_n$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/glasso_Top20_correlation_pairs_N_%s.csv",
                    cur_dir,ident_th))
  ######################################
  
  hic_otu <- unique(c(ew_p$OTU1,ew_p$OTU2,
                      ew_n$OTU1,ew_n$OTU2,
                      ew_mb_p$OTU1,ew_mb_p$OTU2,
                      ew_mb_n$OTU1,ew_mb_n$OTU2,
                      ew_gl_p$OTU1,ew_gl_p$OTU2,
                      ew_gl_n$OTU1,ew_gl_n$OTU2))
  
  write.fasta(seq[hic_otu],names=hic_otu,
              sprintf("%s/high_correlattion_pairs_seq_%s.fasta",save.dir,ident_th))
  
  ################
  #drow co-net
  #slr
  nlay <- slrlist[[i]]
  #optimal inverse covariance -> covariance
  slr.icov <- solve(as.matrix(nlay$est$icov[[getOptInd(nlay)]])) 
  #head(slr.icov)
  #covariance -> correlation
  slr.cor <- cov2cor(slr.icov)
  
  #optimal network
  slr.refit <- as.matrix(getRefit(nlay))
  
  #extract edge weight of optimal network
  slr.cor[which(slr.refit==FALSE)] <- 0
  diag(slr.cor) <- 0
  slr <- abs(slr.cor)
  
  #mb
  nlay2 <- mblist[[i]]
  #
  slr.cor <- as.matrix(symBeta(getOptBeta(nlay2), mode='ave'))
  
  diag(slr.cor) <- 0
  lay_mb <- abs(slr.cor)
  
  #glasso
  nlay3 <- glassolist[[i]]
  #optimal inverse covariance -> covariance
  slr.icov <- solve(as.matrix(nlay3$est$icov[[getOptInd(nlay3)]])) 
  #head(slr.icov)
  #covariance -> correlation
  slr.cor <- cov2cor(slr.icov)
  
  #optimal network
  slr.refit <- as.matrix(getRefit(nlay))
  
  #extract edge weight of optimal network
  slr.cor[which(slr.refit==FALSE)] <- 0
  diag(slr.cor) <- 0
  lay_gl <- abs(slr.cor)
  
  netw <- igraph::graph.adjacency(slr+lay_mb+lay_gl, mode='undirected', diag=FALSE, weighted=TRUE)
  
  am.coord <-layout_nicely(netw)
  rownames(am.coord) <- colnames(nlay$est$data)
  
  slr_p <- drowConet(net=slrlist[[i]],comm=nlay$est$data,NP="posi",col_list = colfng,method="slr",
                     ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  mb_p <- drowConet(net=mblist[[i]],comm=nlay$est$data,NP="posi",col_list = colfng,method="mb",
                    ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  gl_p <- drowConet(net=glassolist[[i]],comm=nlay$est$data,NP="posi",col_list = colfng,method="glasso",
                    ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  
  slr_n <- drowConet(net=slrlist[[i]],comm=nlay$est$data,NP="nega",col_list = colfng,method="slr",
                     ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  mb_n <- drowConet(net=mblist[[i]],comm=nlay$est$data,NP="nega",col_list = colfng,method="mb",
                    ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  gl_n <- drowConet(net=glassolist[[i]],comm=nlay$est$data,NP="nega",col_list = colfng,method="glasso",
                    ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  
  plot_grid(slr_p+theme(legend.position = "none"),slr_n+theme(legend.position = "none"),
            mb_p+theme(legend.position = "none"),mb_n+theme(legend.position = "none"),
            gl_p+theme(legend.position = "none"),gl_n+theme(legend.position = "none"),
            align = "hv",nrow=3,
            labels = c("a","b","c","d","e","f"),
            label_size = 24)
  
  ggsave(sprintf("%s/merge_conet.pdf",cur_dir),h=15,w=10)
  ###########3#####
  
  
  #########
  #compare llr conet degree
  cur_dir2 <- sprintf("%s/%s",cur_dir,"Compare_llr_conet")
  dir.create(cur_dir2)
  
  llCnet_mb_p <- Comp_llCnet(slr=mblist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="posi",method="mb")
  ggsave(plot=llCnet_mb_p[["net"]],sprintf("%s/MB_posinet_deg_llr_Nolatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/MB_kendtest_posinet_deg_llr_Nolatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_mb_p[["cor"]])
  sink()
  
  llCnet_gl_p <- Comp_llCnet(slr=glassolist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="posi",method="glasso")
  ggsave(plot=llCnet_gl_p[["net"]],sprintf("%s/glasso_posinet_deg_llr_Nolatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/glasso_kendtest_posinet_deg_llr_Nolatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_gl_p[["cor"]])
  sink()
  
  llCnet_pOpt <- Comp_llCnet(slr=slrlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="posi",R=names(opt_mod),method="slr")
  ggsave(plot=llCnet_pOpt[["net"]],sprintf("%s/SLRposinet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/SLR_kendtest_posinet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  
  print(llCnet_pOpt[["cor"]])
  sink()
  #####
  
  plot_grid(llCnet_pOpt$net,
            llCnet_mb_p$net,
            llCnet_gl_p$net,
            align = "hv",nrow=1,
            labels = c("a","b","c"),
            label_size = 24)
  
  ggsave(sprintf("%s/conet_Posidegree_vs_llr.pdf",cur_dir2),h=6,w=18)
  
  ###########
  
  llCnet_nOpt <- Comp_llCnet(slr=slrlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="nega",R=names(opt_mod),method="slr")
  ggsave(plot=llCnet_nOpt[["net"]],sprintf("%s/SLR_neganet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/SLR_kendtest_neganet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_nOpt[["cor"]])
  sink()
  
  
  ################3
  cur_dir2 <- sprintf("%s/%s",cur_dir,"Compare_wo_latent")
  dir.create(cur_dir2)
  
  opt_deg <- deg(rtablist[[i]],r.slr=slrlist[[i]],taxa_f,method="slr")
  deg_mb <- deg(rtablist[[i]],r.slr=mblist[[i]],taxa_f,method="mb")
  deg_gl <- deg(rtablist[[i]],r.slr=glassolist[[i]],taxa_f,method="glasso")
  
  dat_p <- data.frame(OTU=opt_deg$OTU,
                      SLR_Degree=opt_deg$posi,
                      MB_Degree=deg_mb[opt_deg$OTU,"posi"],
                      glasso_Degree=deg_gl[opt_deg$OTU,"posi"],
                      taxa_f[opt_deg$OTU,])
  
  dat_n <- data.frame(OTU=opt_deg$OTU,
                      SLR_Degree=opt_deg$nega,
                      MB_Degree=deg_mb[opt_deg$OTU,"nega"],
                      glasso_Degree=deg_gl[opt_deg$OTU,"nega"],
                      taxa_f[opt_deg$OTU,])
  
  write.csv(dat_p[unique(c(dat_p[order(dat_p$SLR_Degree,decreasing = T)[1:20],"OTU"],
                           dat_p[order(dat_p$MB_Degree,decreasing = T)[1:20],"OTU"],
                           dat_p[order(dat_p$glasso_Degree,decreasing = T)[1:20],"OTU"])),],
            sprintf("%s/positive_degreee_top20_%s.csv",cur_dir2,ident_th))
  
  
  write.csv(dat_n[unique(c(dat_n[order(dat_n$SLR_Degree,decreasing = T)[1:20],"OTU"],
                           dat_n[order(dat_n$MB_Degree,decreasing = T)[1:20],"OTU"],
                           dat_n[order(dat_n$glasso_Degree,decreasing = T)[1:20],"OTU"])),],
            sprintf("%s/negative_degreee_top20_%s.csv",cur_dir2,ident_th))
  
  #############
  slr_mb_p <- Comp_deg(net2=slrlist[[i]],
                       net1=mblist[[i]],
                       method2="slr",
                       method1="mb",
                       comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],NP="posi")
  slr_gl_p <- Comp_deg(net2=slrlist[[i]],
                       net1=glassolist[[i]],
                       method2="slr",
                       method1="glasso",
                       comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                       NP="posi")
  
  plot_grid(slr_mb_p[["net"]],
            slr_gl_p[["net"]],
            align = "hv",
            nrow=1,
            labels = c("a","b"),label_size = 24)
  
  ggsave(sprintf("%s/SLR_vs_MB-glasso_deg_posi_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  
  pcomp2 <- slr_mb_p[["net"]]+facet_wrap(~Guild2,nrow=2)
  
  ggsave(plot=pcomp2,
         sprintf("%s/Opt_vs_MB_deg_guilds_posi_%s.pdf",cur_dir2,ident_th),h=9,w=18)
  
  
  pcomp3 <- slr_gl_p[["net"]]+facet_wrap(~Guild2,nrow=2)
  
  ggsave(plot=pcomp3,
         sprintf("%s/Opt_vs_glasso_deg_guilds_posi_%s.pdf",cur_dir2,ident_th),h=9,w=18)
  ###
  
}


