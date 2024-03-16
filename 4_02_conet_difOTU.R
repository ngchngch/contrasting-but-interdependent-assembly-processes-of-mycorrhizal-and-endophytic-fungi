
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

save.dir <- sprintf("%s/02_conet_analysis",current_dir)

dir.create(save.dir)

source("/Volumes/Transcend/data/Sugadaira_sequence/Final_merge/analysis_series/function/conet_functions_230629.R")


colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####

taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("bioinfo_output/sample_process/fungi/OTUseq_0.97.fasta")

######
all_lat <- list.files(path=conet_data_path,
                      pattern = "conet_slr_alllatent",recursive = T,
                      full.names = T)

netlist <- lapply(all_lat,readRDS)

names(netlist) <- sapply(strsplit(all_lat,"/"),function(x)x[11])

occtab <- list.files(path=conet_data_path,
                     pattern = "jsdm_occtable",recursive = T,
                     full.names = T)


tablist <- lapply(occtab,readRDS)

namlist <- lapply(tablist,colnames)

#raw_read table
rtab <- list.files(path=conet_data_path,
                   pattern = "conet_input",recursive = T,
                   full.names = T)

rtablist <- lapply(rtab,readRDS)

names(rtablist) <- sapply(strsplit(rtab,"/"),function(x)x[4])

##jsdm res
lls <- list.files(path=jsdm_data_path,
                   pattern = "restricted_lls_mat",recursive = T,
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
 

head(ll_list[[1]])
     
ll_b2 <- lapply(ll_list,function(x){
  a <- cbind(llr=x$NoBio-x$Full,taxa_f[x$OTU,])
  rownames(a) <- x$OTU
  
  return(a)
})

names(ll_b2) <- sapply(strsplit(lls,"/"),function(x)x[11])

##########

for(i in 1:length(rtablist)){
  show.progress(i,1:length(rtablist))
  #i <-2
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
  ew <- Edge_weight(slr=netlist[[i]],R=names(opt_mod),zero=T)
  
  
  ew1 <- ew[which(taxa_f[ew$OTU1,"Guild"] != "AMF" | taxa_f[ew$OTU2,"Guild"] != "AMF"),]
  ew_p <- ew1[order(ew1$cor)[1:20],]
  
  ew_n <- ew1[order(ew1$cor,decreasing=T)[1:20],]
  
  write.csv(cbind(ew_p,
                  OTU1=taxa_f[ew_p$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_p$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/Top20_correlation_pairs_P_%s_%s.csv",
                    cur_dir,ident_th,names(opt_mod)))
  write.csv(cbind(ew_n,
                  OTU1=taxa_f[ew_n$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ew_n$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/Top20_correlation_pairs_N_%s_%s.csv",
                    cur_dir,ident_th,names(opt_mod)))
  
  ew_r0 <- Edge_weight(slr=netlist[[i]],R="rank0",zero=T)
  ewr0_1 <- ew_r0[which(taxa_f[ew_r0$OTU1,"Guild"] != "AMF" | taxa_f[ew_r0$OTU2,"Guild"] != "AMF"),]
  
  ewr0_p <- ewr0_1[order(ewr0_1$cor)[1:20],]
  ewr0_n <- ewr0_1[order(ewr0_1$cor,decreasing=T)[1:20],]
  
  write.csv(cbind(ewr0_p,
                  OTU1=taxa_f[ewr0_p$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ewr0_p$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/Top20_correlation_pairs_P_%s_%s.csv",
                    cur_dir,ident_th,"rank0"))
  write.csv(cbind(ew_n,
                  OTU1=taxa_f[ewr0_n$OTU1,c("Order","Genus","Guild")],
                  OTU2=taxa_f[ewr0_n$OTU2,c("Order","Genus","Guild")]),
            sprintf("%s/Top20_correlation_pairs_N_%s_%s.csv",
                    cur_dir,ident_th,"rank0"))
  
  hic_otu <- unique(c(ew_p$OTU1,ew_p$OTU2,
           ew_n$OTU1,ew_n$OTU2,
           ewr0_p$OTU1,ewr0_n$OTU2,
           ewr0_n$OTU1,ewr0_n$OTU2))
  
  write.fasta(seq[hic_otu],names=hic_otu,
              sprintf("%s/high_correlattion_pairs_seq_%s.fasta",save.dir,ident_th))
  
  ################
  #drow co-net
  nlay <- netlist[[i]]
  r.slr <- nlay[["rank0"]]
  #optimal inverse covariance -> covariance
  slr.icov <- solve(as.matrix(r.slr$est$icov[[getOptInd(r.slr)]])) 
  #covariance -> correlation
  slr.cor <- cov2cor(slr.icov)
  
  #optimal network
  slr.refit <- as.matrix(getRefit(r.slr))
  sum(slr.refit)
  #extract edge weight of optimal network
  slr.cor[which(slr.refit==FALSE)] <- 0
  diag(slr.cor) <- 0
  slr <- slr.cor
  slr[which(slr <= 0)] <- 0 
  
  netw <- igraph::graph.adjacency(slr, mode='undirected', diag=FALSE, weighted=TRUE)
  
  am.coord <-layout_nicely(netw)
  rownames(am.coord) <- colnames(rtablist[[i]])
    
  drowConet(net=netlist[[i]],R=names(opt_mod),comm=rtablist[[i]],NP="posi",col_list = colfng,
            ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  drowConet(net=netlist[[i]],R="rank0",comm=rtablist[[i]],NP="posi",col_list = colfng,
            ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  drowConet(net=netlist[[i]],R=names(opt_mod),comm=rtablist[[i]],NP="nega",col_list = colfng,
            ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  drowConet(net=netlist[[i]],R="rank0",comm=rtablist[[i]],NP="nega",col_list = colfng,
            ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  
  ###########3#####
  
 
  #########
  #compare llr conet degree
  cur_dir2 <- sprintf("%s/%s",cur_dir,"Compare_llr_conet")
  dir.create(cur_dir2)
  
  llCnet_p0 <- Comp_llCnet(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
              NP="posi",R="rank0")
  ggsave(plot=llCnet_p0[["net"]],sprintf("%s/posinet_deg_llr_Nolatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/kendtest_posinet_deg_llr_Nolatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_p0[["cor"]])
  sink()
  
  llCnet_pOpt <- Comp_llCnet(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                           NP="posi",R=names(opt_mod))
  ggsave(plot=llCnet_pOpt[["net"]],sprintf("%s/posinet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/kendtest_posinet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  
  print(llCnet_pOpt[["cor"]])
  sink()
  
  llCnet_n0 <- Comp_llCnet(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                           NP="nega",R="rank0")
  ggsave(plot=llCnet_n0[["net"]],sprintf("%s/neganet_deg_llr_Nolatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/kendtest_neganet_deg_llr_Nolatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_n0[["cor"]])
  sink()
  
  llCnet_nOpt <- Comp_llCnet(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="nega",R=names(opt_mod))
  ggsave(plot=llCnet_nOpt[["net"]],sprintf("%s/neganet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/kendtest_neganet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_nOpt[["cor"]])
  sink()
  
  
  ################3
  cur_dir2 <- sprintf("%s/%s",cur_dir,"Compare_wo_latent")
  dir.create(cur_dir2)
  
  net_slr <- netlist[[i]]
  r.slr <- net_slr[[names(opt_mod)]]
  
  opt_deg <- deg(rtablist[[i]],r.slr,colnames(r.slr$est$data),taxa_f)
  deg0 <- deg(rtablist[[i]],r.slr=net_slr[["rank0"]],colnames(r.slr$est$data),taxa_f)
  
  dat_p <- data.frame(OTU=opt_deg$OTU,
                     Opt_Degree=opt_deg$posi,
                     Degree0=deg0[opt_deg$OTU,"posi"],
                     taxa_f[opt_deg$OTU,])
  
  dat_n <- data.frame(OTU=opt_deg$OTU,
                      Opt_Degree=opt_deg$nega,
                      Degree0=deg0[opt_deg$OTU,"nega"],
                      taxa_f[opt_deg$OTU,])
  
  write.csv(dat_p[union(dat_p[order(dat_p$Opt_Degree,decreasing = T)[1:20],"OTU"],
                        dat_p[order(dat_p$Degree0,decreasing = T)[1:20],"OTU"]),],
            sprintf("%s/positive_degreee_top20_%s.csv",cur_dir2,ident_th))
  
  
  write.csv(dat_n[union(dat_n[order(dat_n$Opt_Degree,decreasing = T)[1:20],"OTU"],
                        dat_n[order(dat_n$Degree0,decreasing = T)[1:20],"OTU"]),],
            sprintf("%s/negative_degreee_top20_%s.csv",cur_dir2,ident_th))
  
  #############
  
  pnetComp_wolat <- Comp_wo_latent(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                   NP="posi",OptR=names(opt_mod))
  
  ggsave(plot=pnetComp_wolat[["net"]],sprintf("%s/Opt_vs_nolatent_deg_posi_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  
  pcomp2 <- pnetComp_wolat[["net"]]+facet_wrap(~Guild2,nrow=2)
  
  ggsave(plot=pcomp2,
         sprintf("%s/Opt_vs_nolatent_deg_guilds_posi_%s.pdf",cur_dir2,ident_th),h=9,w=18)
  ###
  
  
  nnetComp_wolat <- Comp_wo_latent(slr=netlist[[i]],comm=rtablist[[i]],taxa_f=taxa_f,ll=ll_b2[[i]],
                                     NP="nega",OptR=names(opt_mod))
  
  ggsave(plot=nnetComp_wolat[["net"]],sprintf("%s/Opt_vs_nolatent_deg_nega_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  
  ncomp2 <- nnetComp_wolat[["net"]]+facet_wrap(~Guild2,nrow=2)
  
  ggsave(plot=ncomp2,
         sprintf("%s/Opt_vs_nolatent_deg_guilds_nega_%s.pdf",cur_dir2,ident_th),h=9,w=18)

  
  
  }


