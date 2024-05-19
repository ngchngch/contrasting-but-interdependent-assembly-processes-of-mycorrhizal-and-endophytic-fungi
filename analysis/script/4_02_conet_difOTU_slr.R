
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
library(cowplot)
lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "4_co-occurrence_network_slr"

dir.create(current_dir)
save.dir <- sprintf("%s/02_conet_analysis_slr",current_dir)

dir.create(save.dir)

source("function/conet_functions_240127.R")


colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####
##jSDM result path
jsdm_data_path <- "XXXXXXXXX"


taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

seq <- read.fasta("data_process/05_Merge_sequence_output/Fungi/0.97/0.97/OTUseq_0.97.fasta")

######
all_lat <- list.files(path="4_co-occurrence_network_slr/01_conet_estimate",
                      pattern = "conet_all_slr",recursive = T,
                      full.names = T)

netlist <- lapply(all_lat,readRDS)

names(netlist) <- sapply(strsplit(all_lat,"/"),function(x)x[3])

#####
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

names(ll_b2) <- sapply(strsplit(lls,"/"),function(x)x[3])

##########

for(i in 1:length(netlist)){
  show.progress(i,1:length(netlist))
  #i <-2
  ident_th <- names(netlist)[i]
  ll <- ll_b2[[i]]
  
  cur_dir <- sprintf("%s/%s",save.dir,ident_th)
  dir.create(cur_dir)
  #################
  #BIC with latent models
  bic <-R_BIC(netlist[[i]],netlist[[i]]$rank0$est$data)
  
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
  slr_net <- netlist[[i]]
  opt_slr <- slr_net[[names(opt_mod)]]
  ew <- Edge_weight(net=opt_slr,zero=T,method="slr")
  
  
  ew1 <- ew[which(taxa_f[ew$OTU1,"Guild"] != "AMF" | taxa_f[ew$OTU2,"Guild"] != "AMF"),]
  ew_p <- ew1[order(ew1$cor),]
  
  ew_n <- ew1[order(ew1$cor,decreasing=T),]
  
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
  
  write(c("positive_pair_number",sum(ew>0),
          "negative_pair_number",sum(ew<0)),
        sprintf("%s/correlation_pair_numbers.txt",cur_dir)
  )
  
  ######################################
  
  hic_otu <- unique(c(ew_p$OTU1,ew_p$OTU2,
                      ew_n$OTU1,ew_n$OTU2))
  
  write.fasta(seq[hic_otu],names=hic_otu,
              sprintf("%s/high_correlattion_pairs_seq_%s.fasta",save.dir,ident_th))
  
  ################
  #drow co-net
  #slr
  nlay <- opt_slr
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
  slr <- slr.cor
  
  netw <- igraph::graph.adjacency(slr, mode='undirected', diag=FALSE, weighted=TRUE)
  
  
  am.coord <-layout_with_graphopt(netw,
                                  charge=igraph::degree(netw)/(50*max(igraph::degree(netw))),
                                  mass=30*((igraph::degree(netw)-min(igraph::degree(netw)))/(max(igraph::degree(netw)-min(igraph::degree(netw))))))
  
  rownames(am.coord) <- colnames(nlay$est$data)
  
  # am.coord <- readRDS(sprintf("%s/best_coord.rds",cur_dir))
  slr_p <- drowConet(net=opt_slr,comm=opt_slr$est$data,NP="posi",col_list = colfng,method="slr",module=F,
                     ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  

  slr_p
  
  slr_m <- drowConet(net=opt_slr,comm=opt_slr$est$data,NP="posi",col_list = colfng,method="slr",module=T,
                     ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  
  slr_m
  
  slr_n <- drowConet(net=opt_slr,comm=opt_slr$est$data,NP="nega",col_list = colfng,method="slr",module=F,
                     ord=am.coord,taxa_f=taxa_f,df_root2=df_sl,save.dir=cur_dir,name=ident_th)
  
  cowplot::plot_grid(slr_p$net+theme(legend.position = "none"),
                     slr_n$net+theme(legend.position = "none"),
            align = "hv",nrow=2,
            labels = c("a","b"),
            label_size = 24)
  
  ggsave(sprintf("%s/merge_conet.pdf",cur_dir),w=7,h=13)
  
  ggsave(plot=slr_m$net,sprintf("%s/posi_conet_module.pdf",cur_dir),h=7,w=7)
  
  saveRDS(am.coord,sprintf("%s/best_coord.rds",cur_dir))
  
  #nodule_membermodule=T,
  
  mod_info <- slr_p$data
  write.csv(unique(mod_info[,c("degree","mod","OTU","taxa")]),
            sprintf("%s/module_OTU_list.csv",cur_dir))
  ###########3#####
  opt_deg <- deg(opt_slr$est$data,r.slr=opt_slr,taxa_f,method="slr")
  dat_p <- data.frame(OTU=opt_deg$OTU,
                      SLR_Degree=opt_deg$posi,
                      taxa_f[opt_deg$OTU,])
  
  dat_n <- data.frame(OTU=opt_deg$OTU,
                      SLR_Degree=opt_deg$nega,
                      taxa_f[opt_deg$OTU,])
  
  write.csv(dat_p[dat_p[order(dat_p$SLR_Degree,decreasing = T)[1:20],"OTU"],],
            sprintf("%s/positive_degreee_top20_%s.csv",cur_dir,ident_th))
  
  
  write.csv(dat_n[dat_n[order(dat_n$SLR_Degree,decreasing = T)[1:20],"OTU"],],
            sprintf("%s/negative_degreee_top20_%s.csv",cur_dir,ident_th))
  
  #########
  #compare llr conet degree
  cur_dir2 <- sprintf("%s/%s",cur_dir,"Compare_llr_conet")
  dir.create(cur_dir2)
  
  llCnet_pOpt <- Comp_llCnet(slr=opt_slr,comm=opt_slr$est$data,taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="posi",R=names(opt_mod),method="slr")
  ggsave(plot=llCnet_pOpt[["net"]],sprintf("%s/SLRposinet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  
  sink(sprintf("%s/SLR_kendtest_posinet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  
  print(llCnet_pOpt[["cor"]])
  sink()
  #####
  
  ###########
  
  llCnet_nOpt <- Comp_llCnet(slr=opt_slr,comm=opt_slr$est$data,taxa_f=taxa_f,ll=ll_b2[[i]],
                             NP="nega",R=names(opt_mod),method="slr")
  ggsave(plot=llCnet_nOpt[["net"]],sprintf("%s/SLR_neganet_deg_llr_Optlatent_%s.pdf",cur_dir2,ident_th),
         h=8,w=9)
  sink(sprintf("%s/SLR_kendtest_neganet_deg_llr_Optlatent_%s.txt",cur_dir2,ident_th))
  print(llCnet_nOpt[["cor"]])
  sink()
  
  
  ################3
  g<- plot_grid(llCnet_pOpt$net+
                  labs(x="Number of positive associations",fill="Guild",
                       y='Explanatory power of OTU covariance in jSDM\n[ log-likelihood ratio ("Full" vs "P + S + Sp") ]')+
              geom_blank(aes(y=max))+
                theme_light()+
              theme(legend.position = "none",
                    axis.text = element_text(size=15),
                    axis.title.y = element_text(size=15),
                    axis.title.x = element_text(size=15),
                    panel.background = element_rect(colour = "black",linewidth=1.1)),
            llCnet_nOpt$net+
              labs(x="Number of negative associations",fill="Guild",
                   y='Explanatory power of OTU covariance in jSDM\n[ log-likelihood ratio ("Full" vs "P + S + Sp") ]')+
              geom_blank(aes(y=max))+
              theme_light()+
              theme(legend.position = "none",
                    axis.text = element_text(size=15),
                    axis.title.y = element_text(size=15),
                    axis.title.x = element_text(size=15),
                    panel.background = element_rect(colour = "black",linewidth=1.1)),
            align = "v",hjust=0,vjust=1,
            labels = c("a","b"),
            label_size = 24)
  g
  ggsave(sprintf("%s/merge_llr_conet_degree_%s.pdf",cur_dir2,ident_th),h=7,w=14.5)
}


