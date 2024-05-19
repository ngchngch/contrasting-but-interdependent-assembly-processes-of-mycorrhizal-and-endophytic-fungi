
##read packages
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggtext';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'igraph';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggnetwork';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'graphlayouts';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- "NSM3";library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "2_jSDM"

dir.create(current_dir)

save.dir <- sprintf("%s/03_jsdm_model_comp_functinal_group",current_dir)

dir.create(save.dir)

######
taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]
######

ll <- readRDS(sprintf("%s/jsdm_result_in_google_colab/OTU_0.97/restricted_lls_mat_sPC90Per.rds",current_dir))

a <- readRDS(sprintf("%s/jsdm_result_in_google_colab/OTU_0.97/jsdm_sdist_result_iter100_Msamp5000_sPC90Per.rds",current_dir))

seq <- read.fasta("data_process/05_Merge_sequence_output/Fungi/0.97/OTUseq_0.97.fasta")

#####
#####
ll_otu <- sapply(ll,colSums)

ll_df <- data.frame(OTU=a$species,ll_otu)

####

ll_comp <- data.frame(OTU=ll_df$OTU,taxa_f[ll_df$OTU,],
                      Plant=ll_df$SoSp-ll_df$NoBio,
                      Soil=ll_df$PlSp-ll_df$NoBio,
                      OTU_Cov=ll_df$NoBio-ll_df$Full,
                      Spatial=ll_df$PlSo-ll_df$NoBio)

saveRDS(ll_comp,sprintf("%s/loglike_ratio_factors.rds",save.dir))

ll_comp$Guild2 <- factor(ll_comp$Guild,
                         levels=c(setdiff(unique(ll_comp$Guild),
                                          c("Unassigned","Unidentified")),
                                  "Unassigned","Unidentified"))
top <- 20
factor<- "Plant"
llr_pl <-ll_comp[order(ll_comp[,factor],decreasing=T)[1:top],c("OTU",factor)]
write.csv(cbind(llr_pl,taxa_f[llr_pl$OTU,]),
          sprintf("%s/llr_Top%s_%s.csv",save.dir,top,factor))

write.fasta(seq[llr_pl$OTU],names=llr_pl$OTU,
            sprintf("%s/llr_Top%s_%s.fasta",save.dir,top,factor))

factor<- "Soil"
llr_so <-ll_comp[order(ll_comp[,factor],decreasing=T)[1:top],c("OTU",factor)]
write.csv(cbind(llr_so,taxa_f[llr_so$OTU,]),
          sprintf("%s/llr_Top%s_%s.csv",save.dir,top,factor))

write.fasta(seq[llr_so$OTU],names=llr_so$OTU,
            sprintf("%s/llr_Top%s_%s.fasta",save.dir,top,factor))

factor<- "Spatial"
llr_sp <-ll_comp[order(ll_comp[,factor],decreasing=T)[1:top],c("OTU",factor)]
write.csv(cbind(llr_sp,taxa_f[llr_sp$OTU,]),
          sprintf("%s/llr_Top%s_%s.csv",save.dir,top,factor))

write.fasta(seq[llr_sp$OTU],names=llr_sp$OTU,
            sprintf("%s/llr_Top%s_%s.fasta",save.dir,top,factor))

factor<- "OTU_Cov"
llr_cv <-ll_comp[order(ll_comp[,factor],decreasing=T)[1:top],c("OTU",factor)]
write.csv(cbind(llr_cv,taxa_f[llr_cv$OTU,]),
          sprintf("%s/llr_Top%s_%s.csv",save.dir,top,factor))

write.fasta(seq[llr_cv$OTU],names=llr_cv$OTU,
            sprintf("%s/llr_Top%s_%s.fasta",save.dir,top,factor))


#####################
############################
mj_gld <-table(taxa_f[ll_df$OTU,"Guild"])[table(taxa_f[ll_df$OTU,"Guild"])>5]
mj_gld2 <- setdiff(names(mj_gld),c("Unassigned","Unidentified"))
ll_df2 <- ll_df[which(taxa_f[ll_df$OTU,"Guild"] %in% mj_gld2),]

ll_comp2 <- data.frame(OTU=ll_df2$OTU,taxa_f[ll_df2$OTU,],
                       Host_Plant=ll_df2$SoSp-ll_df2$NoBio,
                       Soil_Fungal_Community=ll_df2$PlSp-ll_df2$NoBio,
                       OTU_Covariance=ll_df2$NoBio-ll_df2$Full,
                       Spatial_Autocorrelation=ll_df2$PlSo-ll_df2$NoBio)

gll <- gather(ll_comp2,Factor,llr,-c(1:(ncol(taxa_f)+1)))


gll$Guild2 <- factor(gll$Guild,
                     levels=c("EcMF","AMF","Endophyte","Other_RAF"))

gll$Factor2 <- gsub("_"," ",gll$Factor)

gll$Factor2 <- factor(gll$Factor2,levels=c("Host Plant",
                                           "Soil Fungal Community",
                                           "Spatial Autocorrelation",
                                           "OTU Covariance"
))


###########

# 正確法
tab <- data.frame(c("EcMF","AMF","Endophyte","Other_RAF"),
                  c(1,2,3,4))
rownames(tab) <- tab[,1]
ll_comp2$Guild2 <- tab[ll_comp2$Guild,2]

set.seed(0)
p_exact = pSDCFlig(ll_comp2$Host_Plant,
                   ll_comp2$Guild2)

p_exact_s = pSDCFlig(ll_comp2$Soil_Fungal_Community,
                     ll_comp2$Guild2)

p_exact_sp = pSDCFlig(ll_comp2$Spatial_Autocorrelation,
                      ll_comp2$Guild2)

p_exact_cov = pSDCFlig(ll_comp2$OTU_Covariance,
                       ll_comp2$Guild2)

res <- data.frame(label=p_exact$labels,Plant=p_exact$p.val,
                  Soil=p_exact_s$p.val,
                  Spatial=p_exact_sp$p.val,
                  Cov=p_exact_cov$p.val)

res$label <- gsub("4","Other_RAF",gsub("3","Endophyte",gsub("2","AMF",gsub("1","EcMF",res$label))))

write.csv(res,sprintf("%s/Steel-Dwass_test.csv",save.dir))

#res <- read.csv(sprintf("%s/Steel-Dwass_test.csv",save.dir),row.names=1)
annos <- list(Plant=c("b","a","ab","a"),
              Soil=c("a","b","ab","b"),
              Space=c("a","ab","bc","c"),
              Cov=c("a","c","a","ab"))

###########
limits <- aggregate(llr ~ Factor2, gll, function(x) c(ymin = quantile(x, 0.05), ymax = max(x)*1.1))
upper_limits <- limits$llr[,2]
names(upper_limits) <- limits$Factor2
gll2 <- cbind(gll,upper_limit=upper_limits[gll$Factor2])


gpl <- ggplot(gll2[gll2$Factor2=="Host Plant",],aes(x=Guild2,y=as.numeric(llr)))+
  geom_boxplot(aes(fill=Guild2),outlier.colour = NA,show.legend = F)+
  geom_blank(aes(y = upper_limit), inherit.aes = FALSE)+
  geom_star(position = position_jitter(width=0.2),alpha=0.7,fill="white",
            starshape=15)+
  scale_fill_manual(values=colfng_gld[levels(gll$Guild2),2])+
  stat_summary(geom = 'text', label =annos$Plant, fun = max, vjust = -1,size=8)+
  labs(x="Guild",title = "Host plant",
       y='Explanatory power of host plant\n[ log-likelihood ratio ("P + S + Sp" vs "S + Sp") ]')+
  
  theme_light()+
  theme(plot.title = element_text(size=28,face="bold",hjust=0.5),
        legend.text = element_text(size=15),
        axis.title.x = element_markdown(size=25),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=45,hjust = 1,vjust = 1),
        panel.background = element_rect(colour = "black",linewidth=1.1))


gpl

factor <- "soil fungal community"
gso <- ggplot(gll2[gll2$Factor2=="Soil Fungal Community",],aes(x=Guild2,y=as.numeric(llr)))+
  geom_boxplot(aes(fill=Guild2),outlier.colour = NA,show.legend = F)+
  geom_blank(aes(y = upper_limit), inherit.aes = FALSE)+
  geom_star(position = position_jitter(width=0.2),alpha=0.7,fill="white",
            starshape=15)+
  scale_fill_manual(values=colfng_gld[levels(gll$Guild2),2])+
  stat_summary(geom = 'text', label =annos$Soil, fun = max, vjust = -1,size=8)+
  labs(x="Guild",title = "Soil fungal community",
       y=sprintf('Explanatory power of %s\n[ log-likelihood ratio ("P + S + Sp" vs "P + Sp") ]',
                 factor))+
  
  theme_light()+
  theme(plot.title = element_text(size=28,face="bold",hjust=0.5),
        legend.text = element_text(size=15),
        axis.title.x = element_markdown(size=25),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=45,hjust = 1,vjust = 1),
        panel.background = element_rect(colour = "black",linewidth=1.1))


gso

factor <- "spatial autocorrelation"
gsp <- ggplot(gll2[gll2$Factor2=="Spatial Autocorrelation",],aes(x=Guild2,y=as.numeric(llr)))+
  geom_boxplot(aes(fill=Guild2),outlier.colour = NA,show.legend = F)+
  geom_blank(aes(y = upper_limit), inherit.aes = FALSE)+
  geom_star(position = position_jitter(width=0.2),alpha=0.7,fill="white",
            starshape=15)+
  scale_fill_manual(values=colfng_gld[levels(gll$Guild2),2])+
  stat_summary(geom = 'text', label =annos$Space, fun = max, vjust = -1,size=8)+
  labs(x="Guild",title = "Spatial autocorrelation",
       y=sprintf('Explanatory power of %s\n[ log-likelihood ratio ("P + S + Sp" vs "P + S") ]',
                 factor))+
  
  theme_light()+
  theme(plot.title = element_text(size=28,face="bold",hjust=0.5),
        legend.text = element_text(size=15),
        axis.title.x = element_markdown(size=25),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=45,hjust = 1,vjust = 1),
        panel.background = element_rect(colour = "black",linewidth=1.1))


gsp


factor <- "OTU covariance"
gcov <- ggplot(gll2[gll2$Factor2== "OTU Covariance",],aes(x=Guild2,y=as.numeric(llr)))+
  geom_boxplot(aes(fill=Guild2),outlier.colour = NA,show.legend = F)+
  geom_blank(aes(y = upper_limit), inherit.aes = FALSE)+
  geom_star(position = position_jitter(width=0.2),alpha=0.7,fill="white",
            starshape=15)+
  scale_fill_manual(values=colfng_gld[levels(gll$Guild2),2])+
  stat_summary(geom = 'text', label =annos$Cov, fun = max, vjust = -1,size=8)+
  labs(x="Guild",title = factor,
       y=sprintf('Explanatory power of %s\n[ log-likelihood ratio ("Full" vs "P + S + Sp") ]',
                 factor))+
  
  theme_light()+
  theme(plot.title = element_text(size=28,face="bold",hjust=0.5),
        legend.text = element_text(size=15),
        axis.title.x = element_markdown(size=25),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18,angle=45,hjust = 1,vjust = 1),
        panel.background = element_rect(colour = "black",linewidth=1.1))


gcov

cowplot::plot_grid(gpl+
            theme(axis.title.x=element_blank(),
                  plot.title = element_text(size=20)),
          gso+
            theme(axis.title.x=element_blank(),
                  plot.title = element_text(size=20)),
          gsp+
            theme(axis.title.x=element_blank(),
                  plot.title = element_text(size=20)),
          gcov+
            theme(axis.title.x=element_blank(),
                  plot.title = element_text(size=20)),
          align = "hv",nrow=1,
          labels=c("a","b","c","d"),label_size = 25)

ggsave(sprintf("%s/llr_boxplot.pdf",save.dir),h=7,w=20)

