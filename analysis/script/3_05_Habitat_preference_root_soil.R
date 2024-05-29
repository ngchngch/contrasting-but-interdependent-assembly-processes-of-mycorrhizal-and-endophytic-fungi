
lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggtext';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggforce';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'cowplot';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- "parallel";library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'glmmTMB';library(package = lib, character.only=TRUE);packageVersion(lib)

####

##set directory

save.dir <- sprintf("%s/05_habitat_preference_root_soil",current_dir)

dir.create(save.dir)

colfng <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")

colfng_gld <- colfng[["Guild"]]

#####
info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")

a <- readRDS("2_jSDM/jsdm_result_in_google_colab/OTU_0.97/jsdm_sdist_result_iter100_Msamp5000_sPC90per.rds")

seq <- read.fasta("data_process/05_Merge_sequence_output/Fungi/0.97/OTUseq_0.97.fasta")

df_rt <- readRDS("0_data_processing/04_rarefaction/OTU97/covrarefy_sqtb_rootF.rds")[[1]]

taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

jsdm_files <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per.rds")

df_soil <- readRDS("0_data_processing/04_rarefaction/OTU97/covrarefy_sqtb_soilF.rds")[[1]]

##
df_rt2 <- df_rt[rownames(jsdm_files$Occ),colnames(jsdm_files$Occ)]

df_root <- df_rt2[rowSums(df_rt2)>0,colSums(df_rt2)>0]/rowSums(df_rt2[rowSums(df_rt2)>0,colSums(df_rt2)>0])

info4 <- info[rownames(df_root),]
####

ll <- readRDS("2_jSDM/jsdm_result_in_google_colab/OTU_0.97/restricted_lls_mat_sPC90per.rds")

llr <- colSums(ll$PlSp - ll$NoBio)

names(llr) <- a$species
  
df <- as.matrix(df_root)
df[df>0] <- 1

df_so <- df_soil/rowSums(df_soil)
df_so2 <- df_so[,colSums(df_so>0)>15]

uniOTU <- intersect(colnames(df_root),colnames(df_so2))

df2 <- df[,uniOTU]

site <- unique(info4$site2)

for(j in 1:length(site)){
  show.progress(j,1:length(site))
snam <- info[which(info$site2 == site[j]),"ID"]

rsnam <- snam[grep("NF",snam)]

plnam <- table(info4[which(rownames(info4) %in% rsnam),"plant"])

for(i in 1:length(plnam)){
  if(plnam[i]==1){#i <-2
    id <- info4[which(rownames(info4) %in% rsnam & info4$plant==names(plnam[i])),"ID"]
    df_sl <- df2[which(rownames(df2) %in% id),]
  }else{
    id <- info4[which(rownames(info4) %in% rsnam & info4$plant == names(plnam[i])),"ID"]
    df_sl <- colSums(df2[which(rownames(df2) %in% id),])
  }

ssnam <- snam[grep("NS",snam)]

df_so_sl <- df_so2[which(rownames(df_so2) %in% ssnam),]


dat <- data.frame(Site=site[j],
                  OTU=uniOTU,
                  plant=names(plnam[i]),
                  Root=df_sl[uniOTU],
                  Root_s_num=sum(rownames(df2) %in% id),
                  Soil=df_so_sl[uniOTU])

if(j==1){
  date <- dat
}else{
  date <- rbind(date,dat)
}
}}


out <- matrix(NA,ncol=4,nrow=length(unique(date$OTU)))
nonconv <- c()
for(k in 1:length(unique(date$OTU))){
show.progress(k , 1:length(unique(date$OTU)))
  otu <- unique(date$OTU)[k]
  
  date2 <- date[which(date$OTU == otu ),]
  
  gres <- glm(Root~scale(log(Soil+1))+plant,offset=log(Root_s_num),
                   data=date2,family = poisson(link = "log"))
  
  gres2 <- summary(gres)
  
  
  out[k,] <- gres2$coefficients[2,]
  
}

dimnames(out) <- list(unique(date$OTU),
                      c("cor","se","z","raw_p"))
  
  
outdat <- as.data.frame(out[which(!rownames(out) %in% nonconv),])

outdat$p_BH <- p.adjust(outdat$raw_p,method = "BH")

pick_otu <- rownames(outdat[outdat$p_BH<0.05,])

uniOTU2 <- setdiff(uniOTU,nonconv)

###
rs <- data.frame(root=colMeans(df_root[,uniOTU2]),
                 root_Occ=colSums(df_root[,uniOTU2]>0),
                 soil=colMeans(df_so2[,uniOTU2]),
                 soil_Occ=colSums(df_so2[,uniOTU2]>0),
                 coef=outdat[uniOTU2,"cor"],
                 llr=llr[uniOTU2],
                 p=outdat[uniOTU2,"p_BH"],
                 sig=sapply(outdat[uniOTU2,"p_BH"],sig,star=F),
                 taxa_f[uniOTU2,])

write.csv(rs,sprintf("%s/habitat_preference_root_soil.csv",save.dir))

##
write.csv(rs[order(rs$root_Occ,decreasing = T)[1:20],-c(5,6,7,8)],
          sprintf("%s/top_20_frequent_in_root_OTU.csv",save.dir))
write.csv(rs[order(rs$soil_Occ,decreasing = T)[1:20],-c(5,6,7,8)],
          sprintf("%s/top_20_frequent_in_soil_OTU.csv",save.dir))
##
rs$Guild2 <- factor(rs$Guild,levels=c(setdiff(unique(rs$Guild),
                                                c("Unassigned","Unidentified")),
                                        "Unassigned","Unidentified"))

rs$signif2 <- gsub("p","*P* ",rs$sig)


rs$signif2 <- gsub("<","< ",rs$signif2)

rs[which(rs$signif2 == ""),"signif2"] <- "N.S." 

write.csv(rs,sprintf("%s/root_soil_ab_Occ.csv",save.dir))

gab <- ggplot(rs,aes(x=log10(soil+1),y=log10(root+1)))+
  geom_star(starshape=15,aes(fill=Guild),size=4)+
  labs(x="log10 (relative abundance in soil +1)",
       y="log10 (relative abundance in root +1)",
       fill="Guild")+
  scale_fill_manual(values=colfng_gld[levels(rs$Guild2),2])+
  #scale_starshape_manual(values=c(1,23,28,15))+
  theme_light()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        strip.text =element_markdown(size=15),
        legend.text =element_markdown(size=15),
        legend.title =element_markdown(size=18),
        panel.background = element_rect(colour = "black",linewidth=1.1))+
  guides(fill=guide_legend(override.aes = list(starshape=15)))+
  coord_cartesian(xlim=c(0,log10(max(c(rs$root,rs$soil))+1)),
                  ylim=c(0,log10(max(c(rs$root,rs$soil))+1)))

gab

ggsave(sprintf("%s/root_soil_abund.png",save.dir),
       h=7,w=9,dpi=300)


goc <- ggplot(rs,aes(x=soil_Occ,y=root_Occ))+
  geom_star(starshape=15,aes(fill=Guild),size=4)+
  labs(x="Occurence in soil",
       y="Occurence in root",
       fill="Guild")+
  scale_fill_manual(values=colfng_gld[levels(rs$Guild2),2])+
  #scale_starshape_manual(values=c(1,23,28,15))+
  theme_light()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        strip.text =element_markdown(size=15),
        legend.text =element_markdown(size=15),
        legend.title =element_markdown(size=18),
        panel.background = element_rect(colour = "black",linewidth=1.1))+
  guides(fill=guide_legend(override.aes = list(starshape=15)))

goc

ggsave(sprintf("%s/root_soil_occurence.png",save.dir),
       h=7,w=9,dpi=300)

merge_plot <- plot_grid(goc+theme(legend.position = "none"),
                        gab+theme(legend.position = "none"),
          labels = c("a","b"),label_size = 24,vjust = 1)

ggsave(plot=merge_plot,sprintf("%s/root_soil_Occurence_Abundance.pdf",save.dir),
       h=6.5,w=13)

ggsave(plot=g_legend(gab),sprintf("%s/legend_root_soil_Occurence_Abundance.pdf",save.dir),
       h=4,w=4)
write.csv(rs[pick_otu,],sprintf("%s/signifOTU_rs_GLM.csv",save.dir))
write.fasta(seq[pick_otu],names=pick_otu,sprintf("%s/signifOTU_rs_GLM.fasta",save.dir))

##########################
#boxplot

mj_gld <- setdiff(names(table(rs$Guild2))[which(table(rs$Guild2)>5)],c("Unassigned","Unidentified"))

rs_blank <- rs 

rs_blank$Guild2 <-"EcMF"


##
# 正確法
tab <- data.frame(mj_gld,
                  c(1:length(mj_gld)))
rownames(tab) <- tab[,1]

set.seed(0)
p_exact = NSM3::pSDCFlig(rs[which(rs$Guild %in% mj_gld),"coef"],
                         rs[which(rs$Guild %in% mj_gld),"Guild2"])


res <- data.frame(label=p_exact$labels,
                  coef=p_exact$p.val)

write.csv(res,sprintf("%s/Steel_Dwass_res.csv",cur_dir))

#res <- write.csv(sprintf("%s/Steel_Dwass_res.csv",save.dir),row.names=1)
#######
annos <- c("b","a","a")

gbox <- ggplot(rs[which(rs$Guild2 %in% mj_gld),],aes(x=Guild2,y=coef))+
  geom_hline(yintercept = 0,color="blue")+
  geom_blank(data=rs_blank,aes(x=Guild2,y=coef+0.1,fill=Guild2),show.legend = F)+
  geom_boxplot(aes(fill=Guild2),outlier.colour = NA,show.legend = F)+
  geom_star(aes(starshape=signif2,size=soil),fill="white",alpha=0.5,
            position = position_jitter(width=0.2),show.legend = F)+
  labs(x="",
       y="Partial regression coefficient in the GLM of root-soil relationship",
       starshape="*P* (FDR)",
       fill="Functional group")+
  stat_summary(geom = 'text', label =annos, fun = max, vjust = -1,size=6)+
  scale_fill_manual(values=colfng_gld[levels(rs$Guild2),2])+
  scale_starshape_manual(values=c(1,23,28,15))+
  #facet_wrap(~Guild2,nrow=2)+
  theme_light()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        strip.text =element_markdown(size=15),
        legend.text =element_markdown(size=15),
        legend.title =element_markdown(size=18))+
  guides(fill=guide_legend(override.aes = list(starshape=15)))
gbox

ggsave(sprintf("%s/coef_boxplot.png",save.dir),h=7,w=4,dpi=300)

########################
#coef vs llr

gscat <- ggplot(rs,aes(x=llr,y=coef))+
  geom_blank(data=rs_blank,aes(x=llr,y=coef+0.1,fill=Guild2),show.legend = F)+
  geom_hline(yintercept = 0,color="blue")+
  geom_star(aes(starshape=signif2,fill=Guild,size=soil),
            position = position_jitter(width=0.2))+
  labs(x='Explanatory power of soil fungal community\n[ log-likelihood ratio ("P + S + Sp" vs "P + Sp") ]',
       y="Partial regression coefficient in the GLM of root-soil relationship",
       starshape="*P* (FDR) of coefficients",
       fill="Guild",size="Relative abundance in soil")+
  scale_fill_manual(values=colfng_gld[levels(rs$Guild2),2])+
  scale_starshape_manual(values=c(1,23,28,15))+
  #facet_wrap(~Guild2,nrow=2)+
  theme_light()+
  theme(axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=18),
        axis.text = element_text(size=15),
        strip.text =element_markdown(size=15),
        legend.text =element_markdown(size=15),
        legend.title =element_markdown(size=18))+
  guides(fill=guide_legend(override.aes = list(size=3,starshape=15)),
         starshape=guide_legend(override.aes = list(size=3)),
         size=guide_legend(override.aes = list(starshape=15)))

gscat

ggsave(sprintf("%s/soil_coef_in_GLM_vs_llr.pdf",save.dir),h=7,w=10)

########
plot_grid(gbox+
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_text(size=17),
                  axis.text.y = element_text(size=15),
                  axis.text.x = element_text(angle=45,size=14,vjust=1,hjust=1),
                  legend.text =element_blank(),
                  legend.title =element_blank(),
                  panel.background = element_rect(colour = "black",linewidth=1.1)),
          gscat+#labs(y="")+
            theme(axis.title.x = element_text(size=17),
                  axis.title.y = element_text(size=17),
                  axis.text.x = element_text(size=15),
                  axis.text.y = element_text(size=15),
                  legend.text =element_markdown(size=15),
                  panel.background = element_rect(colour = "black",linewidth=1.1),
                  legend.title =element_markdown(size=18)),
          nrow=1,align="hv",labels = c("a","b"),label_size = 19,
          rel_widths = c(0.25,0.9)
          )

ggsave(sprintf("%s/rs_glm_vs_llr_box_scat_plot.pdf",save.dir),h=8.8,w=14)

