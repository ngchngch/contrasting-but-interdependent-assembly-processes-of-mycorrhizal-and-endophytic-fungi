##########
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggtext';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'circlize';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ComplexHeatmap';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)

####

##set directory

save.dir <- sprintf("%s/04_preference_heatmap",current_dir)

dir.create(save.dir)

##############
ll <- readRDS("2_jSDM/jsdm_result_in_google_colab/OTU_0.97/restricted_lls_mat_sPC90per.rds")

jmod <- readRDS("2_jSDM/jsdm_result_in_google_colab/OTU_0.97/jsdm_sdist_result_iter100_Msamp5000_sPC90per.rds")

seq <- read.fasta("bioinfo_output/sample_process/fungi/OTUseq_0.97.fasta")

#read randamized result
z_res<- readRDS(sprintf("%s/03_dprime_2dp_randamize/result_inlinux/result/Z_val_2dp_plant.rds",current_dir))

zmat <- readRDS(sprintf("%s/03_dprime_2dp_randamize/result_inlinux/result/2dp_plant_zmat_plant.rds",current_dir))

dfng <- readRDS(sprintf("%s/03_dprime_2dp_randamize/result_inlinux/result/dprime_zval_fungi.rds",current_dir))
dpl <- readRDS(sprintf("%s/03_dprime_2dp_randamize/result_inlinux/result/dprime_zval_plant.rds",current_dir))

#######
froot <- read.csv("/Volumes/Transcend/data/Sugadaira_sequence/Final_merge/analysis_sequence/metadata/FungalRoot_Table_S2.csv",
                  row.names = 1)

#######

rownames(dpl) <- dpl$name
rownames(dfng) <- dfng$name

#######3
ll_otu <- sapply(ll,colSums)

ll_df <- data.frame(OTU=jmod$species,ll_otu)

ll_pl <- data.frame(OTU=ll_df$OTU,taxa_f[ll_df$OTU,],
                    plant=ll_df$SoSp-ll_df$NoBio)
###########3

mat <- zmat

for(i in 1:nrow(zmat)){
  for(j in 1:ncol(zmat)){
    mat[i,j] <- z_res[which(z_res$OTU == rownames(zmat)[i] & z_res$plant == colnames(zmat)[j]),"p_BH"]
  }
}

mat[is.na(mat)] <- 1

mat2 <- t(mat)

#clustering
taxa_f2 <- taxa_f

#zmat[which(zmat>3 | zmat<3)] <- 0
ord <- c()
for(i in c("AMF","EcMF","Endophyte","Pathogen","Mycoparasite","Other_RAF","Unassigned","Unidentified")){
  nam <- rownames(zmat)[which(rownames(zmat)%in%
                                rownames(rowSelect(taxa_f2,"Guild",i)))]
  zmat2 <- zmat[nam,]
  zmat2[is.na(zmat2)] <- 0
  
  if(length(nam)>1){
    a <- hclust( dist(zmat2), "ward.D2" )
    b <- rownames(zmat2) ; names(b) <- 1:nrow(zmat2)
    ord <- c(ord,b[a$order])
  }else{
    ord <- c(ord,nam)
  }
}

# clustering AM Genus & EcM Genus separetly
pg <- substr(colnames(zmat),2,nchar(colnames(zmat))-1)
froot2 <- froot[pg,]
am_genus <- pg[froot2 =="AM"]
ecm_genus <- pg[grep("EcM",froot2)]

zmat_am <- zmat[,which(colnames(zmat) %in% sprintf("*%s*",am_genus))]

c_am <- hclust( dist(t(zmat_am)), "ward.D2" )

d_am <- rownames(t(zmat_am)) ; names(d_am) <- 1:ncol(zmat_am)

zmat_ecm <- zmat[,which(colnames(zmat) %in% sprintf("*%s*",ecm_genus))]

c_ecm <- hclust( dist(t(zmat_ecm)), "ward.D2" )

d_ecm <- rownames(t(zmat_ecm)) ; names(d_ecm) <- 1:ncol(zmat_ecm)

d <- c(d_ecm[c_ecm$order],d_am[c_am$order])

col1 = colorRamp2(c(min(zmat,na.rm=T), 0, max(zmat,na.rm=T)), c("blue", "white", "red"))

col_dfng <- colorRamp2(c(min(dfng$zval,na.rm=T),0,
                         max(dfng$zval,na.rm=T)),
                       c("purple", "white", "darkorange"))
col_dpl <- colorRamp2(c(min(dpl$zval,na.rm=T),
                        max(dpl$zval,na.rm=T)),
                      c("white", "darkgreen"))

col_list <- readRDS("1_basic_discription/02_barplot/fungi_taxa_color.rds")
colfng_gld <- col_list[["Guild"]]
colfng_gld[grep("Endophyte",colfng_gld[,1]),2] <- "#FFFF00"


ha1 <- HeatmapAnnotation(
                        show_annotation_name = F, 
                        p_value = anno_text(dfng[ord,"signif"],which = "column",
                                                                      just="top",location=0.5,
                                                                      gp=gpar(fontsize=17)),
                        
                        dprime_fungi=dfng[ord,"zval"],
                        Functional_group=taxa_f2[ord,"Guild"],
                        col=list(dprime_fungi=col_dfng,
                                 Functional_group=colfng_gld[taxa_f2[ord,"Guild"],2]),
                       
                        gp = gpar(col = "black"))

rha <- rowAnnotation(p_value = anno_text(dpl[d,"signif"],which = "row",just="top",location=0.5,
                                         gp=gpar(fontsize=17)),
                     dprime_plant=dpl[d,"zval"],
                        show_annotation_name = F,
                        col=list(dprime_plant=col_dpl),
                     
                        gp = gpar(col = "black"))



tha <-  HeatmapAnnotation(log_liklihood_ratio = anno_points(ll_pl[ord,"plant"], pch = 16, size = unit(3, "mm"), 
                                                            axis_param = list(
                                                              gp=gpar(fontsize=11)),
                                     #axis_param = list(at = c(0, 2e5, 4e5, 6e5), 
                                      #                 labels = c("0kb", "200kb", "400kb", "600kb")),
                                     width = unit(4, "cm")))


mat2 <- t(mat[ord,d])

ht1 = Heatmap(t(zmat[ord,d]), 
              col=col1,show_row_names = T,
              show_column_names = F,
              row_names_gp = gpar(fontsize = 14),
              top_annotation = tha,
              bottom_annotation = ha1,right_annotation = rha,
              row_labels = gt_render(d),
              cluster_rows = FALSE, cluster_columns = FALSE,
              na_col="gray90",
              heatmap_legend_param = list(title="Two dimensional preference (2DP)"),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(mat2[i, j] < 0.0005) {
                  grid.text("***", x, y,rot=90,vjust = 0.8,
                            gp=gpar(col="white", fontsize=18))
                } else {if(mat2[i, j] < 0.005) {
                  grid.text("**", x, y,rot=90,vjust = 0.8,
                            gp=gpar(col="white", fontsize=18))
                }else{if(mat2[i, j] < 0.025) {
                  grid.text("*", x, y,rot=90,vjust = 0.8,
                            gp=gpar(col="white", fontsize=18))
                }
                }
                }
                }
              )


ht1

png(sprintf("%s/2dp_Zmap.png",save.dir),height = 1700, width = 8000,res=300)
ht1
dev.off()

write.csv(cbind(z_res[which(z_res$p_BH<0.025),c("OTU","plant","zval","p_BH")],
                taxa_f[z_res[which(z_res$p_BH<0.025),c("OTU")],]),
          sprintf("%s/sdp_signif_OTU.csv",save.dir))

write.csv(cbind(dfng[which(dfng$p_BH<0.05),c("name","zval","p_BH")],
                taxa_f[dfng[which(dfng$p_BH<0.05),c("name")],]),
          sprintf("%s/dprime_signif_OTU.csv",save.dir))

write.fasta(seq[union(z_res[which(z_res$p_BH<0.025),c("name")],dfng[which(dfng$p_BH<0.05),c("name")])],
            names=union(z_res[which(z_res$p_BH<0.025),c("name")],dfng[which(dfng$p_BH<0.05),c("name")]),
            sprintf("%s/sig_OTU_seq.fasta",save.dir))

g <- ggplot(z_res,aes(y=p_BH,x=as.numeric(zval)))+
  
  geom_vline(xintercept = 0,color="blue")+
  geom_hline(yintercept = 0.025,color="red")+
  geom_text(y=0.08,x=-4.2,label="0.025",color="red",size=5)+
  geom_point(color="black",shape=1,size=4)+
  labs(x="Two dimensional preference [ *2DP* ]",y="*P* (FDR)")+
  theme_classic()+
  theme(axis.title.x = element_markdown(size=24),
        axis.title.y = element_markdown(size=24),
        axis.text = element_text(size=15))
g
ggsave(sprintf("%s/2DP_Z_P_distribution.png",save.dir),h=8,w=8,dpi=300)


g <- ggplot(dfng,aes(y=p_BH,x=as.numeric(zval)))+
  geom_hline(yintercept = 0.05,color="red")+
  geom_text(y=0.08,x=-4.2,label="0.05",color="red",size=5)+
  geom_point(color="black",shape=1,size=4)+
  labs(x="*d'* (Fungi)",y="*P* (FDR)")+
  theme_classic()+
  theme(axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=15))
g
ggsave(sprintf("%s/dprime_fungi_Z_P_distribution.png",save.dir),h=8,w=8)


g <- ggplot(dpl,aes(y=p_BH,x=as.numeric(zval)))+
  
  geom_hline(yintercept = 0.05,color="red")+
  geom_text(y=0.08,x=-4.2,label="0.05",color="red",size=5)+
  geom_point(color="black",shape=1,size=4)+
  labs(x="*d'* (Plant)",y="*P* (FDR)")+
  theme_classic()+
  theme(axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text = element_text(size=15))
g
ggsave(sprintf("%s/dprime_plant_Z_P_distribution.png",save.dir),h=8,w=8)
###################################