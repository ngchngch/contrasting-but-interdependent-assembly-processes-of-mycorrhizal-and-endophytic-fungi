#install.packages("pROC")

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
lib <- "pROC";library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

current_dir <- "2_jSDM"

dir.create(current_dir)

save.dir <- sprintf("%s/02_jsdm_model_comp_community_level",current_dir)

dir.create(save.dir)

#function
ConvPath <- function(total_ll,compfact,coord_weight=F){
  net_mat <- matrix(NA,nrow=nrow(total_ll),ncol=nrow(total_ll))
  dimnames(net_mat) <- list(total_ll$model,total_ll$model)
  
  for(i in 1:nrow(total_ll)){
    for(j in 1:nrow(total_ll)){#i <- 2;j <- 1
      lls <- total_ll[c(i,j),]
      lls2 <- lls[order(lls$step,decreasing = T),]
      m1_comp <- compfact[[lls2[2,"model"]]]
      m2_comp <- compfact[[lls2[1,"model"]]]
      
      if(abs(total_ll[i,"step"]-total_ll[j,"step"]) == 1 && all(m1_comp %in% m2_comp)){
        ord_ll <- lls[order(lls$step,decreasing = T),"ll"] 
      }else{
        ord_ll <- NA
      }
      net_mat[j,i] <- ord_ll[2]-ord_ll[1]
    }
  }
  
  ig <- igraph::graph.adjacency(net_mat, mode='undirected', diag=FALSE, weighted=TRUE)
  
  tord <- as.numeric(total_ll[rownames(net_mat),"ll"])
  if(!(coord_weight)){
    coord <- cbind(as.numeric(total_ll[rownames(net_mat),"step"])/max(as.numeric(total_ll[rownames(net_mat),"step"])),
                   (rank(tord))/2)
  }else{
    coord <- cbind(as.numeric(total_ll[rownames(net_mat),"step"]),
                   tord)
  }
  
  gig <- ggnetwork(ig,layout=coord,scale=F)
  
  
  return(gig)
}


drowPath <- function(gig,facter=NULL,step,col=col_mod,log=F,auc=F,scale=T){
  
  gig$Model <- step[gig$name,"true_n"]
  
  if(!is.null(facter)){
    gig$col <- ifelse(str_detect(gig$Model,pattern = facter),gig$Model,"others")
    gig$col2 <- factor(gig$col,levels=c(setdiff(unique(gig$col),"others"),"others"))
  }else{
    gig$col <- gig$Model
    gig$col2 <- factor(gig$col,levels=unique(gig$col))
  }
  
  
  col_mod <- col
  
  gig2 <- gig
  gig2$weight2 <- abs(gig2$weight)
  
  gnet <- ggplot(data = gig2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(data = gig2[which(gig2$weight<0),],na.rm = T,
               color = "blue",aes(linewidth=weight2),alpha=0.5,curvature = 0) +
    geom_edges(data = gig2[which(gig2$weight>0),],na.rm=T,show.legend = F,
               color = "red",aes(linewidth=weight2),alpha=0.5,curvature = 0) +
    #geom_nodes(color="black",size=13,shape=15) +
    #geom_nodes(color="white",size=11,shape=15) +
    #geom_nodes(aes(color=name),size=10,shape=15) +
    geom_nodelabel(color = "white",aes(label=Model,fill=col2),size=5,fontface="bold",
                   show.legend = F)+
    labs(color="log-likelihood ratio",y="log-likelihood")+
    scale_fill_manual(values=col_mod[levels(gig2$col2),2])+
    scale_x_continuous(breaks=seq(0,4,1))+
    theme_blank()+
    theme()+
    coord_cartesian(xlim=c(-0.2,1.2))
  
  if(!(scale)){
    gnet <- gnet+theme_classic()+
      coord_cartesian(xlim=c(min(gig2$x)-0.2,max(gig2$x)+0.2),
                      ylim=c(min(gig2$y)-100,max(gig2$y)+100))
  }
  
  if(auc){
    gnet <- ggplot(data = gig2, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(data = gig2[which(gig2$weight<0),],na.rm = T,
                 color = "blue",aes(linewidth=weight2),alpha=0.5,curvature = 0) +
      geom_edges(data = gig2[which(gig2$weight>0),],na.rm=T,
                 color = "red",aes(linewidth=weight2),alpha=0.5,curvature = 0) +
      #geom_nodes(color="black",size=13,shape=15) +
      #geom_nodes(color="white",size=11,shape=15) +
      #geom_nodes(aes(color=name),size=10,shape=15) +
      geom_nodelabel(color = "white",aes(label=Model,fill=col2),size=5,fontface="bold",
                     show.legend = F)+
      labs(color="AUC difference")+
      scale_fill_manual(values=col_mod[levels(gig2$col2),2])+
      scale_x_continuous(breaks=seq(0,4,1))+
      theme_blank()+
      theme()+
      theme_classic()+
      coord_cartesian(xlim=c(min(gig2$x)-1,max(gig2$x)+1),
                      ylim=c(min(gig2$y)-0.01,1))
  }
  
  if(log){
    gnet <- gnet+scale_y_log10()+theme_classic()+
      coord_cartesian(xlim=c(min(gig2$x)-0.2,max(gig2$x)+0.2),
                      ylim=c(min(gig2$y)-100,max(gig2$y)+100))
  }
  return(gnet)
  
}

######
##########################
ll <- readRDS(sprintf("%s/jsdm_result_in_google_colab/OTU_0.97/restricted_lls_mat_sPC90per.rds",current_dir))

a <- readRDS(sprintf("%s/jsdm_result_in_google_colab/OTU_0.97/jsdm_sdist_result_iter100_Msamp5000_sPC90per.rds",current_dir))

mod <- readRDS(sprintf("%s/jsdm_result_in_google_colab/OTU_0.97/pred_rawdata.rds",current_dir))
############################

prab <- a$data$Y

modf <- mod$Full[,c(2,4)]

ROC <- pROC::roc(pres ~ prob,
                 data = modf,
                 ci=TRUE)

png(sprintf("%s/All_ROC_OTU97.png",save.dir),width = 2000,height = 2000,res=300)
plot(ROC,　# roc()で作成したオブジェクト
     identity = TRUE,　# TRUE：AUROC=0.5となるラインを描画
     #print.thres = "best",　# "best"：ROC上の最適なカットポイントを表示
     #print.thres.best.method="closest.topleft",　# "closest.topleft"：左上隅からの距離を利用してカットオフ値を決定
     legacy.axes = TRUE,　# TRUE：x軸の目盛りを1-特異度に合わせて左が0、右が1となるように表示
     xlab="False positive rate",
     ylab="True positive rate",
     main="Cutoff sequence similarity = 97%",
     print.thres.cex  = 1.8,     #  記号の大きさを設定する（標準は１）
     cex.lab  = 1.5,       #  軸の説明の字の大きさを設定する
     cex.axis = 1.2,      #  軸の数字等（ラベル）の大きさを設定する
     cex.main = 1.7
)
text(y=0.1,x=0.2,sprintf("AUC = %s",round(ROC$auc,3)),cex=1.8)
dev.off()


######
#OTUs negativeloglike -> loglike
ll_otu <- -sapply(ll,colSums)

ll_df <- data.frame(OTU=a$species,ll_otu)

step <- rbind(c(model="Full",true_n="Full",step=4),
              cbind(model=c("NoSoil","NoBio","NoPlant","NoSpatial"),
                    true_n=c("P+Sp+Cov",
                             "P+S+Sp",
                             "S+Sp+Cov",
                             "P+S+Cov"),
                    step=3),
              cbind(model=c("PlBi","SoBi","PlSp","SoSp","PlSo","BiSp"),
                    true_n=c("P+Cov",
                             "S+Cov",
                             "P+Sp",
                             "S+Sp",
                             "P+S",
                             "Sp+Cov"),
                    step=2),
              cbind(model=c("Pl","So","Bi","Sp"),
                    true_n=c("P",
                             "S",
                             "Cov",
                             "Sp"),
                    step=1),
              c(model="Null",true_n="Null",step=0))
rownames(step) <- step[,1]


total_ll <- data.frame(model=colnames(ll_otu),
                       step=as.numeric(step[colnames(ll_otu),3]),
                       ll=colSums(ll_otu))

compfact <- list(Full=c("Bi","So","Pl","Sp"),NoBio=c("So","Pl","Sp"),
                 NoSoil=c("Bi","Pl","Sp"),NoPlant=c("Bi","So","Sp"),NoSpatial=c("Bi","So","Pl"),
                 PlSp=c("Pl","Sp"),PlBi=c("Bi","Pl"),PlSo=c("Pl","So"),
                 SoSp=c("So","Sp"),SoBi=c("So","Bi"),BiSp=c("Bi","Sp"),
                 Pl=c("Pl"),Bi=c("Bi"),So=c("So"),Sp=c("Sp"),
                 Null=c())

#############################################################3
gig <- ConvPath(total_ll,compfact)

check_pal(color,10)
col_mod <- cbind(c(unique(step[,"true_n"]),"others"),
                 c(color[c(1,3:8,10:(length(unique(step[,"true_n"]))+2))],"gray80"))
rownames(col_mod) <- col_mod[,1]
#saveRDS(col_mod,"/Volumes/Transcend/data/Sugadaira_sequence/Final_merge/analysis_series/colorlist/jsdm_model_color.rds")

drowPath(gig,col=col_mod,step=step)

ggsave(sprintf("%s/JSDM_Factors_pathway_All.png",save.dir),h=7,w=10)


gig_w <- ConvPath(total_ll,compfact,coord_weight = T)

gway <- drowPath(gig_w,col=col_mod,step=step,scale=F)

gway+labs(y="log-likelihood",x="Number of factors included in the models",
          linewidth="abs(log-likelihood ratio)")+
  theme(legend.title = element_text(hjust=0),
    axis.title.x = element_markdown(size=20),
        axis.title.y = element_markdown(size=20),
        axis.text.y = element_markdown(size=15),
        axis.text.x =element_markdown(size=15) )+
  #scale_y_continuous(labels = scales::label_comma())+
  scale_linewidth_continuous(labels = scales::label_comma())

ggsave(sprintf("%s/JSDM_Factors_pathway_All_weighted_Pll.pdf",save.dir),h=7,w=8)


gig_w$np <- ifelse(gig_w$weight<0,"Negative","Positive")

gnet <- ggplot(data = gig_w, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(data = na.omit(gig_w),na.rm = T,
             aes(color=np),linewidth=3,alpha=0.5,curvature = 0)+
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"))+
  scale_color_manual(values=c("red","blue"))

ggsave(plot=g_legend(gnet),sprintf("%s/NP_legend_pathway.pdf",save.dir),h=3,w=3)

plot(g_legend(gnet))
##
write.csv(total_ll[order(total_ll$ll),],sprintf("%s/ll_ranking.csv",save.dir))
#####################################
