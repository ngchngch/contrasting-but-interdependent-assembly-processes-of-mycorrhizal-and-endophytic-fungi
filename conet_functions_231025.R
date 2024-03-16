#######3
R_BIC <- function(netlist,comm){sapply(netlist,function(x){
  test <- x
  ebic(getRefit(test),comm,
       test$est$loglik[getOptInd(test)],gamma = 0.5)
})}

###################
Edge_weight <- function(net,zero=F,method){
  if(method=="mb"){
    r.slr <- net
    slr.cor <- as.matrix(symBeta(getOptBeta(net), mode='ave'))
  }else{
    if(method=="slr"){
      #optimal inverse covariance -> covariance
      r.slr <- net
    }else{
      if(method=="glasso"){
        r.slr <- net
      }
    }
    
    slr.icov <- solve(as.matrix(r.slr$est$icov[[getOptInd(r.slr)]])) 
    #covariance -> correlation
    slr.cor <- cov2cor(slr.icov)
    
    #optimal network
    slr.refit <- as.matrix(getRefit(r.slr))
    sum(slr.refit)
    #extract edge weight of optimal network
    slr.cor[which(slr.refit==FALSE)] <- 0
    
  }
  diag(slr.cor) <- NA
  
  slr.cor[upper.tri(slr.cor)] <- NA
  
  dimnames(slr.cor) <- list(colnames(r.slr$est$data),
                            colnames(r.slr$est$data))
  
  gcor <- gather(cbind(OTU1=colnames(r.slr$est$data),
               as.data.frame(slr.cor)),
         OTU2,cor,-1)
  if(zero){
    gcor2 <- gcor[which(gcor$cor!=0),]  
    gcor2 <- na.omit(gcor2)
  }else{
    gcor2 <- na.omit(gcor)
  }
  
  
  return(gcor2)  
  }

###################
deg <- function(comm,r.slr,taxa,method){
  OTU <- colnames(r.slr$est$data)
  if(method=="mb"){
    slr.cor <- symBeta(getOptBeta(r.slr), mode='ave')
  }else{
    #optimal inverse covariance -> covariance
    slr.icov <- solve(as.matrix(r.slr$est$icov[[getOptInd(r.slr)]])) 
    #covariance -> correlation
    slr.cor <- cov2cor(slr.icov)
    
    #optimal network
    slr.refit <- as.matrix(getRefit(r.slr))
    sum(slr.refit)
    #extract edge weight of optimal network
    slr.cor[which(slr.refit==FALSE)] <- 0
  }
  
  diag(slr.cor) <- 0
  
  #positive network
  slr.posi <- slr.cor
  slr.posi[which(slr.posi <= 0)] <- 0
  #negativenetwork
  slr.nega <- slr.cor
  slr.nega[which(slr.nega >= 0)] <- 0
  
  ig.slr.posi <- igraph::graph.adjacency(slr.posi, mode='undirected', diag=FALSE, weighted=TRUE)
  ig.slr.nega <- igraph::graph.adjacency(-slr.nega, mode='undirected', diag=FALSE, weighted=TRUE)
  
  deg_p <- igraph::degree(ig.slr.posi)
  deg_n <- igraph::degree(ig.slr.nega)
  deg_df <- data.frame(OTU=OTU,
                       Occ=apply(comm[,OTU],2,function(x)sum(x>0)),
                       posi=deg_p,
                       nega=deg_n,
                       taxa[OTU,])
  return(deg_df)
}

###########################
drowConet <- function(net,comm,taxa_f,NP,ord=NA,curve=0,class = "Guild",method,module=T, 
                      df_root2,save.dir,name,col_list){
  
  dir.create(sprintf("%s/%s",save.dir,name))
  ##########
  if(method=="mb"){
    r.slr <- net
    slr.cor <- symBeta(getOptBeta(net), mode='ave')
  }else{
    if(method=="slr"){
      #optimal inverse covariance -> covariance
      r.slr <- net
    }else{
      if(method=="glasso"){
        r.slr <- net
        }
    }
    
    slr.icov <- solve(as.matrix(r.slr$est$icov[[getOptInd(r.slr)]])) 
    #covariance -> correlation
    slr.cor <- cov2cor(slr.icov)
    
    #optimal network
    slr.refit <- as.matrix(getRefit(r.slr))
    sum(slr.refit)
    #extract edge weight of optimal network
    slr.cor[which(slr.refit==FALSE)] <- 0
    
  }
  
  diag(slr.cor) <- 0
  
  slr.abs <- abs(slr.cor)
  slr <- slr.cor[rowSums(slr.abs)>0,colSums(slr.abs)>0]
  if(NP=="posi"){
    slr[which(slr <= 0)] <- 0  
    Ecol <- "blue"
  }else{
    if(NP=="nega"){
      slr[which(slr >= 0)] <- 0
      Ecol <- "red"
      slr <- -slr
    }
    
  }
  
  netw <- igraph::graph.adjacency(slr, mode='undirected', diag=FALSE, weighted=TRUE)
  #########################
  igraph::V(netw)$degree <- igraph::degree(netw)
  memb <- igraph::cluster_louvain(netw)$membership
  
  memb[which(memb %in% names(table(memb)[table(memb)==1]))] <- "no_module"
  
  igraph::V(netw)$mod <- memb
  
  igraph::V(netw)$OTU <- colnames(comm)[colSums(slr.abs)>0]
  
  #######################################################
  #select graph layout
  if(is.null(ncol(ord))){
    am.coord <-layout_with_stress(netw)
    rownames(am.coord) <- igraph::V(netw)$OTU
    
  }else{
      am.coord <-ord
  }
  
  g.slr1 <- ggnetwork(netw,scale=F)
  g.slr1$endOTU <- unique(g.slr1[which(g.slr1$x == g.slr1[i,"xend"] & g.slr1$y == g.slr1[i,"yend"]),"OTU"])
  
  g.slr1$endOTU <- NA
  
  for(j in 1:nrow(g.slr1)){
    g.slr1[j,"endOTU"] <- unique(g.slr1[which(g.slr1$x == g.slr1[j,"xend"] & g.slr1$y == g.slr1[j,"yend"]),
                                        "OTU"])
  }
  
  g.slr <- g.slr1
  
  g.slr$x <- am.coord[g.slr1$OTU,1]
  g.slr$y <- am.coord[g.slr1$OTU,2]
  g.slr$xend <- am.coord[g.slr1$endOTU,1]
  g.slr$yend <- am.coord[g.slr1$endOTU,2]
  
  
  
  #change here##
    #"family","genus","Order","Phylum","Guild"
  lim <- 30 #max 74
  kingdom <- "Fungi"
  node <- g.slr
  method <- "SLR"
  #####
  ##unique(taxa.ITS[,"family"])
  
  #col <- readRDS("")
  
  node$taxa <-NA
  for(i in unique(node[,"OTU"])){#i <- "1101.27569.12496__T593__S_0079"
    node[which(node[,"OTU"] %in% i), "taxa"] <- taxa_f[i,class]
  }
  
  
  node <- as.data.frame(node)
  
  colfng <- col_list[[class]]
  
  
  node$OTU2 <- factor(node$taxa,levels=c("AMF","EcMF","Endophyte",
                                          "Mycoparasite","Other_RAF","Unassigned","Unidentified")) 
  
  if(module==T){
    gnet <- ggplot(data = node, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = Ecol,alpha=0.5,curvature = curve) +
      geom_nodes(color="black",aes(size=ifelse(degree==0,0,degree*3.5)),
                 show.legend = F) +
      geom_nodes(aes(color=mod,size=ifelse(degree==0,0,degree*3)),
                 show.legend = F) +
      scale_color_manual(values=setdiff(color,colfng[levels(node$OTU2),2]))+
      ggnewscale::new_scale_color()+
      geom_nodes(color="black",aes(size=ifelse(degree==0,0,degree*1.1)),
                 show.legend = F) +
      geom_nodes(aes(color = OTU2,size=degree)) +
      scale_color_manual(values=colfng[levels(node$OTU2),2])+
      scale_size_continuous(range=c(1,max(node$degree+1)))+
      #geom_nodetext(color = "black",aes(label =OTU),size = 3, fontface = "bold") +
      #labs(title = paste("SLR_th_",sample.th,"OTU colored", class,  "cooccurence network_",plant,"r=",optR),
      #    color=class) + 
      labs(color="Functional group",size="Degree")+
      theme_blank()+
      theme(legend.text = element_markdown(size=13),
            legend.key.size=unit(1, "mm")) +
      #scale_color_manual(values=col3,
      #                  guide=guide_legend(ncol=2)) +
      guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
  }else{
    gnet <- ggplot(data = node, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = Ecol,alpha=0.5,curvature = curve) +
      geom_nodes(color="black",aes(size=ifelse(degree==0,0,degree*1.5)),
                 show.legend = F) +
      geom_nodes(aes(color = OTU2,size=degree)) +
      scale_color_manual(values=colfng[levels(node$OTU2),2])+
      scale_size_continuous(range=c(1,max(node$degree+1)))+
      #geom_nodetext(color = "black",aes(label =OTU),size = 3, fontface = "bold") +
      #labs(title = paste("SLR_th_",sample.th,"OTU colored", class,  "cooccurence network_",plant,"r=",optR),
      #    color=class) + 
      labs(color="Functional group",size="Degree")+
      theme_blank()+
      theme(legend.text = element_markdown(size=13),
            legend.key.size=unit(1, "mm")) +
      #scale_color_manual(values=col3,
      #                  guide=guide_legend(ncol=2)) +
      guides(color=guide_legend(ncol=1,override.aes = list(size=5)))
  }
 
  
  #gnet1 <- gnet + labs(title = paste(class,"level co-occurence network")) 
  gnet2 <- g_legend(ggplot(data = node, aes(x = x, y = y, xend = xend, yend = yend)) +
                      geom_edges(color = Ecol,alpha=0.5) +
                      geom_nodes(aes(color = OTU2,size=degree)) +
                      scale_color_manual(values=colfng[levels(node$OTU2),2])+
                      scale_size_continuous(range=c(4,10))+
                      #geom_nodetext(color = "black",aes(label =OTU),size = 3, fontface = "bold") +
                      #labs(title = paste("SLR_th_",sample.th,"OTU colored", class,  "cooccurence network_",plant,"r=",optR),
                      #    color=class) + 
                      labs(color="Functional group",size="Degree")+
                      theme_blank()+
                      theme(legend.text = element_markdown(size=13),
                            legend.key.size=unit(1, "mm")) +
                      #scale_color_manual(values=col3,
                      #                  guide=guide_legend(ncol=2)) +
                      guides(color=guide_legend(ncol=1,override.aes = list(size=5))))
  #gnet2 <- gnet1 + legend()
  
  
  
  ggsave(gnet, filename=sprintf("%s/%s/%s_%s_net_SpiecEasi_%s.pdf",save.dir,name,method,NP,name), h=8, w=10)
  ggsave(gnet2, filename=sprintf("%s/%s/%s_%s_net_legend_%s.pdf",save.dir,name,method,NP,name), h=8, w=10)
  
  return(list(net=gnet,data=node))
}
###########################
drow_ModsBIC <- function(bic){
  bic_df <- data.frame(rank=as.numeric(gsub("rank","",names(bic))),
                       BIC=bic,Min=ifelse(bic==min(bic),"min","" ))
  
  g <- ggplot(bic_df,aes(x=rank,y=BIC))+
    geom_line(linewidth=3,color="darkgreen")+
    geom_star(aes(starshape=Min,fill=Min,size=Min),starstroke=1,
              show.legend = F)+
    scale_starshape_manual(values=c(13,1))+
    scale_size_manual(values = c(3,6))+
    scale_fill_manual(values=c("gray10","red"))+
    labs(y="Extended BIC",x="Number of latent valuables")+
    theme_classic()+
    theme(axis.title.x = element_markdown(size=20),
          axis.title.y = element_markdown(size=20),
          axis.text = element_text(size=18))
  
  return(g)
}

#######################
#compare conet degree - llr in jsdm

Comp_llCnet <- function(slr,comm,method,#comm ==raw read teable
                        taxa_f,ll,NP="posi",R){
  
  if(NP=="posi"){
    negaposi <- "Positive"
  }else{
    if(NP=="nega"){
      negaposi <- "Negative"
    }else{
      print("NP = posi/nega")
    }
  }
  
  
  if(method=="mb"){
    dlist <- deg(comm,r.slr=slr,taxa_f,method=method)
    latent <- "mb without"
  }else{
    if(method=="slr"){
      dlist <- deg(comm,r.slr=slr,taxa_f,method=method)
      
      if(R=="rank0"){
        latent <- "without"
      }else{
        latent <- sprintf("with %s",gsub("rank","",R))
      }
      
    }else{
      if(method=="glasso"){
        dlist <- deg(comm,r.slr=slr,taxa_f,method=method)
        latent <- "glasso without"
      }
    }
  }
  
  
  
  df_llcor <- data.frame(OTU=rownames(ll),
                         llr=as.numeric(ll[rownames(ll),"llr"]),
                         taxa_f[rownames(ll),],
                         tag_rank=dlist[rownames(ll),NP],
                         max=max(dlist[rownames(ll),c("posi","nega")])
                         )
  
  r1_ct <- cor.test(df_llcor$llr,df_llcor$tag_rank,method = "kendal")
  
  df_llcor$Guild2 <- factor(df_llcor$Guild,levels=c(setdiff(unique(df_llcor$Guild),
                                                            c("Unassigned","Unidentified")),
                                                    "Unassigned","Unidentified"))
  g <- ggplot(df_llcor,aes(x=tag_rank,y=llr))+
    geom_star(aes(fill=Guild2),starshape=15,size=3)+
    
    scale_fill_manual(values=colfng_gld[levels(df_llcor$Guild2),2])+
    labs(x=sprintf("%s Degree in Co-occurrence network",negaposi),
         y="Explanatory power of OTU Covariance in jSDM\n[ log likelihood ratio (Full model vs model without OTU co-variance) ]",
        fill="Functional group"
  )+
    theme_classic()+
    theme(plot.title = element_markdown(),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text = element_text(size=15))
  return(list(net=g,cor=r1_ct))
  
}

#########################################
Comp_deg <- function(comm,net1,net2,method1,method2,#comm ==raw read teable
                           taxa_f,ll,NP="posi"){
  
  
  if(NP=="posi"){
    negaposi <- "Positive"
  }else{
    if(NP=="nega"){
       negaposi <- "Negative"
    }else{
      print("NP = posi/nega")
    }
  }
  
  dlist0 <- deg(comm,r.slr=net1,taxa_f,method=method1)
  dlist <- deg(comm,r.slr=net2,taxa_f,method=method2)
  
  
  df_llcor <- data.frame(OTU=rownames(ll),
                         llr=as.numeric(ll[rownames(ll),"llr"]),
                         taxa_f[rownames(ll),],
                         rank0=dlist0[rownames(ll),NP],
                         Optrank=dlist[rownames(ll),NP]
  )
  
  ltdeg <- data.frame(OTU=df_llcor$OTU,resid_deg=df_llcor$Optrank,
                           lat_deg=df_llcor$rank0,
                           llr=as.numeric(ll[df_llcor$OTU,"llr"]),
                           taxa_f[df_llcor$OTU,])
  
  
  ltdeg$Guild2 <- factor(ltdeg$Guild,levels=c(setdiff(unique(ltdeg$Guild),
                                                                c("Unassigned","Unidentified")),
                                                        "Unassigned","Unidentified"))
  
  g <- ggplot(ltdeg,aes(x=lat_deg,y=resid_deg))+
    geom_abline(intercept = 0,slope = 1,color="blue",alpha=0.4)+
    geom_star(aes(fill=Guild2),size=4,starshape=15)+
    
    scale_fill_manual(values=colfng_gld[levels(ltdeg$Guild2),2])+
    labs(y=sprintf("Number of %s Co-occurrence with removing latent factors\n[ Degree in %s ]",
                   negaposi,method2),
         x=sprintf("Number of %s Co-occurrence without removing latent factors\n[ Degree in %s ]",
                   negaposi,method1),
         fill="Functional group")+
    
    theme_light()+
    theme(plot.title = element_markdown(),
          legend.title=element_text(size=18),
          legend.text=element_text(size=15),
          strip.text = element_text(size=15),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text = element_text(size=15))+
    guides(fill=guide_legend(override.aes = list(size=5)))
  
  return(list(net=g,tab=ltdeg))
}


###############################
