
lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- "ggtext";library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'cowplot';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "1_basic_discription"

dir.create(current_dir)

save.dir <- sprintf("%s/02_barplot",current_dir)

dir.create(save.dir)
##

info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_root <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_rootF.rds")[[1]]
df_soil <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_soilF.rds")[[1]]
taxa <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

####
taxa_f2 <- taxa
taxa_f3 <- cbind(taxa_f2,Genus2=NA)
taxa_f3[,"Genus2"] <- apply(taxa_f2,1,function(x){#x <- taxa_f[1,]
  if(length(grep("Unidentified",x["Genus"]))==0){sprintf("*%s*",x["Genus"])}else{
    x["Genus"]
  }})

taxa_f4 <- apply(taxa_f3,2,function(x){#x <- taxa_f3[1,]
  x[grep("Unidentified",x)]<-"Unidentified"
  return(x)
})


dimnames(taxa_f4) <- dimnames(taxa_f3)

##plant count barplot
unique(info$plant)
info_sl <- info[which(!info$plant %in% c("-")),]

info_sl$palnt2 <- NA

info_sl[,"plant2"] <- ifelse(info_sl[,"plant"] == "Unidentified",
                             "Unidentified",
                             sprintf("*%s*",info_sl[,"plant"]))

df_sl <- df_root[which(rownames(df_root) %in% rownames(info_sl)),]
info_sl2 <- info_sl[rownames(df_sl),]


df_pl <- data.frame(table(info_sl2[,"plant2"]))

df_pl$x <- factor(df_pl$Var1,levels=c(setdiff(df_pl[order(df_pl$Freq,decreasing=T),"Var1"],"Unidentified"),"Unidentified"))
df_pl$x2 <- factor(df_pl$Var1,levels=c("Unidentified",setdiff(df_pl[order(df_pl$Freq),"Var1"],"Unidentified")))


gg_plant <- ggplot(df_pl[which(df_pl$Var1 != "Unidentified"),],aes(x=x2,y=Freq))+
  geom_bar(aes(color=x,fill=x),stat="identity")+
  labs(x="Plant (Genus)",y="No. of samples")+ 
  theme_light()+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_markdown(size=12),
        axis.title.x = element_text(size=15),axis.title.y = element_markdown(size=15),
        legend.position = "none"
  )+
  scale_color_manual(values=c(color[1:(length(unique(df_pl$Var1))-1)],"gray30"))+
  scale_fill_manual(values=c(color[1:(length(unique(df_pl$Var1))-1)],"gray30"))+
  coord_flip()

########
colpl <- cbind(levels(df_pl$x2),c(color[1:(length(unique(df_pl$Var1))-1)],"gray30"))
rownames(colpl) <- colpl[,1]
saveRDS(colpl,sprintf("%s/plantcolor.rds",save.dir))

#####Fungi
df_sl2 <- df_root[which(rownames(df_root) %in% rownames(info_sl[which(info_sl[,"plant2"] != "Unidentified"),])),]

df<- df_sl2/rowSums(df_sl2)

df_so <-df_soil/rowSums(df_soil)

##relative abundance of each taxa
#Root
write.csv(colSums(Taxa.mat(df,taxa_f4,"Order"))/sum(df),
          sprintf("%s/root_Order_relativeabundance.csv",save.dir))

write.csv(colSums(Taxa.mat(df,taxa_f4,"Genus"))/sum(df),
          sprintf("%s/root_Genus_relativeabundance.csv",save.dir))

gsum <- colSums(Taxa.mat(df,taxa_f4,"Guild"))/sum(df)

write.csv(c(gsum,
           Any_GUild=sum(gsum[setdiff(names(gsum),c("Unassigned","Unidentified"))]),
          Unass_Unident=sum(gsum[c("Unassigned","Unidentified")])),
          sprintf("%s/root_Guild_relativeabundance.csv",save.dir))

#Root
write.csv(colSums(Taxa.mat(df_so,taxa_f4,"Order"))/sum(df_so),
          sprintf("%s/soil_Order_relativeabundance.csv",save.dir))

write.csv(colSums(Taxa.mat(df_so,taxa_f4,"Genus"))/sum(df_so),
          sprintf("%s/soil_Genus_relativeabundance.csv",save.dir))

gsum2 <- colSums(Taxa.mat(df_so,taxa_f4,"Guild"))/sum(df_so)

write.csv(c(gsum2,
            Any_GUild=sum(gsum2[setdiff(names(gsum2),c("Unassigned","Unidentified"))]),
            Unass_Unident=sum(gsum2[c("Unassigned","Unidentified")])),
          sprintf("%s/soil_Guild_relativeabundance.csv",save.dir))


#select colors for fungaltaxa
class <- "Genus2";lim<-30
rt_gns <- setdiff(colnames(Taxa.mat(df,taxa_f4,class)),"Unidentified")[1:lim]
so_gns <- setdiff(colnames(Taxa.mat(df_so,taxa_f4,class)),"Unidentified")[1:lim]
colfng_gns <- rbind(cbind(taxa=c(union(rt_gns,so_gns),"others","Unidentified"),
                          color=c(color[1:length(union(rt_gns,so_gns))],"gray80","gray40")))

rownames(colfng_gns) <- colfng_gns[,1]
                    
class <- "Order";lim<-10
rt_ord <- setdiff(colnames(Taxa.mat(df,taxa_f4,class)),"Unidentified")[1:lim]
so_ord <- setdiff(colnames(Taxa.mat(df_so,taxa_f4,class)),"Unidentified")[1:lim]
colfng_ord <- rbind(cbind(taxa=c(union(rt_ord,so_ord),"others","Unidentified"),
                          color=c(color[1:length(union(rt_ord,so_ord))],"gray80","gray40")))

rownames(colfng_ord) <- colfng_ord[,1]

class <- "Guild"
rt_gld <- setdiff(colnames(Taxa.mat(df,taxa_f4,class)),c("Unassigned","Unidentified"))
so_gld <- setdiff(colnames(Taxa.mat(df_so,taxa_f4,class)),c("Unassigned","Unidentified"))
colfng_gld <- rbind(cbind(taxa=c(union(rt_gld,so_gld),"Unassigned","Unidentified"),
                          color=c(color[1:length(union(rt_gld,so_gld))],"gray80","gray40")))

rownames(colfng_gld) <- colfng_gld[,1]


saveRDS(list(Genus=colfng_gns,
             Order=colfng_ord,
             Guild=colfng_gld),sprintf("%s/fungi_taxa_color.rds",save.dir))


##
##
taxa_rank <- function(df,taxa,class){
  tax_tab <- Taxa.mat(df[,which(colnames(df) %in% rownames(taxa))],taxa,class)
  nam <- colnames(tax_tab)[order(colSums(tax_tab),decreasing = T)]
  return(nam)
}
ggBar <- function(df,taxa,info,class,x="sample",select_taxa=NA,lim=0,
                  spNames=NA,Unident_Names="Unidentified",
                  color,order=NA,col_selected=F){
  tax_tab <- Taxa.mat(df[,which(colnames(df) %in% rownames(taxa))],taxa,class)
  if(lim!=0){
    if(is.na(spNames[1])){
      select_taxa2 <- setdiff(select_taxa,Unident_Names)[1:lim]
      tax_tab3 <- tax_tab[,which(colnames(tax_tab) %in% select_taxa2)]
      tax_tab4 <- cbind(tax_tab3,
                        Unidentified=tax_tab[,Unident_Names],
                        others=rowSums(tax_tab[,setdiff(colnames(tax_tab),colnames(tax_tab3))]))
    }else{
      select_taxa2 <- setdiff(select_taxa,c(spNames,Unident_Names))[1:lim]
      tax_tab3 <- tax_tab[,which(colnames(tax_tab) %in% select_taxa2)]
      tax_tab4 <- cbind(tax_tab3,
                        spName=tax_tab[,spNames],
                        Unidentified=tax_tab[,Unident_Names],
                        others=rowSums(tax_tab[,setdiff(colnames(tax_tab),colnames(tax_tab3))]))
      
      colnames(tax_tab4) <- c(colnames(tax_tab3),spNames,"Unidentified","others")
    }
    
  }else{
    if(is.na(spNames[1])){
      select_taxa2 <- setdiff(select_taxa,Unident_Names)
      tax_tab3 <- tax_tab[,which(colnames(tax_tab) %in% select_taxa2)]
      tax_tab4 <- cbind(tax_tab3,
                        Unidentified=tax_tab[,Unident_Names])
    }else{
      select_taxa2 <- setdiff(select_taxa,c(spNames,Unident_Names))
      tax_tab3 <- tax_tab[,which(colnames(tax_tab) %in% select_taxa2)]
      tax_tab4 <- cbind(tax_tab3,
                        spName=tax_tab[,spNames],
                        Unidentified=tax_tab[,Unident_Names])
      
      colnames(tax_tab4) <- c(colnames(tax_tab3),spNames,"Unidentified")
    }
    
  }
  
    
    
  #apply(df_f3[setdiff(rownames(tax_tab),rownames(info)),],1,function(x){taxa_f[names(x[x>0]),]})
  
  #class <- "ASV"#"Order"#"Family"#"Genus"
  
  #tax_tab4 <- log10(tax_tab4_)*tax_tab4_/rowSums(tax_tab4_)
  gtab <- gather(cbind(info[rownames(tax_tab4),],as.data.frame(tax_tab4)),class,abundance,-c(1:(ncol(info))))
  
  if(!col_selected){
    colfng <- cbind(taxa=taxa,color=c(color[1:(length(taxa)-2)],"gray80","gray40"))
    rownames(colfng) <- colfng[,1]
  }else{
    colfng <- color
    colnames(colfng) <- c("taxa","color")
  }

  
  gtab$z <- factor(gtab$class,levels=colnames(tax_tab4))
  
  d <- vegdist(tax_tab4)
  
  a <- hclust(d)
  
  if(x=="sample"){
    
    if(!is.na(order[1])){
      gtab$ID2 <- factor(gtab$ID,levels=order)
    }else{
      nam <- rownames(tax_tab4)[a$order]
      gtab$ID2 <- factor(gtab$ID,levels=nam)
      return(nam)
    }
    
    
    gg.bar.asv <- ggplot()+
      geom_bar(data=gtab, 
               aes(x=ID2, y=abundance, fill=z,color=z),position="fill",
               stat='identity',width=0.95) +
      #scale_y_log10()+
      labs(title = paste("Top", lim, class,  "Bar plot"),
           fill = class, color = class,y = "The proportion of sequencing read counts",
           x="Root-tip samples")+ 
      theme_light()+
      theme(axis.text.x = element_blank(),axis.text.y = element_text(size=10),
            legend.text = element_markdown(),
            legend.key.size=unit(3, "mm")#,strip.background = element_blank()
            #,strip.text = element_blank()
      )+
      #facet_wrap(~sampling,scale="free")+
      coord_cartesian(ylim = c(0.04,0.98))+
      scale_fill_manual(values=c(colfng[levels(gtab$z),2]),
                        guide=guide_legend(ncol=1))+
      scale_color_manual(values=c(colfng[levels(gtab$z),2]),
                         guide=guide_legend(ncol=1))
    
  }else{
    if(x=="plant"){
      gtab$plant2 <- ifelse(gtab$plant == "Unidentified",
             "Unidentified",
             sprintf("*%s*",gtab$plant))
      
        gtab$plant3 <- factor(gtab$plant2,levels=order)
     
        
        gg.bar.asv <- ggplot()+
          geom_bar(data=gtab, 
                   aes(x=plant3, y=abundance,color=z, fill=z),position="fill",
                   stat='identity',width=0.9) +
          #scale_y_log10()+
          labs(title = paste("Top", lim, class,  "Bar plot"),
               fill = class, color = class,y = "The proportion of sequencing read counts",
               x="Plant (Genus)")+ 
          theme_light()+
          theme(axis.text.x = element_markdown(size=15),
                axis.text.y = element_markdown(size=15),
                axis.title = element_text(size=20),
                legend.text = element_markdown(),
                legend.key.size=unit(4, "mm")#,strip.background = element_blank()
                #,strip.text = element_blank()
          )+
          #facet_wrap(~sampling,scale="free")+
          coord_cartesian(ylim = c(0.04,0.98))+
          coord_flip()+
          scale_fill_manual(values=c(colfng[levels(gtab$z),2]),
                            guide=guide_legend(ncol=1))+
          scale_color_manual(values=c(colfng[levels(gtab$z),2]),
                             guide=guide_legend(ncol=1))
      }
    }
  
 
  
  
  return(gg.bar.asv)
 
}

sampleOrder <- function(df,taxa,class){
  n_class <- length(class)
  tab <- matrix(1,nrow=nrow(df),ncol=1)
  for(i in 1:n_class){
    tab <- cbind(tab,Taxa.mat(df,taxa,class[i]))
  }
  
  d <- vegdist(tab)
  
  a <- hclust(d)
  
  nam <- rownames(df)[a$order]
  return(nam)
}

####
##plant vs fungi
rtaxa_g <- taxa_rank(df,taxa_f4,"Genus2")
rtaxa_o <- taxa_rank(df,taxa_f4,"Order")
rtaxa_gld <- taxa_rank(df,taxa_f4,"Guild")

class <- "Genus2"
gns_pl <- ggBar(df,taxa_f4,x="plant",select_taxa = rtaxa_g,lim=30,
                info,class,color=colfng_gns,
                col_selected=T,order=levels(df_pl$x2))

class <- "Order"
ord_pl <- ggBar(df,taxa_f4,x="plant",select_taxa = rtaxa_o,lim=10,
                info,class,color=colfng_ord,
                col_selected=T,order=levels(df_pl$x2))

class <- "Guild"
gld_pl <- ggBar(df,taxa_f4,x="plant",select_taxa = rtaxa_gld,
                info,class,spNames="Unassigned",
                color=colfng_gld,col_selected=T,order=levels(df_pl$x2))

plot_grid(NULL,NULL,NULL,NULL,
          gns_pl+
            theme(axis.text.y = element_blank(),
                  axis.text.x = element_markdown(size=13),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size=13),
                  plot.title = element_blank(),
                  plot.margin = unit(c(0,0.1,0,0), "cm"),
                  legend.position = "none"),
          ord_pl+
            theme(axis.text.y = element_blank(),
                  axis.text.x = element_markdown(size=10),
                  axis.title.y = element_blank(),
                  plot.margin = unit(c(0,0.1,0,0), "cm"),
                  axis.title.x = element_text(size=12),
                  plot.title = element_blank(),
                  legend.position = "none"),
          gld_pl+
            theme(axis.text.y = element_blank(),
                  axis.text.x = element_markdown(size=10),
                  axis.title.y = element_blank(),
                  plot.margin = unit(c(0,0.1,0,0), "cm"),
                  axis.title.x = element_text(size=12),
                  plot.title = element_blank(),
                  legend.position = "none"),
          gg_plant+
            theme(axis.text.y = element_blank(),
                  axis.text.x = element_text(size=10),
                  axis.title.y = element_blank(),
                  plot.margin = unit(c(0,0.1,0,0), "cm"),
                  axis.title.x = element_text(size=12)),
          # plot(g_legend(gns_pl+labs(fill="Genus",color="Genus")+
          #           theme(legend.text = element_markdown(size=12)))),
          nrow=2,
          align = "hv",
          labels = c("","","","","a","b","c","d"),label_size = 24,vjust = -0.5,
          rel_heights = c(0.1,0.1,0.1,0.1,1,1,1,1),
          rel_widths = c(1,1,1,0.5,1,1,1,0.5))

ggsave(sprintf("%s/Merge_Barplot_Guild_with_Plant.png",save.dir),w=12,h=12)

#################################33
mdf <- matrix(0,
              ncol=length(union(colnames(df),colnames(df_so))),
              nrow=2)
colnames(mdf) <- union(colnames(df),colnames(df_so))

mdf[1,colnames(df)] <- colSums(df)
mdf[2,colnames(df_so)] <- colSums(df_so)


taxa_g <- taxa_rank(mdf,taxa_f4,"Genus2")
taxa_o <- taxa_rank(mdf,taxa_f4,"Order")
taxa_gld <- taxa_rank(mdf,taxa_f4,"Guild")

sample_order <- sampleOrder(df,taxa_f4,c("Order","Genus2","Guild"))

class <- "Genus2"

gns_rt <- ggBar(df,taxa_f4,info,class,lim=30,
                select_taxa = taxa_g,color=colfng_gns,
                col_selected=T,order=sample_order)
#
class <- "Order"
ord_rt <- ggBar(df,taxa_f4,info,class,select_taxa = taxa_o,
                lim=10,color=colfng_ord,col_selected=T,order=sample_order)
#

class <- "Guild"
gld_rt <- ggBar(df,taxa_f4,info,class,spNames="Unassigned",select_taxa = taxa_gld,
                color=colfng_gld,col_selected=T,order=sample_order)

#merge graphs

plot_grid(NULL,ord_rt+labs(y="The proportion of sequencing read counts\n(Order)")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  legend.text = element_markdown(size=10),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  legend.position = "none",
                  plot.title = element_blank(),
                  axis.title.x = element_blank()),
          gns_rt+labs(y="The proportion of sequencing read counts\n(Genus)")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  legend.text = element_markdown(size=10),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  legend.position = "none",
                  plot.title = element_blank(),
                  axis.title.x = element_blank()),
          
          gld_rt+labs(y="The proportion of sequencing read counts\n(Guild)")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  legend.text = element_markdown(size=10),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  legend.position = "none",
                  plot.title = element_blank()),
          rel_heights = c(0.1,1,1,1),
          #scale=c(0.9,0.9,0.7),
          nrow=4,align="hv",
          labels=c("","a","b","c"),vjust = -0.5,
          label_size = 25)

ggsave(sprintf("%s/root_barplot.png",save.dir),h=12,w=8.5)

lg_ord <- g_legend(ord_rt+theme(legend.title = element_text(size=15),
                                legend.text = element_markdown(size=12)))

ggsave(plot=lg_ord,sprintf("%s/legend_Order_rs_barplot.pdf",save.dir),h=4,w=4)

lg_gns <- g_legend(gns_rt+labs(fill="Genus",color="Genus")
                   +theme(legend.title = element_text(size=15),
                          legend.text = element_markdown(size=12))
                   +guides(fill=guide_legend(ncol=2),
                           color=guide_legend(ncol=2)))

ggsave(plot=lg_gns,sprintf("%s/legend_Genus_rs_barplot.pdf",save.dir),h=4,w=4)

lg_gld <- g_legend(gld_rt+
                     labs(fill="Guild",
                          color="Guild")+
                     theme(legend.title = element_text(size=15),
                                legend.text = element_markdown(size=12)))


ggsave(plot=lg_gld,sprintf("%s/legend_Guild_rs_barplot.pdf",save.dir),h=4,w=4)

##soil

sample_order2 <- sampleOrder(df_so,taxa_f4,c("Order","Genus2","Guild"))


#
class <- "Genus2"
gns_so <- ggBar(df_so,taxa_f4,info,class,lim=30,select_taxa = taxa_g,
                color=colfng_gns,col_selected = T,order=sample_order2)

#
class <- "Order"
ord_so <- ggBar(df_so,taxa_f4,info,class,lim=10,select_taxa = taxa_o,
                color=colfng_ord,col_selected = T,order=sample_order2)

#
class <- "Guild"
gld_so <- ggBar(df_so,taxa_f4,info,class,spNames="Unassigned",select_taxa = taxa_gld,
                color=colfng_gld,col_selected = T,order=sample_order2)

#merge graphs
plot_grid(NULL,ord_so+labs(y="The proportion of sequencing read counts\n(Order)")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  
                  legend.text = element_markdown(size=10),
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  legend.position = "none",
                  plot.title = element_blank(),
                  axis.title.x = element_blank()),
          gns_so+labs(y="The proportion of sequencing read counts\n(Genus)")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  legend.text = element_markdown(size=10),
                  legend.position = "none",
                  plot.title = element_blank(),
                  
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  axis.title.x = element_blank()),
          
          gld_so+labs(y="The proportion of sequencing read counts\n(Guild)",
                      x="Sampling positions")+
            theme(axis.text.y = element_text(size=11),
                  axis.title = element_text(size=13),
                  legend.text = element_markdown(size=10),
                  legend.position = "none",
                  
                  plot.margin = unit(c(0,0,0,0), "cm"),
                  plot.title = element_blank()),
          rel_heights = c(0.1,1,1,1),
          #scale=c(0.9,0.9,0.7),
          nrow=4,align="hv",
          labels=c("","d","e","f"),vjust = -0.5,
          label_size = 25)

ggsave(sprintf("%s/soil_barplot.png",save.dir),h=12,w=8)

