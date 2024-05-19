
lib <- 'Matrix';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- "ggtext";library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'cowplot';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggstar';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'ggplot2';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory
current_dir <- "3_additional_analysis"

dir.create(current_dir)
save.dir <- sprintf("%s/02_plot_Mantel",current_dir)

dir.create(save.dir)
##

mg_all <- list.files(sprintf("%s/01_Mantel_jac",current_dir),
                     pattern = "Mantel_",full.names = T)

lmg_all <-lapply(mg_all,readRDS)


names(lmg_all) <- sapply(strsplit(mg_all,"_"),function(x){gsub(".rds","",x[6])})

xlist <- lmg_all
  
####
sink(sprintf("%s/mantel_results.txt",save.dir))
xlist
sink()
###

max.x <-max(sapply(xlist,function(z){max(na.omit(z$mantel.res)[,1])}))
min.y <-min(sapply(xlist,function(z){min(na.omit(z$mantel.res)[,3])}))
max.y <-max(sapply(xlist,function(z){max(na.omit(z$mantel.res)[,3])}))
  
g <- list(NULL)

for(i in 1:length(xlist)){#i <- 1
    
    lmg3 <- xlist[[i]]
    dat <- na.omit(as.data.frame(lmg3$mantel.res))
    
    dat$sig <- ifelse(dat$`Pr(corrected)`<0.05,"*P* < 0.05","N.S.")
    
    g[[i]] <- ggplot(dat,aes(y=Mantel.cor,x=class.index))+
      geom_hline(yintercept = 0,color="red")+
      geom_line()+
      geom_star(aes(starshape=sig,fill=sig),size=3,show.legend = F)+
      scale_starshape_manual(values = c(1,15))+
      scale_fill_manual(values=c("black","white"))+
      labs(x="Spatial distance (m)",
           y="Mantel's correlation",
           title=names(xlist)[i])+
      theme_classic()+
      theme(plot.title = element_text(size=19,hjust=0.5),
            axis.title = element_text(size=20),
            axis.text = element_text(size=18),
            legend.text =element_markdown(size=15),
            legend.title =element_markdown(size=18))+
      coord_cartesian(ylim=c(min.y,max.y),xlim=c(0,max.x))
  }
  
names(g) <- names(xlist)
cow_g <- plot_grid(g[["All"]]+labs(title="All fungal community")+
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_text(size=20)),
            g[["AMF"]]+labs(title="AM fungi")+
              theme(axis.title = element_blank()),
            g[["EcMF"]]+labs(title="EcM fungi")+
              theme(axis.title = element_blank()),
            g[["Endophyte"]]+labs(title="Endophyte")+
              theme(axis.title = element_blank()),
            g[["Pathogen"]]+labs(title="Plant pathogen")+
              theme(#axis.title.y = element_blank(),
                         axis.title.x = element_text(size=20)),
            g[["Mycoparasite"]]+labs(title="Mycoparasite")+
              theme(axis.title = element_blank()),
            g[["Nematophagous"]]+labs(title="Nematophagous")+
              theme(axis.title = element_blank()),
            g[["Other"]]+labs(title="Other root-associated fungi")+
              theme(axis.title = element_blank()),
            nrow=2,
            align = "hv",labels=c("a","b","c","d","e","f","g","h"),
            label_size = 19
            )
 
  
  ggsave(plot=cow_g,sprintf("%s/Mantel_merge.png",save.dir),h=9,w=18,dpi=300)
  
