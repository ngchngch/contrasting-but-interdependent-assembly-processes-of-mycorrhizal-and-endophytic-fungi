lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)


##set directory

save.dir <- sprintf("%s/01_rarecurve",current_dir)

dir.create(save.dir)


info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_root <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_rootF.rds")[[1]]
df_soil <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_soilF.rds")[[1]]
taxa <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")


taxa_f <- cbind(taxa,Genus2=NA)
taxa_f[,"Genus2"] <- apply(taxa,1,function(x){#x <- taxa_f[1,]
  if(length(grep("Unidentified",x["Genus"]))==0){sprintf("*%s*",x["Genus"])}else{
    x["Genus"]
  }})

info2 <- rowSelect(info,"target",sp_Names = "fungi")
info3 <- rowSelect(info2,"plant",sp_Names = c("Unidentified","-"),invert = T)
###3
##############
df_root2 <- df_root[which(rownames(df_root) %in% rownames(info3)),
                    which(colnames(df_root) %in% rownames(rowSelect(taxa_f,"Kingdom","Fungi")))]
df_root2 <- df_root2[,colSums(df_root2)>0]

df_soil2 <- df_soil[,which(colnames(df_soil) %in% rownames(rowSelect(taxa_f,"Kingdom","Fungi")))]
df_soil2 <- df_soil2[,colSums(df_soil2)>0]
###

length(color)
par(mar = c(6.0, 12.0, 4.1, 2))
par(mgp = c(2,1,1))
png(sprintf("%s/rarecurve_Root.png",save.dir),
    height = 3200, width = 3000, res = 600)
rarecurve(df_root2[sample(rownames(df_root2),size=200),],label=F,col=color,
                   xlab="Read counts",y="Number of fungal OTUs",
          cex.lab  = 1.2, 
          cex.axis = 1.1)

dev.off()
par(mar = c(6.0, 12.0, 4.1, 2))
par(mgp = c(2,1,1))
png(sprintf("%s/rarecurve_Soil.png",save.dir),
    height = 3200, width = 3000, res = 600)
rarecurve(df_soil2,label=F,col=color,
          cex.lab  = 1.2, 
          cex.axis = 1.1,
                   xlab="Read counts",y="Number of fungal OTUs")
dev.off()
############3
sp <- readRDS(sprintf("%s/analysis_in_linux/01_rarecurve/specaccum_object.rds",current_dir))

sp_r <- sp[["root"]]
sp_s <- sp[["soil"]]

png(sprintf("%s/spacccurve_Root.png",save.dir),
    height = 3200, width = 3000, res = 600)
par(mar = c(6.0, 6.0, 4.1, 2))
par(mgp = c(4,1.2,0))
plot(sp_r, ci.type="poly", col="blue", lwd=2, ci.lty=0,
     ci.col="lightblue",cex.lab  = 1.5, 
     cex.axis = 1.2,
     xlab="Number of Samples",ylab="Number of fungal OTUs")
dev.off()

png(sprintf("%s/spacccurve_Soil.png",save.dir),
    height = 3200, width = 3000, res = 600)
par(mar = c(6.0, 6.0, 4.1, 2))
par(mgp = c(4,1.2,0))
plot(sp_s, ci.type="poly", col="red", lwd=2, ci.lty=0,
     ci.col="orange",cex.lab  = 1.5, 
     cex.axis = 1.2,
     xlab="Number of Samples",ylab="Number of fungal OTUs")
dev.off()

png(sprintf("%s/spacccurve_Root_Soil.png",save.dir),
    height = 3200, width = 3000, res = 600)
par(mar = c(6.0, 6.0, 4.1, 2))
par(mgp = c(4,1.2,0))
plot(sp_r, ci.type="poly", col="blue", lwd=2, ci.lty=0,
     ci.col="lightblue",cex.lab  = 1.5, 
     cex.axis = 1.2,xlim=c(0,1650),y=c(0,3000),
     xlab="Number of Samples",ylab="Number of fungal OTUs")
par(new=T)
plot(sp_s, ci.type="poly", col="red", lwd=2, ci.lty=0,
     ci.col="orange",cex.lab  = 1.5, 
     cex.axis = 1.2,xlim=c(0,1650),y=c(0,3000),
     xlab="Number of Samples",ylab="Number of fungal OTUs")
dev.off()

