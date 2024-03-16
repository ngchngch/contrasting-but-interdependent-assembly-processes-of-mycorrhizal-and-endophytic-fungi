##read packages
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'rgdal';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

save.dir <- sprintf("%s/05_make_sample_info",current_dir)

dir.create(save.dir)

#=read_data=============#
sinfo1 <- readRDS(sprintf("%s/02_root_plant_annotation/sample_info_with_plant.rds",current_dir))

sinfo <- sinfo1[grep("seedling",sinfo1[,"sampling"],invert=T),]
taxa_f <- readRDS("bioinfo_output/230414_dada2annotation/fungi_annnotaion_results_dada2.rds")


site <- read.csv("metadata/sample_point_info.csv",row.names=3,header=FALSE)

#====sampling point coordinate===========#

colnames(site) <- c("x","y")
rownames(site) <- gsub(" ","",rownames(site))
#plot(site[,1],site[,2])

info3_5 <- sinfo[-union(grep("-",sinfo[,"site"]),grep("BLANK|NegaCon",sinfo[,"site"])),]
info4 <- cbind(info3_5,site2=NA)

info4[,"site2"] <- paste("p",formatC(as.numeric(info3_5[,"site"]),width=3,flag="0"),sep="")

xy <- matrix(unlist(apply(info4,1,function(x){#x <- info4[1,]
  y <- x["site2"]
  site[which(rownames(site) %in% y),]
})),ncol=2,byrow=T)

colnames(xy) <- c("x","y")
rownames(xy) <- rownames(info4)
info5 <- data.frame(info4,xy)
rownames(info5) <- info5[,"ID"]

#####

sp2 <- SpatialPoints(info5[,c("y","x")],proj4string=CRS("+proj=longlat +datum=WGS84"))

a <- spTransform(sp2, CRS=CRS("+init=epsg:2451"))
info6 <- data.frame(info5,matrix(as.numeric(formatC(a@coords,format="f",digits=3)),ncol=2,dimnames=list(rownames(info5),c("y_m","x_m"))))

info7 <- info6[which(info6[,"plant"] !="contami"),]
write.csv(info7,sprintf("%s/comp_sample_info.csv",save.dir))

saveRDS(info7,sprintf("%s/comp_sample_info.rds",save.dir))
