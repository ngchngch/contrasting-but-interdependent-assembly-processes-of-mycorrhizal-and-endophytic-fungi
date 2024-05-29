
##set directory
current_dir <- "0_data_processing"

dir.create(current_dir)

save.dir <- sprintf("%s/06_Guild_annotation",current_dir)

dir.create(save.dir)

##read guild
raw_gld <- read.csv("Raw_data/metadata/Fungaltrait_modified_230929.csv")

taxa <- readRDS("data_process/06_Taxonomy_annotaion/fungi_annnotaion_results_dada2.rds")
#apply(raw_gld,2,unique)

prim <- c("mycoparasite","root_endophyte","ectomycorrhizal","arbuscular_mycorrhizal")
pr_pp <- c("plant_pathogen")

#exclude ErM necause No host sample 
sec <- c("mycoparasite","root-associated","root_endophyte",
         "root_endophyte_dark_septate",
         "nematophagous","arbuscular_mycorrhizal")

sc_pp <- c("plant_pathogen")

rt_pp <- c("root_pathogen")

pr_ap <- c("animal_parasite")

sc_ap <- c("animal_parasite")

rt_ap <- c("nematophagous")

##modefied fungaltrait for root fungi
gld <- raw_gld
gld[which(!raw_gld$primary_lifestyle %in% prim),"primary_lifestyle"] <- NA
gld[which(raw_gld$primary_lifestyle %in% pr_pp & raw_gld$Plant_pathogenic_capacity_template %in% rt_pp),
    "primary_lifestyle"] <- "root_pathogen"

gld[which(raw_gld$primary_lifestyle %in% pr_ap & raw_gld$Animal_biotrophic_capacity_template %in% rt_ap),
    "primary_lifestyle"] <- "nematophagous"


gld[which(!raw_gld$Secondary_lifestyle %in% sec),"Secondary_lifestyle"] <- NA
gld[which(raw_gld$Secondary_lifestyle %in% sc_pp & raw_gld$Plant_pathogenic_capacity_template %in% rt_pp),
    "Secondary_lifestyle"] <- "root_pathogen"


gld[which(raw_gld$Secondary_lifestyle %in% sc_ap & raw_gld$Animal_biotrophic_capacity_template %in% rt_ap),
    "Secondary_lifestyle"] <- "nematophagous"

gld[which(gld$Secondary_lifestyle == "root_endophyte_dark_septate" & 
            gld$primary_lifestyle == "root_endophyte"),"primary_lifestyle"] <- "root_endophyte_dark_septate"
gld_tab <- unique(cbind(gld[,1:5],
                 Guild=ifelse(!is.na(gld[,"primary_lifestyle"]),
                              gld[,"primary_lifestyle"],
                              ifelse(!is.na(gld[,"Secondary_lifestyle"]),
                                     gld[,"Secondary_lifestyle"],NA))))


gld_tab2 <- gld_tab[!is.na(gld_tab$Guild),]
rownames(gld_tab2) <- gld_tab2$GENUS
       

    
taxa_g <- as.data.frame(taxa)

gld2 <-ifelse(taxa_g$Genus %in% gld_tab2$GENUS,
                              gld_tab2[taxa_g$Genus,"Guild"],
                              "Unassigned")
taxa_g$Guild <- ifelse(str_detect(taxa_g$Genus,pattern = "Unidentified|Incertae_sedis"),
                       "Unidentified",gld2)

table(taxa_g$Guild)

#rename guild
taxa_g2 <- taxa_g

taxa_g2[which(taxa_g$Guild %in% c("root_endophyte","root_endophyte_dark_septate")),"Guild"] <- "Endophyte"
taxa_g2[grep("root_pathogen",taxa_g$Guild),"Guild"] <- "Pathogen"
taxa_g2[grep("nematophagous",taxa_g$Guild),"Guild"] <- "Nematophagous"
taxa_g2[grep("root-associated",taxa_g$Guild),"Guild"] <- "Other_RAF"
taxa_g2[grep("mycoparasite",taxa_g$Guild),"Guild"] <- "Mycoparasite"
taxa_g2[grep("ectomycorrhizal",taxa_g$Guild),"Guild"] <- "EcMF"

taxa_g2[grep("arbuscular_mycorrhizal",taxa_g$Guild),"Guild"] <- "AMF"
taxa_g2[which(taxa_g$Phylum == "Glomeromycota"),"Guild"] <- "AMF"

table(taxa_g2$Guild)


saveRDS(taxa_g2,sprintf("%s/dada2_taxonomy_with_fungaltrait.rds",save.dir))
write.csv(taxa_g2,sprintf("%s/dada2_taxonomy_with_fungaltrait.csv",save.dir))

########
write.csv(unique(taxa_g2[which(taxa_g2$Kingdom=="Fungi"),c("Kingdom","Phylum","Class","Order","Family","Genus","Guild")]),sprintf("%s/Unique_genus_guild_for_manual_check.csv",save.dir))

