
##set directory
current_dir <- "0_data_processing"

dir.create(current_dir)

save.dir <- sprintf("%s/07_Guild_manual",current_dir)

dir.create(save.dir)



taxa <- readRDS(sprintf("%s/06_Guild_annotation/dada2_taxonomy_with_fungaltrait.rds",current_dir))

#Unasigned Genus
taxa_list <- read.csv(sprintf("%s/06_Guild_annotation/Manual/Changed_Unique_genus_guild_for_manual_check.csv",current_dir),
                      header = T)

manual_taxa <- taxa_list[which(taxa_list$Manual_guild_change=="changed"),]
rownames(manual_taxa) <- manual_taxa$Genus

taxa[which(taxa$Genus %in% manual_taxa$Genus),"Guild"] <- manual_taxa[taxa[which(taxa$Genus %in% manual_taxa$Genus),"Genus"],"Guild"]


######3
saveRDS(taxa,sprintf("%s/taxa_list_mod.rds",save.dir))
write.csv(taxa,sprintf("%s/taxa_list_mod.csv",save.dir))

