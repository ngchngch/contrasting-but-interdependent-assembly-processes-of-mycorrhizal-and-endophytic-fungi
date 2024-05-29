##read packages
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)

####

##set directory
current_dir <- "3_additional_analysis"

dir.create(current_dir)

save.dir <- sprintf("%s/01_Mantel_jac",current_dir)

dir.create(save.dir)
#======read data=======================#

info <- readRDS("0_data_processing/05_make_sample_info/comp_sample_info.rds")
df_rt <- readRDS("0_data_processing/04_rarefaction/covrarefy_sqtb_rootF.rds")[[1]]
taxa_f <- readRDS("0_data_processing/07_Guild_manual/taxa_list_mod.rds")

jsdm_files <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per.rds")

##
Oc <- df_rt[rownames(jsdm_files$Occ),]

#####

Occ <- as.matrix(Oc[rowSums(Oc)>0,colSums(Oc)>0])
Occ[Occ>0] <- 1

################

D_eco <- vegdist(Occ,method="jaccard")

D_geo <- dist(info[intersect(rownames(info),rownames(Occ)),c("x_m","y_m")])

correlog <- mantel.correlog(D_eco,D_geo,r.type="spearman",nperm=10000,mult="BH")

saveRDS(correlog,sprintf("%s/Mantel_All.rds",save.dir))

