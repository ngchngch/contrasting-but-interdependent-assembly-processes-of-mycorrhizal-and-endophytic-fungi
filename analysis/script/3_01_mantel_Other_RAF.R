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

tax_sl <- taxa_f[colnames(Occ),"Guild"]
mj_gld <- setdiff(names(table(tax_sl))[table(tax_sl)>5],c("Unassigned","Unidentified"))

print(mj_gld)

###
Guild <- "Other_RAF"

otu <-rownames(taxa_f[which(taxa_f[,"Guild"] == Guild),])
df2 <- Occ[intersect(rownames(info),rownames(Occ)),which(colnames(Occ) %in% otu)]
df2 <- df2[rowSums(df2)>0,colSums(df2)>0]

D_eco <- vegdist(df2,method="jaccard")

D_geo <- dist(info[intersect(rownames(info),rownames(df2)),c("x_m","y_m")])

correlog <- mantel.correlog(D_eco,D_geo,r.type="spearman",nperm=10000,mult="BH")

saveRDS(correlog,sprintf("%s/Mantel_%s.rds",save.dir,Guild))

