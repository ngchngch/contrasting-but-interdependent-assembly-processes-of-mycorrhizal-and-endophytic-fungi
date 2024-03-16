
##set directory
setwd("/Volumes/Transcend/data/Sugadaira_sequence/Final_merge/analysis_sequence")


##read original functiions
source('01_1_function.R')


current_dir <- "4_co-occurrence_network_slr"

dir.create(current_dir)


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
lib <- 'SpiecEasi';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)

lib <- 'seqinr';library(package = lib, character.only=TRUE);packageVersion(lib)

##set directory

save.dir <- sprintf("%s/01_conet_estimate",current_dir)

dir.create(save.dir)


#########3
#OTU97
dir.create(sprintf("%s/OTU_97",save.dir))

save.dir2 <-  sprintf("%s/OTU_97",save.dir)

#make input df
raw_df <- readRDS("0_data_processing/03_make_merge_seqtab/OTU97/No_rarefy_sqtb_rootF.rds")
jsdm_files <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per.rds")

input <- raw_df[rownames(jsdm_files$Occ),colnames(jsdm_files$Occ)]

slr <- spiec.easi(input, method='slr', r=0:20, lambda.min.ratio=1e-2,
                  nlambda=50, pulsar.params=list(rep.num=100, ncores=7))

bic <-R_BIC(slr,input)

slr_opt <-slr[[which.min(bic)-1]]

saveRDS(slr,sprintf("%s/conet_all_slr_OTU97.rds",save.dir2))
saveRDS(slr_opt,sprintf("%s/conet_opt_slr_OTU97.rds",save.dir2))

#OTU93
dir.create(sprintf("%s/OTU_93",save.dir))

save.dir3 <-  sprintf("%s/OTU_93",save.dir)

##

raw_df2 <- readRDS("0_data_processing/03_make_merge_seqtab/OTU93/No_rarefy_sqtb_rootF.rds")
jsdm_files2 <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per_OTU93.rds")

#
input2 <- raw_df2[rownames(jsdm_files2$Occ),colnames(jsdm_files2$Occ)]

slr2 <- spiec.easi(input2, method='slr', r=0:20, lambda.min.ratio=1e-2,
                   nlambda=50, pulsar.params=list(rep.num=100, ncores=7))

bic2 <-R_BIC(slr2,input2)


saveRDS(slr2,sprintf("%s/conet_all_slr_OTU93.rds",save.dir3))
saveRDS(slr_opt2,sprintf("%s/conet_opt_slr_OTU93.rds",save.dir3))
