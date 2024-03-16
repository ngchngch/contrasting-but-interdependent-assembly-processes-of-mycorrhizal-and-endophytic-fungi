
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
raw_df <- readRDS("0_data_processing/03_make_merge_seqtab/No_rarefy_sqtb_rootF.rds")
jsdm_files <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per.rds")

input <- raw_df[rownames(jsdm_files$Occ),colnames(jsdm_files$Occ)]

mb <- spiec.easi(input, method='mb', lambda.min.ratio=1e-2,
           nlambda=20, pulsar.params=list(rep.num=50,ncores=4))

glasso <- spiec.easi(input, method='glasso', lambda.min.ratio=1e-2,
           nlambda=20, pulsar.params=list(rep.num=50,ncores=4))

slr <- spiec.easi(input, method='slr', r=0:20, lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=20, ncores=4))

bic <-R_BIC(slr,input)

slr_opt <- spiec.easi(input, method='slr', r=which.min(bic)-1, lambda.min.ratio=1e-2,
                  nlambda=20, pulsar.params=list(rep.num=50, ncores=4))

saveRDS(slr,sprintf("%s/conet_all_slr_OTU97.rds",save.dir2))
saveRDS(slr_opt,sprintf("%s/conet_opt_slr_OTU97.rds",save.dir2))
saveRDS(mb,sprintf("%s/conet_mb_OTU97.rds",save.dir2))
saveRDS(glasso,sprintf("%s/conet_glasso_OTU97.rds",save.dir2))

#OTU93
dir.create(sprintf("%s/OTU_93",save.dir))

save.dir3 <-  sprintf("%s/OTU_93",save.dir)

##
raw_df2 <- readRDS("bioinfo_output/230523_makeOTUtables/0.93/seqOTUtab_th0.93.rds")

jsdm_files2 <- readRDS("2_jSDM/01_prep_jsdm_files/jsdm_files_sPC90per_OTU93.rds")

#
input2 <- raw_df2[rownames(jsdm_files2$Occ),colnames(jsdm_files2$Occ)]

mb93 <- spiec.easi(input2, method='mb', lambda.min.ratio=1e-2,
                 nlambda=20, pulsar.params=list(rep.num=50,ncores=4))

glasso93 <- spiec.easi(input2, method='glasso', lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=50,ncores=4))


slr2 <- spiec.easi(input2, method='slr', r=0:20, lambda.min.ratio=1e-2,
                  nlambda=20, pulsar.params=list(rep.num=20, ncores=4))

bic2 <-R_BIC(slr2,input2)

slr_opt2 <- spiec.easi(input2, method='slr', r=which.min(bic2)-1, lambda.min.ratio=1e-2,
                      nlambda=20, pulsar.params=list(rep.num=50, ncores=4))

saveRDS(slr2,sprintf("%s/conet_all_slr_OTU93.rds",save.dir3))
saveRDS(slr_opt2,sprintf("%s/conet_opt_slr_OTU93.rds",save.dir3))
saveRDS(mb93,sprintf("%s/conet_mb_OTU93.rds",save.dir3))
saveRDS(glasso93,sprintf("%s/conet_glasso_OTU93.rds",save.dir3))
