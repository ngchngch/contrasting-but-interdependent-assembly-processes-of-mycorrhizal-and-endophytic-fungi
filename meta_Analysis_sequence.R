
##set directory
setwd("/Volumes/Transcend/data/Sugadaira_sequence/Final_merge/analysis_sequence")


##read original functiions
source('01_1_function.R')

##data_processing
current_dir <- "bioinfo_output"

dir.create(current_dir)

#=====ASV -> OTU97/93

source("Script/0_00_OTU_convert.R")

##data_processing
current_dir <- "0_data_processing"

dir.create(current_dir)

#=========pick_up vividplantea sequence for MEGA

source("Script/0_01_plant_seq_extract.R")

#MEGA phylogenetic tree & BLAST search

#=========Root plant anotation

source("Script/0_02_root_plant_annotation.R")

#========merge multiple sequence result

source("Script/0_03_No_rarefy_merge_seqtab.R")

raw_sqtb_path <- "0_data_processing/03_make_merge_seqtab" 

#=========covarage-based rarefaction

source("Script/0_04_covrarefy.R")

rarefy_sqtb_path <- "0_data_processing/04_rarefaction"

#==========make sample infomation with root plant

source("Script/0_05_make_sample_info.R")

#==========Guild annotation

source("Script/0_06_make_guild_table.R")

#manual guild correction

source("Script/0_07_guild_modified.R")


#==========covert to OTU cut-off threshold 0.93

source("Script/0_08_OTU_convert_0.93.R")

##basic description
current_dir <- "1_basic_discription"

dir.create(current_dir)

#==========rarefaction&species acumuration_curve
#analysis in multi core Linux PC
#caluculate species acumulation cureve

source("XXXXXXXXXXXXXXXXXX")

#drow graphs

source("Script/1_01_rarefaction_curve.R")

#===========Barplot
source("Script/1_02_barplot.R")

##joint Species distribution modeling
current_dir <- "2_jSDM"

dir.create(current_dir)

#=========prepare files for jsdm in google colab

source("Script/2_01_prep_jsdm_files.R")

#analysis in google Colab


jsdm_data_path <- sprintf("%s/jsdm_result_in_google_colab",current_dir)

#=========analysis jsdm result

# factor influencing community assembly

source("Script/2_02_sjsdm_pathway_OTU97.R")

source("Script/2_02_sjsdm_pathway_OTU93.R")

source("Script/2_02_sjsdm_pathway_pH_OTU97.R")

# factor influencing ecological process of fungal fuctional group 

source("Script/2_03_jsdm_loglikeratio_comp_OTU97.R")

source("Script/2_03_jsdm_loglikeratio_comp_OTU93.R")

##additional analysis each factor
#in this section, only OTU 0.97 data were used

current_dir <- "3_additional_analysis"

dir.create(current_dir)

#=======Mantel correlogram

#whole fungal community
source("Script/3_01_mantel_all.R")

#each guild
source("Script/3_01_mantel_AMF.R")

source("Script/3_01_mantel_EcMF.R")

source("Script/3_01_mantel_Endophyte.R")

source("Script/3_01_mantel_Mycoparasite.R")

source("Script/3_01_mantel_Nematophagous.R")

source("Script/3_01_mantel_Other_RAF.R")

source("Script/3_01_mantel_Pathogen.R")

#drow mantel correlogram panels

source("Script/3_02_mantel_plot.R")

#==========Host preference

#analysis in multi core linux PC
source("Script/3_03_2dp_dprime_randamize.R")

source("Script/3_04_plot_preference_heatmap.R")

#==========Habitat preference

source("Script/3_05_Habitat_preference_root_soil.R")

##co-occurrence network

current_dir <- "4_co-occurrence_network_slr"

dir.create(current_dir)

#=======SPIEC EASI
#analysis in multi core linux PC

source("Script/4_01_conet_SLR")

conet_data_path <- sprintf("%s/analysis_in_liunx/230620_conet_difOTU_9397",current_dir)

source("Script/4_02_conet_difOTU_slr.R")

source("Script/4_03_conet_slr_jSDM_compare.R")
