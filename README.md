# RootFungalAssembly
Scripts for

"Mycorrhizal and endophytic fungi structure contrasting but
interdependent assembly processes"

Mikihito Noguchi1,2 and Hirokazu Toju3,4

1Center for Ecological Research, Kyoto University, Otsu, Shiga 520-2133, Japan

2Research Fellow of Japan Society for the Promotion of Science

3Laboratory of Ecosystems and Coevolution, Graduate School of
Biostudies, Kyoto University, Kyoto 606-8501, Japan

4Center for Living Systems Information Science (CeLiSIS), Graduate
School of Biostudies, Kyoto University, Kyoto 606-8501, Japan

Correspondence:
Mikihito Noguchi: noguchi.mikihito.57s@st.kyoto-u.ac.jp
Hirokazu Toju: toju.hirokazu.4c@kyoto-u.ac.jp

DOI: https://biorxiv.org/cgi/content/short/2024.02.17.580831v1

# Content
at first, demultiplexed fastq files were convert to sample-OTU matrix & OTU taxonomy annotation information in data_process (script in data_process/script)

the processed data were analyzed in anlysis (script in analysis/script)

Original functions in function

Some output of analyses are in Output
##Fungal OTU-sample matrix
#sequance similarity threshold: 93%
Output/Fungi/seqOTUtab_93.rds
#sequance similarity threshold: 97%
Output/Fungi/seqOTUtab_97.rds

##Plant OTU-sample matrix
Output/Plant/seqOTUtab.rds

##Fungal OTU taxonomy annotation with manual modified guild information
Output/Fungi/taxa_list_mod.csv

##Plant OTU taxonomy annotation
Output/Plant/OTUseq_0.97.5nn.tsv

##sample information with root plant annotation
Output/sample_information/comp_sample_info.csv


Metadata in Raw_data/metadata

In detail, 

######################################################################

@directory data_process
@@below scripts are in data_process/script

############ processed by each run ###################################

##adaptor trimming 
#Fungi
script02_Cutadaptor_ITS_3p62plF1_ITS_4unR1.sh
#Plant
script02_Cutadaptor_ITS1F-Kyo1_ITS2-KYO2.sh

##Filtering & Trimming by dada2
script03_FilterTrimming.sh

##Denoiding by dada2
script04_Denoising.sh

#######################################################################
## Decontamination > Merge run data > OTU clustering
231007_sample_process.R

## Taxonomy annotation
#Fungi 
dada2_annotation.R
#Plant
x_nn_BLAST.sh

#######################################################################

@directory anlysis
@@below scripts are in anlysis/script

#######################################################################
###analysis_00_data_process

##pick_up vividplantea sequence for MEGA
0_01_plant_seq_extract.R
#MEGA phylogenetic tree & BLAST search > relurt in 01_extract_pseq/Manual/Plant_seq_in_ITS1_BLAST.csv

##Root plant anotation
0_02_root_plant_annotation.R

##separate samples by root vs soil
0_03_No_rarefy_merge_seqtab.R

##covarage-based rarefaction
0_04_covrarefy.R

##make sample information with root plant
0_05_make_sample_info.R

##Guild annotation
0_06_make_guild_table.R
#manual modified guild annotation in 06_Guild_annotation/Manual/Changed_Unique_genus_guild_for_manual_check.csv

##manual guild correction
0_07_guild_modified.R

###analysis_01_basic description

##rarefaction & species acumulation cureve
1_01_rarefaction_curve.R

##Barplot
1_02_barplot.R


###analysis_02_joint species distribution modeling

##prepare files for jsdm in google colab
2_01_prep_jsdm_files.R

##jSDM
#parametor fitting & model comparison in google Colab
2_sjSDM.ipynb

#factor influencing community assembly
2_02_sjsdm_pathway_OTU97.R
2_02_sjsdm_pathway_OTU93.R
2_02_sjsdm_pathway_pH_OTU97.R

#factor influencing ecological process of fungal fuctional group 
2_03_jsdm_loglikeratio_comp_OTU97.R
2_03_jsdm_loglikeratio_comp_OTU93.R

###additional analysis each factor
##in this section, only OTU 0.97 data were used

##Mantel correlogram
#whole fungal community
3_01_mantel_all.R

#each guild
3_01_mantel_EcMF.R
3_01_mantel_Endophyte.R
3_01_mantel_Mycoparasite.R
3_01_mantel_Nematophagous.R
3_01_mantel_Other_RAF.R
3_01_mantel_Pathogen.R

#drow mantel correlogram panels
3_02_mantel_plot.R

##Host preference
#calculate preference
3_03_2dp_dprime_randamize.R

#graphics
3_04_plot_preference_heatmap.R

##Habitat preference
3_05_Habitat_preference_root_soil.R

###co-occurrence network
#co-occurence inference by SPIECEASI
4_01_conet_SLR

#graphics & additional anlysis
4_02_conet_difOTU_slr.R
