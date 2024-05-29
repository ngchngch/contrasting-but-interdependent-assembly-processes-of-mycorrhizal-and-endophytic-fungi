#!/bin/bash

####################################################################
## 											
## ---------- 	   Denoising sequence by DADA2      ------------- ##
##
## 					2022. 02. 28. by Fujita
####################################################################

expid=XXXXXXXXXXX #analyze by each ruin

target="Fungi" ## Fungi or Plant
inputdir=${expid}/${target}/03_FilterTrimming_fastaFiles

outputdir=04_Denoising/${target}
thread=32
vsearchpath=`which vsearch`
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="04_Denoising"

####################################################################

## ============= Remove and Make directories ==================== ##

## --Remove older directory

if [ -d ${piplineID}_fastaFiles ]; then
	rm -r ${piplineID}
fi

if ls log_and_script/${piplineID}_log* >/dev/null 2>&1 ; then
	rm log_and_script/log${piplineID}
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}

## ------------------------------------------------------------- ##
## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt

####################################################################

cat <<RRR > log_and_script/script${piplineID}.R

########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
inputdir <- "$inputdir" 
outputdir <- "$outputdir"
minident=${minident}
vsearchpath="${vsearchpath}"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="${piplineID}"

start <- Sys.time()
print(start); cat("\n")
########################################################################

RRR

cat <<'RRR' >> log_and_script/script${piplineID}.R

library(dada2)
library(seqinr)

dir.create(outputdir)


takeSamplenames <- function(x, sep=''){#x=list.files(inputdir)
  
  split <- do.call(rbind, strsplit(x, sep) )
  check.unique <- sapply(apply(split, 2, table), length)
  uniquename <- which(check.unique>1)
  if(length(uniquename)==1  ) {
    samplename <- split[,uniquename]
  }else{
    samplename <- apply(split[,uniquename], 1, paste, collapse="_")
  }
  return(samplename)
}
########################################################################

# Learn forward error rates
cat( sprintf("Start learnerror from %s...\n", Sys.time()))
errF <- invisible( learnErrors(inputdir, nbases=1e8, multithread=TRUE) )

pdf(sprintf('%s/plotErrors.pdf', outputdir))
print(plotErrors(errF, nominalQ=TRUE))
dev.off()

# Infer sequence variants
derepFs <- derepFastq(inputdir, verbose = TRUE)

# Name the derep-class objects by the sample names
samplenames <- takeSamplenames(list.files(inputdir), sep='__')


names(derepFs) <- samplenames
dadaFs <- invisible( dada(derepFs, err = errF, multithread = TRUE) )

# Construct sequence table and write to disk
st.all <- makeSequenceTable(dadaFs)
saveRDS(st.all, sprintf('%s/stall_no_rmchimera_%s.rds', outputdir,expid))

##################################################################################################

##################################################################################################

finish <- Sys.time()
print(finish-start); cat("\n")
##################################################################################################
sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

RRR

Rscript log_and_script/script${piplineID}.R 2>> log_and_script/log${piplineID}.txt


