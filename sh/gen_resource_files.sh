#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

$RSCRIPT_BIN $GRPATH/CHIGP/R/generateResourceFiles.R --prefix="mifsud_" --cSNP_file="../DATA/support/cSNPs_w_ENSG.e75.bed" --interaction_file="../DATA/chic/mifsud_et_al.pm.for.chicp.25_09_2015b.tab" --pchic.thresh=5 --res_frag_file='../DATA/support/Digest_Human_HindIII.bed' --region_bed_file='../DATA/support/0.1cM_regions.b37.bed' out_dir='../DATA/RDATA/'
