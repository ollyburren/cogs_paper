#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

OUTDIR="${GRPATH}cogs_paper/DATA/out/"
PPIPATH="${OUTDIR}ppi/*.ppi"
SCOREPATH="${OUTDIR}gene_score/"
if [ ! -e $SCOREPATH ];
	then
	mkdir -p $SCOREPATH;
	echo "Creating $SCOREPATH";
fi

for i in `\ls $PPIPATH`; do 
	echo "Processing $i"
	ofile=`basename $i .ppi`.tab;
	if [ ! -e "${SCOREPATH}$ofile" ];
		then
		$RSCRIPT_BIN $GRPATH/CHIGP/R/computeGeneScore.R --pmi_file=$i --out_file=${SCOREPATH}$ofile --csnps=../DATA/RDATA/mifsud_csnps.by.ld.RData --int=../DATA/RDATA/mifsud_interactions.RData --frags=../DATA/RDATA/mifsud_frags.by.ld.RData --target.gene.cSNPs.only=1
	fi
done
