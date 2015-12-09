#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`;
WINDOW_SIZE=2e5;
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

OUTDIR="${GRPATH}cogs_paper/DATA/out/"
PPIPATH="${OUTDIR}ppi/*.ppi"
SCOREPATH="${OUTDIR}gene_score/"
SUPPORTPATH="${GRPATH}cogs_paper/DATA/support/"
if [ ! -e $SCOREPATH ];
	then
	mkdir -p $SCOREPATH;
	echo "Creating $SCOREPATH";
fi

for i in `\ls $PPIPATH`; do 
	echo "Processing $i"
	disease=`basename $i .ppi`;
	outfile="$SCOREPATH${disease}_closest_single.tab"
	if [ ! -e "$outfile" ];
		then
		$RSCRIPT_BIN $GRPATH/cogs_paper/R/closest_gene.R  --ppi_file=$i --out_file=$outfile --region_bed=$SUPPORTPATH/0.1cM_regions.b37.bed --frags=${GRPATH}cogs_paper/DATA/RDATA/hindIII_gene_assignment.RData
	fi
done
