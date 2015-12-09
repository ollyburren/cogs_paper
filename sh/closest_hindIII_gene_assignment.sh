#!/bin/bash
## get location of Rscript
RSCRIPT_BIN=`which Rscript`;
WINDOW_SIZE=2e5;
if [ ! -e $RSCRIPT_BIN  ]
 then
	echo "Cannot find Rscript aborting\n"
fi

DATADIR="${GRPATH}cogs_paper/DATA/";
SUPPORTPATH="${DATADIR}/support/";
outfile="${DATADIR}RDATA/hindIII_gene_assignment.RData";


if [ ! -e "$outfile" ];
	then
	$RSCRIPT_BIN $GRPATH/cogs_paper/R/closest_gene_hindIII_assignment.R --h3.file=${SUPPORTPATH}Digest_Human_HindIII.bed --gene.file=${SUPPORTPATH}Homo_sapiens.GRCh37.75.genes.gtf --out.file=$outfile
fi

