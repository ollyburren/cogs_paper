#!/bin/bash
PPI_DIR=$GRPATH/cogs_paper/DATA/out/ppi/

for i in `\ls $PPI_DIR`; do
if [ -d "$PPI_DIR$i" ]
	then
		if [ ! -e "$PPI_DIR$i.ppi" ]
			then
			echo "Collating $i";
			## get a header
			hfile=`\ls $PPI_DIR/$i/*.ppi | head -1`;
			head -1 $hfile > $PPI_DIR$i.ppi;
			## Note that sometimes R returns standard form which causes issues with TABIX formatted files perl here should fix that
			cat $PPI_DIR/$i/*.ppi | perl -F"\t" -lane 'print join("\t",@F[0],(map{sprintf("%d",$_)}@F[1..2]),@F[3..$#F]) if $F[0] ne "chr"' | sort -k1,1 -k2,2n | sed -e 's/^chr//'   >> $PPI_DIR/$i.ppi;
	fi
fi
done
