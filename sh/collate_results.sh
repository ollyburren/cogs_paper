## format is different between test types and actual analysis
RESULTSDIR="${GRPATH}cogs_paper/DATA/out/gene_score/";
OUTDIR="${GRPATH}cogs_paper/DATA/out/";

for d in `\ls $RESULTSDIR/*.tab | sed -e 's/^.*\/\([^_.]*\)[_.].*/\1/' | uniq`;do
	echo "Processing $RESULTSDIR$d";
	for i in `\ls $RESULTSDIR$d*`;do
		ccount=$(head -1 $i | awk '{print NF}');
		ftype=$(basename $i .tab | sed -e 's/^[^_]*\_\(.*\)$/\1/');
		echo "Processing $ftype";
		if [ "$ftype" == "$d" ]; then
			cut -f1-7 $i |  awk  'BEGIN { OFS = "\t"} {if($1 != "disease"){print "cogs",$0}}' | sed -e 's/\.ppi//' >> "${OUTDIR}$d.tab"
		else
			awk -v ftype=$ftype 'BEGIN { OFS = "\t"} {if($8 == "FALSE"){print ftype,$0}}' $i | cut -f1-8 >> "${OUTDIR}$d.tab"
		fi
	done
done
