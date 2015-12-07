## This is Rscript takes raw files and filters poorly QC'ed SNPs and then reformats

library(data.table)
library(parallel)

options(stringsAsFactors=FALSE,mc.cores=8)

GRPATH<-Sys.getenv("GRPATH")
DATA_DIR<-Sys.getenv("T1DIMP_DIR")
cases<-5913
controls<-8829

files<-list.files(path=DATA_DIR,pattern="*.RData",full.names=TRUE)

t1d<-mclapply(files,function(f){
        message(paste("Processing",f))
        data.table(get(load(f)))
})

t1d<-rbindlist(t1d)
t1d.f<-subset(t1d,qc.check!='FAIL' & !is.na(p.meta) & !is.na(position))
t1d.f$details<-with(t1d.f,paste(rsid,cases,controls,sep=":"))
out<-t1d.f[,.(chromosome,position,details,p.meta),]
setnames(out,c('chr','end','name','pval'))
out$start<-as.character(as.numeric(out$end)-1)
setcolorder(out,c('chr','start','end','name','pval'))
write.table(out,file=file.path(GRPATH,'/cogs_paper/DATA/gwas/T1D.raw.bed'),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

### next we create a tabix file for further processing

#sort -k1,1 -k2,2n "$GRPATH/cogs_paper/DATA/gwas/T1D.raw.bed" > "$GRPATH/cogs_paper/DATA/gwas/T1D.bed"
#~/GIT_REPOS/htslib/bgzip "$GRPATH/cogs_paper/DATA/gwas/T1D.bed"
#~/GIT_REPOS/htslib/tabix -p bed "$GRPATH/cogs_paper/DATA/gwas/T1D.bed.gz"
#rm "$GRPATH/cogs_paper/DATA/gwas/T1D.raw.bed"
