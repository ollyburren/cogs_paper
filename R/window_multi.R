## This assigns SNPs to genes based on a window around the cannonical TSS start site for a gene
## it allows for SNPs to belong to more than one gene

library(data.table)
library(GenomicRanges)

options(stringsAsFactors=FALSE)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")

DATA_DIR<-file.path(GRPATH,'cogs_paper/DATA/')

source(file.path(GRPATH,'CHIGP/R/common.R'))

numeric.args=c('offset')

args<-list(
	offset = 2e5,
	ppi_file = file.path(DATA_DIR,'out/ppi/T1D.ppi'),
	out_file = file.path(DATA_DIR,'out/gene_score/T1D_window_multi.tab'),
	region_bed = file.path(DATA_DIR,'support/0.1cM_regions.b37.bed'),
	frag_bed = file.path(DATA_DIR,'support/Digest_Human_HindIII.bed'),
	tss_site = file.path(DATA_DIR,'support/tss.e75.transcripts.bed')
)

if(!interactive())
  args<-getArgs(verbose=TRUE,numeric=numeric.args)

prefix<-args[['prefix']]
offset<-args[['offset']]
ppi.file=args[['ppi_file']]
out.file<-args[['out_file']]
region.bed<-args[['region_bed']]
frag.bed<-args[['frag_bed']]
tss.site<-args[['tss_site']]
disease<-sub("\\.ppi$","",basename(ppi.file))


## we want the same file format as usual.


tss<-fread(tss.site)
setnames(tss,c('chr','start','end','details'))
det<-data.frame(do.call("rbind",strsplit(tss$details,":")))
names(det)<-c('gname','ensg','bt','enst','strand')
tss<-cbind(tss,det)

## split into forward/reverse genes

fwd<-subset(tss,strand=='+')
rev<-subset(tss,strand=='-')
fwd.fiveprime<-fwd[,list(tss=min(start)),by=ensg]
rev.fiveprime<-rev[,list(tss=max(start)),by=ensg]

tss.can<-rbind(fwd.fiveprime,rev.fiveprime)

setkey(tss,ensg)
setkey(tss.can,ensg)
m<-unique(tss.can[tss])

## next load in all the possible restriction frags

frags<-fread(frag.bed,header=FALSE)
fgr<-with(frags,GRanges(seqnames=Rle(V1),ranges=IRanges(start=V2,end=V3)))
mgr<-with(m,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start-offset+1,end=end+offset-1),ensg=ensg))

ol<-as.matrix(findOverlaps(mgr,fgr))

ga<-cbind(m[ol[,1]],frags[ol[,2]])
res<-ga[,list(chr=unique(chr),start=min(V2),end=max(V3)),by=ensg]

## we should be able to calculate gene scores more easily by loading in ppi directly.

test.ppi<-fread(ppi.file,header=TRUE)
ppi.gr<-with(test.ppi,GRanges(seqnames=Rle(chr),ranges=IRanges(start=end,end=end),ppi=ppi))

## next read in ld blocks

ld.blocks<-fread(region.bed,header=FALSE)
setnames(ld.blocks,c('chr','start','end','det'))
ld.gr<-with(ld.blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
ol<-as.matrix(findOverlaps(ppi.gr,ld.gr))
ppi.gr$ld<-integer(length=nrow(test.ppi))
ppi.gr[ol[,1]]$ld<-ol[,2]

## next we merge ppi with gene data

r.gr<-with(res,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ensg=ensg))
ol<-as.matrix(findOverlaps(r.gr,ppi.gr))
foo<-cbind(mcols(r.gr[ol[,1]]),mcols(ppi.gr[ol[,2]]))
bar<-data.table(as.data.frame(foo))
test<-bar[,list(sppi=sum(ppi)),by="ensg,ld"]
## final gene score involves 1-prod(1-ld.block.ppi)

final<-test[,list(gene.score=1-prod(1-sppi)),by=ensg]
final<-subset(final,!is.na(gene.score))
fi<-res[final]
det<-data.table(det)
setkey(det,ensg)
det<-unique(det)
fi<-fi[det]
fi$disease<-disease
mhc.genes<-subsetByOverlaps(r.gr,GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=35e6)))$ensg
fi$overlaps_mhc<-fi$ensg %in% mhc.genes
fi$enst<-NULL
fi$start<-NULL
fi$end<-NULL
setnames(fi,c('ensg','region_chr','geneScore','geneName','biotype','strand','disease','overlaps_mhc'))
setcolorder(fi,c('disease','ensg','geneName','biotype','strand','region_chr','geneScore','overlaps_mhc'))
#out.file<-file.path(out.dir,paste0(paste(disease,'window_multi',sep="_"),'.tab'))
write.table(fi,file=out.file,sep="\t",row.names=FALSE,quote=FALSE)






