## assigns SNPs to genes based on proximity and computes a gene score. where genes overlap uses random
## sampling to assign. Note that allows only a 1-1 mapping between a SNP and 
## a gene.

## depends on ./closest_gene_hindIII_assignment.R

library(data.table)
library(GenomicRanges)


options(stringsAsFactors=FALSE)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")

DATA_DIR<-file.path(GRPATH,'cogs_paper/DATA/')

source(file.path(GRPATH,'CHIGP/R/common.R'))

args<-list(
	ppi_file = file.path(DATA_DIR,'out/ppi/T1D.ppi'),
	out_file = file.path(DATA_DIR,'out/gene_score/T1D_closest_single.tab'),
	region_bed = file.path(DATA_DIR,'support/0.1cM_regions.b37.bed'),
	frags = file.path(DATA_DIR,'RDATA/hindIII_gene_assignment.RData')
)

if(!interactive())
  args<-getArgs(verbose=TRUE)

ppi.file=args[['ppi_file']]
out.file<-args[['out_file']]
region.bed<-args[['region_bed']]
frags<-args[['frags']]
disease<-sub("\\.ppi$","",basename(ppi.file))


## we want the same file format as usual.


## next load in all the possible restriction frags

fgr<-get(load(frags))


## we should be able to calculate gene scores more easily by loading in ppi directly.

test.ppi<-fread(ppi.file,header=TRUE)
ppi.gr<-with(test.ppi,GRanges(seqnames=Rle(chr),ranges=IRanges(start=end,end=end),ppi=ppi,rs=rsid))

## remove snps that overlap the MHC region

mhc.snps<-subsetByOverlaps(ppi.gr,GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=35e6)))$rs
ppi.gr<-subset(ppi.gr,!rs %in% mhc.snps)
#test.ppi<-subset(test.ppi,!rs %in% mhc.snps)

## next read in ld blocks

ld.blocks<-fread(region.bed,header=FALSE)
setnames(ld.blocks,c('chr','start','end','det'))
ld.gr<-with(ld.blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ld.id=1:nrow(ld.blocks)))
ol<-as.matrix(findOverlaps(ppi.gr,ld.gr))
ppi.gr$ld<-integer(length=length(ppi.gr))
ppi.gr[ol[,1]]$ld<-ld.gr[ol[,2],]$ld.id

## next we merge ppi with gene data

foo<-mergeByOverlaps(ppi.gr,fgr)

test<-data.table(as.data.frame(foo[,c('ppi','ld','rs','ensg')]))

test<-test[,list(sppi=sum(ppi)),by="ensg,ld"]

## next for each 

## final gene score involves 1-prod(1-ld.block.ppi)

final<-test[,list(gene.score=1-prod(1-sppi)),by=ensg]
final<-subset(final,!is.na(gene.score))
det<-data.table(as.data.table(mcols(h.gr)))
setkey(det,ensg)
det<-unique(det)
fi<-det[final]
fi$overlaps_mhc<-FALSE
fi$disease<-disease
fi$order<-NULL
fi$id<-NULL


setnames(fi,c('ensg','geneName','region_chr','biotype','strand','geneScore','overlaps_mhc','disease'))
setcolorder(fi,c('disease','ensg','geneName','biotype','strand','region_chr','geneScore','overlaps_mhc'))
write.table(fi,file=out.file,sep="\t",row.names=FALSE,quote=FALSE)
print(paste("Success",out.file,"written"))






