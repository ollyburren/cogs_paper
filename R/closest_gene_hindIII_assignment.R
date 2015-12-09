library(GenomicRanges)
library(data.table)


options(stringsAsFactors=FALSE)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")

DATA_DIR<-file.path(GRPATH,'cogs_paper/DATA/')

source(file.path(GRPATH,'CHIGP/R/common.R'))


args<-list(
  h3.file = file.path(DATA_DIR,'support/Digest_Human_HindIII.bed'),
  gene.file = file.path(DATA_DIR,'support/Homo_sapiens.GRCh37.75.genes.gtf'),
  out.file =  file.path(DATA_DIR,'/RDATA/hindIII_gene_assignment.RData')
)

if(!interactive())
  args<-getArgs(verbose=TRUE)


h3<-fread(args[['h3.file']],header=FALSE)
setnames(h3,c('chr','start','end','id'))
## get rid of y chromosome
h3<-subset(h3,chr!='Y')

genes<-fread(args[['gene.file']],header=FALSE)
setnames(genes,c('chr','bt','type','start','end','phase','strand','score','details'))
genes<-subset(genes,chr %in% unique(h3$chr))
genes$ensg<-gsub(".*(ENSG[0-9]+).*","\\1",genes$details)
genes$geneName<-gsub('.*gene_name "([^"]+)".*',"\\1",genes$details)

det<-genes[,.(geneName,ensg,chr,bt,strand)]
setnames(det,"strand","xstrand")

g.gr<-with(genes,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),ensg=ensg))
h.gr<-with(h3,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=id))
seqlevels(g.gr)<-seqlevels(h.gr)

gl<-sapply(split(end(h.gr),seqnames(h.gr)),max)
s<-Seqinfo(seqnames=names(gl),seqlength=gl,isCircular=rep(FALSE,length(gl)),genome="GRCh37")
seqinfo(h.gr)<-s
seqinfo(g.gr)<-s

## first find h3 that overlap actual genes and assign

ol<-mergeByOverlaps(h.gr,g.gr)
## find the fragments that overlap more than one gene
doubs<-unique(ol[duplicated(ol$id),]$id)
idx<-which(ol$id %in% doubs)
doubs<-ol[idx,]
## for doubs we randomly select a gene to assign to a fragment
spdoubs<-split(doubs$ensg,doubs$id)
d<-sapply(seq_along(spdoubs),function(i){
	ret<-sample(spdoubs[[i]],1)
	names(ret)<-names(spdoubs)[i]
	ret
})

sings<-ol[-idx,]
tmp<-sings$ensg
names(tmp)<-sings$id
olgenes<-c(tmp,d)
olgenes<-olgenes[order(as.numeric(names(olgenes)))]


intergenic<-gaps(g.gr)
intergenic<-intergenic[strand(intergenic)=="*",]
intergenic$iid<-1:length(intergenic)

## divide each of these by two

lower<-intergenic
upper<-intergenic
mid<-intergenic
## divide each into two and compute the middle region

end(lower)<-end(lower)-floor(width(lower)/2)
start(upper)<-start(upper) + ceiling(width(intergenic)/2)
start(mid)<-end(lower)
end(mid)<-end(lower)
nd<-h.gr[!h.gr$id %in% names(olgenes),]
## work out which fragments overlap midpoints

both.gr<-subsetByOverlaps(nd,mid)

upperb<-subsetByOverlaps(upper,both.gr)
lowerb<-subsetByOverlaps(lower,both.gr)
upperb$iid<-paste0('U',upperb$iid)
lowerb$iid<-paste0('L',lowerb$iid)

t.gr<-c(upperb,lowerb)

cv<-coverage(both.gr)
t.grl<-split(t.gr,seqnames(t.gr))
prop<-do.call('c',lapply(seq_along(cv),function(i){
	seqname<-names(cv)[i]
	t<-Views(cv[[seqname]],start=start(t.grl[[seqname]]), end=end(t.grl[[seqname]]), names=t.grl[[seqname]]$iid)
	p<-viewSums(t)/width(t)
	names(p)<-names(t)
	p
}))

u<-prop[grep("^U",names(prop))]
names(u)<-sub("^U","",names(u))
l<-prop[grep("^L",names(prop))]
names(l)<-sub("^L","",names(l))
foo<-merge(l,u,by.x=0,by.y=0)
row.names(foo)<-foo$Row.names
foo$Row.names<-NULL
names(foo)<-c('l','u')
pick.lower<-row.names(foo[foo$l>=foo$u,])

lowerb<-subset(lowerb,iid %in% paste0('L',pick.lower))

lowerolb<-mergeByOverlaps(lowerb,both.gr)


pick.upper<-row.names(foo[foo$l<foo$u,])

upperb<-subset(upperb,iid %in% paste0('U',pick.upper))

upperolb<-mergeByOverlaps(upperb,both.gr)

## remove these asignations from the remaining

#lower<-subset(lower,!iid %in% as.numeric(rownames(foo)))
#upper<-subset(upper,!iid %in% rownames(foo))

## next compute uppers and lowers



## remove these from further calculations
nd<-subset(nd,!id %in% union(upperolb$id,lowerolb$id))

lower.h<-mergeByOverlaps(lower,nd)
##add in frags that have frags that span midpoint that we have assigned
names(lowerolb)<-names(lower.h)
lower.f<-rbind(lower.h,lowerolb)


upper.h<-mergeByOverlaps(upper,nd)
names(upperolb)<-names(upper.h)
upper.f<-rbind(upper.h,upperolb)

## add gene annotations

ol<-as.matrix(findOverlaps(g.gr,lower.f$lower,maxgap=1L))
lower.f$ensg<-character(length=nrow(lower.f))
lower.f[ol[,2],]$ensg<-g.gr[ol[,1],]$ensg

lowergenes<-lower.f$ensg
names(lowergenes)<-lower.f$id

ol<-as.matrix(findOverlaps(g.gr,upper.f$upper,maxgap=1L))
upper.f$ensg<-character(length=nrow(upper.f))
upper.f[ol[,2],]$ensg<-g.gr[ol[,1],]$ensg

uppergenes<-upper.f$ensg
names(uppergenes)<-upper.f$id

final<-c(uppergenes,lowergenes,olgenes)

final<-final[order(as.numeric(names(final)))]

final.l<-split(final,names(final))

h.gr$ensg<-unlist(final.l[as.character(h.gr$id)])
h.gr<-subset(h.gr,ensg!="")
m<-as.data.frame(mcols(h.gr))
m$order<-1:nrow(m)
fob<-merge(m,det,by.x="ensg",by.y="ensg",all.x=TRUE)
fob<-fob[order(fob$order),]
mcols(h.gr)<-DataFrame(fob)

## using the above gr we can for each pmi assign SNPs to a gene and compute gene score

save(h.gr,file=args[['out.file']])








