library(data.table)
library(ggplot2)
library(hexbin)
library(magrittr)
library(reshape2)

options(stringsAsFactors=FALSE)

if(!interactive())
  GRPATH<-Sys.getenv("GRPATH")

DATA_DIR<-file.path(GRPATH,'cogs_paper/DATA/')
PLOT_DIR<-file.path(GRPATH,'cogs_paper/DATA/plots/')
#source(file.path(GRPATH,'cogs_paper/R/enrichment_functions.R'))

load.f<-function(f){
  fi<-fread(f,header=FALSE)
  setnames(fi,c('method','disease','ensg','geneName','biotype','strand','chr','score'))
  fi
}


extract.genesets <- function(str,myid,gmt,minsize=10) {
  gmt.sub<-subset(gmt,grepl(str,gmt$set))
  nm <- gmt.sub[,"set"]
  genes <- gmt.sub[,"genes"] %>%
    strsplit(.,split="&") %>% lapply(.,unique) %>% lapply(.,intersect,myid)
  use <- sapply(genes,length) > minsize
  structure( genes[ use ], names=nm[ use ])
}

valtest_genesetlist <- function(test.id, test.val, setlist) {
  use <- test.id %in% unlist(setlist)
  tmp <- lapply(setlist, function(set) valtest_geneset(test.id[use],test.val[use],set))
  names(tmp) <- names(setlist)
  return(tmp)
}
valtest_geneset<-function(test.id,test.val,set){
  ## message("test: ",length(test))
  ## message("set: ",length(set))
  ## message("universe: ",length(universe))
  ## m <- sum(universe %in% set)
  ## n <- length(universe) - m
  ## k <- sum(test %in% universe)
  ## q <- sum(universe %in% test & universe %in% set)
  inset <- which(test.id %in% set)
  kk <- t.test(test.val[inset], test.val[-inset])
  qnorm(kk$p.value/2,lower.tail=FALSE) * sign(kk$statistic)
}

plot_bc<-function(r){
  lapply(seq_along(r),function(i){
    x<-r[[i]]
    P <- abs(x) %>% pnorm(., lower.tail=FALSE)*2 %>% p.adjust()
    thr <- min(abs(x)[ P<0.01 ])
    use <- apply(x>=thr,1,any)
    x <- x[rownames(x) %in% names(use)[use],]
    x<-melt(x)
    names(x)<-c('pathway','method','zscore')
    x$pathway<-gsub("HALLMARK\\_","",x$pathway)
    ggplot(x,aes(x=pathway,y=zscore,fill=method)) + geom_bar(stat = "identity", position="dodge") + geom_abline(intercept=thr,slope=0) + coord_flip() + ggtitle(names(r)[i]) + theme_bw()
  })
}

gsea <- function(gene.scores,geneset="hall",...) {
  ## pathwayshallmark.ensg.msigdb.v5.0.symbols.gmt
  gmt_file<-file.path(DATA_DIR,"support/ensg.msigdb.v5.0.symbols.gmt")
  gmt<-read.delim(gmt_file,sep="\t",stringsAsFactors=FALSE)
  names(gmt)<-c('set','url','genes')
  gmt.lists <- lapply(list(hall="^HALLMARK",
                           reactome="^REACTOME",
                           c1="^C1"),
                      extract.genesets,unique(gene.scores$ensg),gmt)
  ## for each type compute the enrichment scores over all hallmark pathways
  tests<-names(gene.scores)
  tests<-tests[!tests %in% c('ensg','disease')]
  RESULTS <- matrix(NA,length(gmt.lists[[geneset]]),length(tests),
                    dimnames=list(names(gmt.lists[[geneset]]),tests))
  for(t in tests) {
    RESULTS[,t] <- with(gene.scores,
                          valtest_genesetlist(ensg, rank(gene.scores[[t]]), gmt.lists[[geneset]]) %>% unlist())
    
  }
  RESULTS
}

## get files 

fls<-list.files(path=file.path(DATA_DIR,'out'),pattern="*.tab",full.names = TRUE)

dt<-lapply(fls,load.f)

## what about ragged lists where we include everything ?

gmt_file<-file.path(DATA_DIR,"support/ensg.msigdb.v5.0.symbols.gmt")
gmt<-read.delim(gmt_file,sep="\t",stringsAsFactors=FALSE)
names(gmt)<-c('set','url','genes')


thresh<-0.5

pway_set<-'^KEGG'

non.filt<-do.call("rbind",lapply(dt,function(d){
  tmp<-subset(d,biotype=="protein_coding")
  t.results<-split(tmp,tmp$method)
  re<-do.call("rbind",lapply(t.results,function(t){
    thresh.t<-subset(t,score>=thresh)
    hmark<-extract.genesets(pway_set,unique(t$ensg),gmt)
    out<-data.frame(
      pway=names(hmark),
      count=lapply(hmark, function(h) sum(h %in% thresh.t$ensg)) %>% unlist(),
      gsea.Z=with(t,valtest_genesetlist(ensg, rank(score), hmark) %>% unlist()),
      method=unique(t$method),
      disease=unique(t$disease)
    )
    P<-abs(out$gsea.Z) %>% pnorm(., lower.tail=FALSE)*2
    out$padj<-p.adjust(abs(out$gsea.Z) %>% pnorm(., lower.tail=FALSE)*2)
    out$p<-P
    out$logp<--log(P)
    out
  }))
  rownames(re)<-NULL
  re
}))
subset(non.filt,padj<0.01 & method=="cogs" & gsea.Z>0)




## only show pathways that are significant in at least one method
ggplot(subset(non.filt,padj<0.01 & gsea.Z>0) ,aes(x=count,y=-log(padj),colour=disease)) + 
  geom_point() + facet_grid(.~method) + theme_bw() + xlab("Genes in set with gene score > 0.5") + ylab("-log(p.adj) enrichment") + ggtitle("KEGG")


## for a given pathway extract all the genes from it that are above threshold
pathway<-'KEGG_P53_SIGNALING_PATHWAY'
do.call("rbind",lapply(dt,function(d){
  hmark<-extract.genesets(paste0("^",pathway),unique(tmp$ensg),gmt) %>% unlist()
  subset(d,biotype=="protein_coding" & method=="cogs" & ensg %in% hmark & score>=thresh)
}))


## next we want to filter so we only compare genes that have scores across platforms

filt<-lapply(dt,function(d){
  disease<-unique(d$disease)
  tmp<-melt(d,id.vars=c('ensg','method'),measure.vars = 'score')
  tmp<-dcast(tmp,ensg~method)
  tmp<-tmp[rowSums(is.na(as.matrix(tmp[,2:ncol(tmp)])))==0,]
  tmp$disease<-disease
  tmp
})

names(filt)<-sub("\\.tab","",basename(fls))

r<-lapply(filt,gsea)

pdf(file=file.path(PLOT_DIR,'hallmark.pdf'))
lapply(seq_along(r),function(i){
  x<-r[[i]]
  P <- abs(x) %>% pnorm(., lower.tail=FALSE)*2 %>% p.adjust()
  thr <- min(abs(x)[ P<0.01 ])
  use <- apply(x>=thr,1,any)
  x <- x[rownames(x) %in% names(use)[use],]
  x<-melt(x)
  names(x)<-c('pathway','method','zscore')
  x$pathway<-gsub("HALLMARK\\_","",x$pathway)
  ggplot(x,aes(x=pathway,y=zscore,fill=method)) + geom_bar(stat = "identity", position="dodge") + geom_abline(intercept=thr,slope=0) + coord_flip() + ggtitle(names(r)[i]) + theme_bw()
})
dev.off()
