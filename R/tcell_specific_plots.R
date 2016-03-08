library(ggplot2)

options(stringsAsFactors=FALSE)

data.dir<-'/home/oliver/JAVIERRE_GWAS/out/blockshifter_two_tail/'
gwas.meta<-read.csv('/home/oliver/JAVIERRE_GWAS/support/gwas_manifest.csv')

###barplot first

f<-list.files(path=data.dir,pattern='*.txt',full.names = TRUE)

all.bs<-do.call('rbind',lapply(f,read.delim))

### get CD4 act vs non act

bp.d<-subset(all.bs,test=='Total_CD4_NonActivated' & control=='Total_CD4_Activated')
bp.d$bp.val<-with(bp.d,-log(p.emp.twotail) * -sign(z))
bp.d<-bp.d[order(bp.d$bp.val,decreasing = TRUE),]

## add category

gw.lu<-split(gwas.meta$category,gwas.meta$label)
bp.d$category<-do.call('c',gw.lu[bp.d$gwas])
bp.d$gwas<-factor(bp.d$gwas,levels = bp.d$gwas)

ggplot(bp.d,aes(x=gwas,y=bp.val,fill=category)) + geom_bar(stat="identity") + theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("Non Act. CD4+ <-log(p)> Act. CD4+")

## SCATTER PLOT
## Due to some issues with the way I submit jobs some of the comparisons
## are back to front the code below fixes this by reversing the sign of the z score
bp.control.d<-subset(all.bs,test=='Control' | control=='Control')
idx<-which(bp.control.d$test=='Control')
bp.control.d[idx,]$test<-bp.control.d[idx,]$control
## for those that are zero set to 1/perm
zidx<-which(bp.control.d$p.emp.twotail==0)
bp.control.d[zidx,]$p.emp.twotail<-1/bp.control.d[zidx,]$perm
bp.control.d[idx,]$z<-bp.control.d[idx,]$z * -1
bp.control.d$score<-with(bp.control.d,-log(p.emp.twotail) * sign(z))

library(reshape2)

dm<-melt(bp.control.d,id=c('gwas','test'),measure.vars ='score')
p.d<-dcast(dm,gwas ~ test)
p.d$category<-do.call('c',gw.lu[p.d$gwas])
p.d$gwas<-factor(p.d$gwas,levels = p.d$gwas)


##ggplot(p.d,aes(x=Total_CD4_Activated,y=Total_CD4_NonActivated,color=category,pch=category,label=gwas)) + geom_point(size=5) + geom_text(size=5,hjust=-0.2, vjust=0, angle=0, show_guide=FALSE)


ggplot(data=p.d,aes(x=Total_CD4_Activated,y=Total_CD4_NonActivated ,pch=category,colour=category,label=gwas)) +
  geom_point(size=3) + theme_bw() +  xlab("< Ery + Mega vs CD4+ Act >") +
  ylab("< Ery + Mega vs CD4+ Non Act >") +
  geom_text(size=5,hjust=-0.2, vjust=0, angle=0, show_guide=FALSE) +
  scale_colour_manual(name = NULL,values = c(Autoimmune="red",Blood="black",Metabolic="darkgreen",Other="blue",T2D="yellow")) + scale_shape_discrete(name=NULL) + coord_flip() + ylim(c(-10,10)) +
  geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) + guides(size=FALSE,text=FALSE) +
  theme(legend.position=c(0.85, 0.3), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))
