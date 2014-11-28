library(ggplot2)
library(qvalue)
library(DESeq2)

path<-"/home/lpantano/iRNAseq/znf257/znf257.DE/"
setwd(path)
id<-read.table("ensemb2id",sep="\t",header=T)
id<-id[id[,4]!="",]
gensum<-read.table("~/crickhome/iRNAseq/annotation/gene.symmary.sep",header=T,sep="\t")

 #### read files #############
list.s<-vector()
table<-data.frame()
for (s in dir(path,paste(sep="","genes$"))){
  #print(s)
  name<-strsplit(s,"[::.::]")[[1]][1]
  list.s<-c(list.s,name)
  t<-read.table(paste(sep="",path,s))
  names(t)<-c("others",name)
  if (length(list.s)==1){table<-t}else{
    table<-cbind(table,t[,2])
    names(table)[ncol(table)]<-name
  }
}

total<-apply(table[,2:9],1,sum)
tab<-as.matrix(table[total>=30,2:9])
row.names(tab)<-table[total>=30,1]


 #### DESeq2 #############
design<-data.frame(row.names=colnames(tab),pop=c("CHB","CHB","JPT","JPT","JPT","JPT","JPT","JPT"),sex=c("Male","Male","Female","Female","Female","Female","Female","Male"),condition=c("STD","HET","STD","STD","STD","HET","HET","HET"))
dse <- DESeqDataSetFromMatrix(countData = tab, colData = design,design = ~ sex+pop+condition)
dse <- DESeq(dse,fitType="parametric")
res <- results(dse)
save(dse,file=paste(sep="",path,"dse"))
resultsNames(dse)

selected<-as.data.frame(res[res[,5]<=0.2 & !is.na(res[,4]),])
write.table(selected,paste(sep="",path,"DE_output.txt"),sep="\t",row.names=T)
res["ENSG00000197134",]
head(res[order(res$pvalue),])
tab.n<-counts(dse,normalized=TRUE)

as.data.frame.DataFrame <- selectMethod("as.data.frame", "DataFrame")
selected<-as.data.frame(res[res[,4]<=0.01 & !is.na(res[,4]),])
selected$gene<-row.names(selected)
selected<-merge(selected[,c(6,1:5)],gensum,by=1,all.x=T)
write.table(selected[order(selected$padj),],paste(sep="",path,"top0-01.txt"),sep="\t",row.names=F)

selected<-as.data.frame(res[res[,5]<=0.1 & !is.na(res[,4]),])
write.table(row.names(selected),paste(sep="",path,"topFDR0-1.txt"),sep="\t",col.names=F,quote=F,row.names=F)
write.table(row.names(res),paste(sep="",path,"topFDR0-1.node.txt"),sep="\t",col.names=F,quote=F,row.names=F)

selected<-as.data.frame(res[res[,5]<=0.35 & !is.na(res[,4]),])
selected.f<-merge(selected[,c(6,2)],id[,c(1,4)],by=1)
write.table(selected.f[order(selected[,2],decreasing=T),c(3,2)],"selected_gsea_35.rnk",col.names=F,sep="\t",quote=F,row.names=F)
 
 #### qvalue analysis #############
qobj <- qvalue(res$pvalue[!is.na(res$pvalue)],fdr.level=0.1)
qwrite(qobj,filename=paste(sep="",path,"qvalue.slides"))
t<-qsummary(qobj)
pi0<-t$pi0
pdf(paste(sep="",path,"qvalue.pdf"))
qplot.log<-try(qplot(qobj),silent=TRUE)
dev.off()


 #### plot selected genes #############
genes<-selected[selected[,6]<=0.01,1]
for (i in genes){
  exp<-log2(tab.n[i,]+1)
  pv<-format(res[i,4],scientific=TRUE,digits=3)
  fdr<-format(res[i,5],scientific=TRUE,digits=3)
  d<-data.frame(exp=exp,genotype=design$condition,sex=design$sex)
  p<-ggplot(d, aes(x=factor(genotype), y=exp,fill=factor(genotype))) +
    geom_boxplot(outlier.size = 0) + 
    geom_jitter(position=position_jitter(width=0.4),aes(factor(genotype), exp,colour=factor(sex))) +
    scale_color_manual("",values=c("blue2","green2","blue3","green3"))+
    scale_fill_manual("",values=c("LightBlue2","DarkSeaGreen3","LightBlue3","DarkSeaGreen4"))+
    theme_bw(base_size = 12, base_family = "") +
    labs(list(title=paste(selected[selected$gene==i,7],"\npv:",pv,"\nfdr:",fdr),y="log2(norm_counts)",x=""))
  pdf(paste(sep="",path,i,".pdf"))
  print(p)
  dev.off()
}

 #### plot exon levels #############
g<-read.table("znf257.junction.tab",row.names=1,header=T)
l<-data.frame(row.names=as.character(names(sizeFactors(dse))),lib=sizeFactors(dse))
d<-colData(dse)
g.n<-data.frame(exp=t(g),lib=l[as.character(names(g)),],sex=d[as.character(names(g)),2],genotype=d[as.character(names(g)),3])
g.n$norm<-g.n[,1]*(1.2/g.n$lib)
p<-ggplot(g.n, aes(x=factor(genotype), y=e1e2ZNF257,fill=factor(genotype))) +
  geom_boxplot(outlier.size = 0) + 
  geom_jitter(position=position_jitter(width=0.4),aes(factor(genotype), e1e2ZNF257,colour=factor(sex))) +
  scale_color_manual("",values=c("blue2","green2","blue3","green3"))+
  scale_fill_manual("",values=c("LightBlue2","DarkSeaGreen3","LightBlue3","DarkSeaGreen4"))+
  theme_bw(base_size = 12, base_family = "") +
  labs(list(title="First ee junstion of ZNF257",y="log2(norm_counts)",x=""))
pdf(paste(sep="","exon1exon2ZNF257.pdf"))
print(p)
dev.off()


 #### GSEA parsed #############
gsea.sum<-read.table("/home/lpantano/iRNAseq/annotation/msigdb.v3.1.symbols.genes",sep="\t",row.names=1)

resrank<-read.table("/home/lpantano/znf257/znf257.DE/gsea.top/msigdb.GseaPreranked.1381400117458/edb/results.lastrank.tab",sep="\t",row.names=1)
resr<-read.table("/home/lpantano/znf257/znf257.DE/gsea.top/msigdb.GseaPreranked.1381400117458/edb/results.parsed.tab",header=F,row.names=2)
top<-resr[resr$V8<=0.01 & resr$V4>0,]
#dse should be loaded
res <- results(dse)
rank.last.sign<-max(which.max(res$pvalue[res$pvalue<=0.07]))
list.here<-vector()
for (row in row.names(top)){#GO TRHOUGH EACH PATHWAY
        print(row)
        gens<-unlist(strsplit(as.character(gsea.sum[row,]),","))
      if (length(gens[!is.na(gens)])>0){
        conver<-merge(id[,2:1],gens,by=1)
        conver<-conver[!duplicated(conver[,2]),]
        row.names(conver)<-conver[,2]
        table<-conver[,2:1]#GENES AND ESEMBLE IDS
        npval<-vector()
        selected<-data.frame()
        
        #SELECT ONLY THE GENES IN THE PATHWAY     
        selected<-res[as.character(conver[,2]),]
        selected<-selected[!is.na(selected[,1]),]
        selected<-as.data.frame(selected)
        npval<-sum(selected$pval<=0.05,na.rm=T)
        selected$name<-row.names(selected)
        #names(selected)[4:5]<-paste(labs[idx],c("pval","FDR"),sep=".")
        #table<-merge(table,selected[,c(6,4,5)],by=1)
      
        #GET QUANTILES 0.5 0.10 0.25 0.35 OF PVALUE AND FDR
        #quan.pval<-quantile(selected[,4],c(.10,.25,.35),na.rm =TRUE)
        #quan.fdr<-quantile(selected[,5],c(.10,.25,.35),na.rm =TRUE)
        
        true.rank<-resrank[row,1]<=rank.last.sign
        #GET ONLY PATHWAYS WHERE LAST PVALUE GENE<=LAST RANKED GENE
        #AND Q25<=0.05
        #true.quan<-quan.pval[2]<=0.05
        #print(resrank[row,1])
        #print(true.rank)
        #if (true.quan & true.rank){
        if (true.rank){
          list.here<-c(list.here,row)
          
        }
      }
}
