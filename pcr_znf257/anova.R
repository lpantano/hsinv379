data<-read.csv("ExpressionDataDEG.csv",dec=",",header=T,sep=";")

tab<-data.frame(row.names=names(data)[3:13])
for (g in names(data)[3:13]){
    exp1<-data[,g]
    g1<-rep(1,5)  
    g2<-rep(0,5)  
    y1<-c(g1,g2)
    lm1<-lm(y1~exp1)
    g1<-rep(1,4)
    y2<-c(0,g1,g2)
    lm2<-lm(y2~exp1)
    #anova(lm1,lm2)
    tab[g,1]<-AIC(lm1)
    tab[g,2]<-AIC(lm2)
}

tab$dif<-tab[,1]-tab[,2]
names(tab)<-c("std","inv","dif")
write.table(tab,"best.model",sep="\t",col.names=T,quote=F)
