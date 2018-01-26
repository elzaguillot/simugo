
## script to analyses bunches fo simulation created in a specific folder

## analysis of single simulations

library(dplyr)
library(ggplot2)
require(gridExtra)
require(goseq)

foldername="simusept/"

tot2=matrix(0,nrow=1000,ncol=23)
#colnames(saveOutput)=c("acc1","acc2","acc3","cor1","cor2","cor3","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2","FN3","TN3","TP3","FP3")
i=1
for( l in list.files(foldername)){
    filename=paste(foldername,l,sep="")
    load(paste("simusept/"cool,l,sep=""))
    simugenes$go_id=gsub(":","",simugenes$go_id)
    simugenes$gene_id=gsub(":","",simugenes$gene_id)
    simugo$go_id=gsub(":","",simugo$go_id)
    rownames(simugo)=gsub(":","",simugo$go_id)
    lengthbygo=tapply(simugenes$genelength,simugenes$go_id,mean,na.rm=T)
    lengthbygo2=tapply(simugenes$genelength,simugenes$go_id,median,na.rm=T)
    lengthbygosd=tapply(simugenes$genelength,simugenes$go_id,sd,na.rm=T)
    lbgo=as.data.frame(lengthbygo)
    names(lbgo)[1]="genelength"
    lbgo$go_id=rownames(lbgo)
    lbgo$medgenelength=lengthbygo2
    lbgo$sdgenelength=lengthbygosd    
    result=dotestnewdf(simugo,simugenes)
    gofinal=result[[1]]
    names(gofinal)[1]="go_id"
    param=result[[3]]
    mysimu <- gofinal %>% inner_join(simugo,by="go_id") %>% inner_join(lbgo,by="go_id")
    summary(mysimu)
#    g1 <- ggplot(mysimu,aes(x=genelength,y=pval1,col=selected))+geom_point()+geom_smooth()+ggtitle(paste(param,sep=" , ",collapse=" - "))
#    g2 <- ggplot(mysimu,aes(x=genelength,y=pval2,col=selected))+geom_point()+geom_smooth()+ggtitle(paste(i))
#    g3 <- ggplot(mysimu,aes(x=genelength,y=log(nbgenes),col=selected))+geom_point()+geom_smooth()+ggtitle(paste(i))
#    grid.arrange(g1,g2,g3,ncol=3)
#    invisible(readline(prompt="Press [enter] to continue"))
    ## newdata <- mutate(data,EnsemblGeneID=transformID(EnsemblGeneID))  %>% filter(!is.null(EnsemblGeneID ))%>% distinct(EnsemblGeneID,.keep_all=T)
    assayed.genes <- unique(simugenes$gene_id)
    de.genes <- simugenes$gene_id[(simugenes$selected==1)]
    gene.vector=as.integer(assayed.genes%in%de.genes[[1]])
    pwf1 <- nullp(gene.vector,bias.data=simugenes$genelength) # new_lengthData)
    rownames(pwf1)=simugenes$gene_id
    go_map=data.frame("genes"=simugenes$gene_id,"category"=simugenes$go_id)
    go_map=as.list(tapply(simugenes$go_id,simugenes$gene_id,as.vector))
    goseq_go=goseq(pwf1,gene2cat=go_map)
    names(goseq_go)[1]="go_id"
    mysimu <- gofinal %>% inner_join(simugo,by="go_id")%>% inner_join(lbgo,by="go_id") %>% inner_join(goseq_go,by="go_id")
#    over_represented_pvalue
#    g1 <- ggplot(mysimu,aes(x=genelength,y=pval1,col=selected))+geom_point()+geom_smooth()+ggtitle("New method")+ylab("corrected pvalue")+xlab("log(gene length")+theme_minimal()
#    g2 <- ggplot(mysimu,aes(x=genelength,y=pval2,col=selected))+geom_point()+geom_smooth()+ggtitle("no correction")+ylab("corrected pvalue")+xlab("log(gene length")+theme_minimal()
#    g3 <- ggplot(mysimu,aes(x=genelength,y=over_represented_pvalue,col=selected))+geom_point()+geom_smooth()+ggtitle("goseq")+ylab("corrected pvalue")+xlab("log(gene length)")+theme_minimal()
    ##  i=i+1
 #   if(i<20){ #save juste a few
 #   pdf(paste("exple_pval_goseqetal",i,".pdf",sep=""),width=12,height=5)
#    grid.arrange(g1,g2,g3,ncol=3)
#    dev.off()}
    ##    invisible(readline(prompt="Press [enter] to continue"))
    tot2[i,]=docompare2(mysimu,param)
    i=i+1
}

mydf=data.frame(tot)
names(mydf)=c("acc1","acc2","acc3","cor1","cor2","cor3","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2","FN3","TN3","TP3","FP3","gosel","param2","param3","param4")
for(i in 1:23)
    mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3

#library(goseq)

ggplot(mydf,aes(acc1/tot,acc2/tot))+geom_point()+geom_abline(slope=1,intercept=0,col="red")


ggplot(mydf,aes(acc1/tot))+geom_histogram()

pdf("simu_mevsnocorrection.pdf",width=12,height=5)
ggplot(mydf,aes(x=acc1/tot,y=acc2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("percent of accurate calls with new method")+ylab("percent of accurate calls without correction")
dev.off()

pdf("simu_mevsgoseq.pdf",width=12,height=5)
ggplot(mydf,aes(x=acc1/tot,y=acc3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("percent of accurate calls with new method")+ylab("percent of accurate calls in goseq")
dev.off()

pdf("simu_accuratecalls.pdf",width=12,height=5)
g1=ggplot(mydf,aes(x=acc1/tot,y=acc3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("proportion of accurate calls with new method")+ylab("proportion of accurate calls in goseq")+ggtitle("comparison with GOSEQ")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+labs(col="proportion of go selected")+xlim(c(0.75,1))+ylim(c(0.75,1))+theme_minimal()
g2=ggplot(mydf,aes(x=acc1/tot,y=acc2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("proportion of accurate calls with new method")+ylab("porportion of accurate calls without correction")+ggtitle("comparison with no correction")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+labs(col="proportion of go selected")+xlim(c(0.75,1))+ylim(c(0.75,1))+theme_minimal()
grid.arrange(g1,g2,ncol=2)
dev.off()

pdf("simu_accuratecalls_vsgosel.pdf",width=12,height=5)
g1=ggplot(mydf,aes(x=acc3/tot,y=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("proportion of accurate calls with new method")+ylab("proportion of accurate calls in goseq")+ggtitle("comparison with GOSEQ")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+labs(col="proportion of go selected")+xlim(c(0.75,1))+ylim(c(0,0.1))+theme_minimal()
g2=ggplot(mydf,aes(x=acc1/tot,y=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("proportion of accurate calls with new method")+ylab("porportion of accurate calls without correction")+ggtitle("comparison with no correction")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+labs(col="proportion of go selected")+xlim(c(0.75,1))+ylim(c(0,0.1))+theme_minimal()
grid.arrange(g1,g2,ncol=2)
dev.off()


ggplot(mydf,aes(cor1,cor3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("correlation new method")+ylab("correlation goseq")

ggplot(mydf,aes(cor1,cor2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("correlation new method")+ylab("correlation no bias")


mysimu <- gofinal %>% inner_join(simugo,by="go_id")%>% inner_join(lbgo,by="go_id") %>% inner_join(lbgo,by="go_id")


#pdf("goseq_vs_new_TPFN.pdf")
pdf("goseq_vs_new_TPFN.pdf",width=12,height=5)
g1=ggplot(mydf,aes(x=TP1/tot,y=TP3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE POSITIVES")+ylim(c(0,0.5))
g2=ggplot(mydf,aes(x=FP1/tot,y=FP3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE POSITIVES")+ylim(c(0,0.1))
g3=ggplot(mydf,aes(x=TN1/tot,y=TN3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE NEGATIVES")
g4=ggplot(mydf,aes(x=FN1/tot,y=FN3/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE NEGATIVES")+xlim(c(0,0.5))
grid.arrange(g1,g2,g3,g4,ncol=4)
dev.off()

pdf("old_vs_new_TPFN.pdf",width=12,height=5)
g1=ggplot(mydf,aes(x=TP1/tot,y=TP2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE POSITIVES")+ylim(c(0,0.5))
g2=ggplot(mydf,aes(x=FP1/tot,y=FP2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE POSITIVES")+ylim(c(0,0.1))
g3=ggplot(mydf,aes(x=TN1/tot,y=TN2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE NEGATIVES")
g4=ggplot(mydf,aes(x=FN1/tot,y=FN2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE NEGATIVES")+xlim(c(0,0.5))
grid.arrange(g1,g2,g3,g4,ncol=4)
dev.off()

save(mydf,file="30nov_simu.Rdata")

load("30nov_simu.Rdata")
