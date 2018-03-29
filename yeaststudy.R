## script to simulate and analyse go enrichment from yeast data in ensembl

## what is the distribution of genelength depending on GO terms

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

## Load data from ensembl
require(biomaRt)

datachoice=1

## potentially two dataset within biormart
if(datachoice==1)
{
    mart = useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
    ensembl = useDataset("scerevisiae_gene_ensembl", mart = mart)
}
else ##or
{
    mart = useMart("fungi_mart",host="fungi.ensembl.org")
    ensembl = useDataset("scerevisiae_eg_gene", mart = mart)
}
                      
## get GO terms
#Q1 = getBM(attributes =
#              c("go_id","ensembl_gene_id"),#
#             mart=ensembl)

## get gene length
#Q2 = getBM(attributes=c("transcript_length","ensembl_gene_id"),
#          mart=ensembl)

Q = getBM(attributes=c("transcript_length","ensembl_gene_id","go_id","transcript_length"), mart=ensembl)

names(Q)[4]="gene_length"

## merge both
#Q = merge(Q1,Q2)

## save this datset to not have to fetch it again
save(Q,file="go_yeast_150318_fungimart.txt")


load(file="go_yeast_150318_fungimart.txt")

#install.packages("tidyverse")
#install.packages("gridExtra")

library(tidyverse)
library(dplyr)
library(gridExtra)


######## DATA EXPLORATION ##############


#filter by the number of genes per GO
#Q3=list('go_id'=names(which(table(Q[,2])>10)))


## only keep GO with more than 100 terms
Q3=list('go_id'=names(which(table(Q[,3])>10)))
maxnbgenes=max(table(Q[,2]))
Qclean=Q[Q$go_id %in% unlist(Q3),]




Qtable=as_tibble(Qclean) 
Qclean %>% group_by(ensembl_gene_id, go_id) %>% filter(row_number()==1) %>% filter(go_id!="")

length(unique(Qclean$go_id))
## 795
length(unique(Q$go_id))
## 5345
length(unique(Qclean$ensembl_gene_id))
## 6684

## avergae nb of go per gene
summary(as.vector(unlist(table(Qclean$ensembl_gene_id))))
## one gene with far too many go  ENSG00000007372 with 7209 go
hist(log(sort(as.vector(unlist(table(Qclean$ensembl_gene_id))),decreasing=T))) #remove it, actually show it is NOT an outlier
## most gene with 1 go, few genes with ~2000 go

hist(log(sort(as.vector(unlist(table(Qclean$ensembl_gene_id))),decreasing=T)))

## link gene length and nb of go ?


## work on transcript length
minny=tapply(Qclean$transcript_length,Qclean$go_id,mean)
stdy=tapply(Qclean$transcript_length,Qclean$go_id,sd)
miny=tapply(Qclean$transcript_length,Qclean$go_id,min)
maxy=tapply(Qclean$transcript_length,Qclean$go_id,max)
q95=tapply(Qclean$transcript_length,Qclean$go_id,quantile,0.95)
q05=tapply(Qclean$transcript_length,Qclean$go_id,quantile,0.05)


gotable=data.frame(tapply(Qclean$transcript_length,Qclean$go_id,mean))
gotable$nbgenes=tapply(Qclean$transcript_length,Qclean$go_id,length)
colnames(gotable)[1]="meangl"
colnames(gotable)[2]="nbgenes"
hist(log(gotable$nbgenes))
## poisson log distributed
gotable$rank=rank(gotable$meangl)
    ##sort(gotable$meangl,index.return=T)$ix
plot(gotable$rank,gotable$meangl)
gotable$go_id=rownames(gotable)


## remove the outliers with too many genes
##gotable=gotable[-which(log(gotable$nbgenes)>8),]
#gotable=gotable[-which(gotable$go_id==" "),]

low_gl=gotable[gotable$meangl %in% sort(gotable$meangl)[1:20],]

top_gl=gotable[gotable$meangl %in% sort(gotable$meangl,decreasing=T)[1:20],]


Qtop=Qclean[Qclean$go_id %in% rownames(top_gl),]
Qlow=Qclean[Qclean$go_id %in% rownames(low_gl),]

library(ggplot2)

Qtot<-merge(Qclean,gotable,by="go_id")

Qtot %>% filter(nbgenes>20) %>% dim
Qtot %>% dim

ggplot(Qtot,aes(x=log(transcript_length+0.1),col=as.factor(rank),lwd=rank*1.0/1000))+geom_density()+guides(col=F)+  scale_color_manual(values=rainbow(4000)[Qtot$rank])#+geom_density(data=Qlow,mapping=aes(x=log(transcript_length+0.1),col=go_id))

ggplot(Qtop,aes(x=log(transcript_length+0.1),col=go_id))+geom_density()+guides(col=F)
ggplot(Qlow,aes(x=log(transcript_length+0.1),col=go_id))+geom_density()+guides(col=F)


ggplot(Qtop,aes(x=log(transcript_length+0.1),col=rank))+geom_density()+guides(col=F)+geom_density(data=Qlow,mapping=aes(x=log(transcript_length+0.1),col=go_id))
    
ggplot(Qlow,aes(x=log(transcript_length+0.1),col=go_id))+geom_density()+guides(col=F)

## is there a relationship between meangl and the size of go?
pdf("nbgenes_vs_meangl_pergo_emp_yeast.pdf",height=5,width=5)
ggplot(gotable,aes(x=log(meangl),y=log(nbgenes)))+geom_point()+geom_smooth(col="red")+ylab("nb of genes per GO category")+xlab("mean( log(gene length)) per GO")+theme_minimal()
dev.off()

pdf("histgo_vs_meangl_pergo_emp_yeast.pdf",height=5,width=5)
ggplot(gotable,aes(x=log(meangl)))+geom_histogram(fill="red")+ylab("nb of GO category")+xlab("mean( log(gene length)) per GO")+  theme_minimal() 
dev.off()

cor(gotable$meangl,gotable$nbgenes,method="pearson")
##-0.026 pearson  so not really
cor(gotable$meangl,gotable$nbgenes,method="spearman")
## 0.066 spearman

mygo=names(sort(minny))
ngo=length(mygo)

#png('go,png')
plot(tan(sqrt(sort(minny))),col='red',type='l',lwd=2,main='dist of gene length',xlab='GO',ylab='gene length')#ylim=c(0,8000),
points(q95[mygo],col='purple',pch='-',lwd=2)
points(q05[mygo],col='purple',pch='-',lwd=2)
##points(stdy[mygo],col='green',pch='-')
#points(miny[mygo],col='blue',pch='-')
#points(maxy[mygo],col='purple',pch='-')
#dev.off()

rsq=numeric(ngo)
rsq2=numeric(ngo)
cof=numeric(ngo)
Tt=hist(Qclean$transcript_length,plot=F)
plot(Tt$mids,Tt$density,type='l',lwd=6,log='xy',xlab="transcript length",ylab="frequency per GO term",main="freq of gene length in go terms (log scale)")

load("dotest.R")

for(i in 1:ngo)
{
    Tt=hist(Qclean[which(Qclean$go_id==mygo[i]),]$transcript_length,plot=F)
    ##    points(T$mids,T$density,col=rainbow(100)[i%%100],pch="-",cex=5)
    ##    lines(T$mids,T$density,col=rainbow(100)[i%%100],lwd=5)
    y=lm((Tt$density)~(Tt$mids))
    rsq[i]=summary(y)['adj.r.squared']
    cof[i]=y$coefficients[2]
    y=lm(log(T$density+0.000000001)~log(T$mids))
    rsq2[i]=summary(y)['adj.r.squared']    
    ##    myp=table((Q[which(Q$go_id==mygo[i]),3]%%500)*500)
    ##    myp=myp/sum(myp)
    ##    points(myp,col=rainbow(50)[i%%50],lty=2,pch='-',lwd=5)
    ##    lines(names(myp),myp,col=rainbow(50)[i%%50],lty=2,lwd=1)
}
lines(Tt$mids,Tt$density,type='l',lwd=3)


sorted <- sort(unlist(rsq),index.return=TRUE)$ix
plot(unlist(rsq)[sorted],col=rainbow(50),ylim=c(0,1))
points(unlist(rsq2)[sorted],col=rainbow(50))

sorted <- sort(unlist(rsq2),index.return=TRUE)$ix
plot(unlist(rsq2)[sorted],col=rainbow(300),ylim=c(0,1))
points(unlist(rsq)[sorted],col=rainbow(300))

plot(unlist(rsq2)[sorted],col=abs(unlist(rsq)[sorted]*100),ylim=c(0,1))


############# DO SIMULATION ###########

simugenes=Qtot[,c(1,3,2)]
names(simugenes)[3]="genelength"
names(simugenes)[2]="gene_id"
simugo=gotable[,c(1,2,4)]

source("dosimu.R")

for(i in 1:1000)
{
    doempsimunew(1,simugenes,simugo,0,10000,foldername="yeastsimu15mar/")
}


############## DO TESTING #############

source("dotest.R")

foldername="yeastsimu15mar/"
tot=matrix(0,nrow=1000,ncol=17)
tot2=matrix(0,nrow=1000,ncol=17)
tot3=matrix(0,nrow=1000,ncol=17)
##foldername="yeastsimu/"
##tot2=matrix(0,nrow=50000,ncol=17)
i=0
for( filename in list.files(path=foldername))
{
    print(paste("processing",paste(foldername,filename,sep="")))
    print(i)
    i=i+1
    t=dotestnew(paste(foldername,filename,sep=""))
    tot[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    t=dotestJensen(paste(foldername,filename,sep=""))
    tot2[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
#    t=dotestgoseq(paste(foldername,filename,sep=""))
#    tot3[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))    
    ##    t=dotestPerm(paste(foldername,filename,sep=""))
    ##    tot3[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    ##   save(tot,tot2,tot3,file="emptot.Rdata")
#    plotTFNP(tot)
}
save(tot,file="emptot.R")


plotTFNP <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=TP1/tot,y=TP2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE POSITIVES")+ylim(c(0,0.5))
    g2=ggplot(mydf,aes(x=FP1/tot,y=FP2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE POSITIVES")+ylim(c(0,0.1))
    g3=ggplot(mydf,aes(x=TN1/tot,y=TN2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=FN1/tot,y=FN2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" new method")+ylab("goseq")+ggtitle("FALSE NEGATIVES")+xlim(c(0,0.5))
    grid.arrange(g1,g2,g3,g4,ncol=4)
}

plotTFNP2 <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=TP1/tot,y=FP1/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" TP")+ylab("FP")+ggtitle("TRUE POSITIVES")+ylim(c(0,0.5))
    g2=ggplot(mydf,aes(x=TP2/tot,y=FP2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" TP")+ylab("FP")+ggtitle("FALSE POSITIVES")+ylim(c(0,0.1))
    g3=ggplot(mydf,aes(x=TN1/tot,y=FN1/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" TN")+ylab("FN")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=TN2/tot,y=FN2/tot,col=gosel))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("TN")+ylab("FN")+ggtitle("FALSE NEGATIVES")+xlim(c(0,0.5))
    grid.arrange(g1,g2,g3,g4,ncol=4)
}


plotTFNP3 <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=TP1/tot+FP1/tot,y=gosel,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FP+TP")+ylab("gosel")+ggtitle("Test 1")
    g2=ggplot(mydf,aes(x=TP2/tot+FP2/tot,y=gosel,col=TP2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FP+TP")+ylab("gosel")+ggtitle("Test 2")
    g3=ggplot(mydf,aes(x=TN1/tot+FN1/tot,y=gosel,col=TN1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FN+TN")+ylab("gosel")+ggtitle("Test 1")
    g4=ggplot(mydf,aes(x=TN2/tot+FN2/tot,y=gosel,col=TN2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("FN+TN")+ylab("gosel")+ggtitle("Test 2")
    grid.arrange(g1,g2,g3,g4,ncol=4)
}


plotTFNP4 <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=FP1/tot,y=gosel,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FP")+ylab("gosel")+ggtitle("Test 1")
    g2=ggplot(mydf,aes(x=FP2/tot,y=gosel,col=TP2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FP")+ylab("gosel")+ggtitle("Test 2")
    g3=ggplot(mydf,aes(y=gosel,x=FN1/tot,col=TN1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab(" FN")+ylab("gosel")+ggtitle("Test 1")
    g4=ggplot(mydf,aes(y=gosel,x=FN2/tot,col=TN2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+xlab("FN")+ylab("gosel")+ggtitle("Test 2")
    grid.arrange(g1,g2,g3,g4,ncol=4)
}


plotTFN <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=FN1/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g2=ggplot(mydf,aes(x=TN1/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g3=ggplot(mydf,aes(x=FN2/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=TN2/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    g5=ggplot(mydf,aes(x=FP1/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g6=ggplot(mydf,aes(x=TP1/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g7=ggplot(mydf,aes(x=FP2/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g8=ggplot(mydf,aes(x=TP2/tot,y=gosel,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4)
}

plotTFNparam2 <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=FN1/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g2=ggplot(mydf,aes(x=TN1/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g3=ggplot(mydf,aes(x=FN2/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=TN2/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    g5=ggplot(mydf,aes(x=FP1/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g6=ggplot(mydf,aes(x=TP1/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g7=ggplot(mydf,aes(x=FP2/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g8=ggplot(mydf,aes(x=TP2/tot,y=param2,col=param2))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4)
}

plotTFNparam3 <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=FN1/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g2=ggplot(mydf,aes(x=TN1/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g3=ggplot(mydf,aes(x=FN2/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=TN2/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    g5=ggplot(mydf,aes(x=FP1/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g6=ggplot(mydf,aes(x=TP1/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g7=ggplot(mydf,aes(x=FP2/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g8=ggplot(mydf,aes(x=TP2/tot,y=param3,col=param3))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE NEGATIVES")
    grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,ncol=4)
}

plotTop <- function(tot){
    mydf=data.frame(tot)
    names(mydf)=c("acc1","acc2","gosel","param2","param3","cor1","cor2","FN1","TN1","TP1","FP1","FN2","TN2","TP2","FP2")
    for(i in 1:17)
        mydf[,i]=as.numeric(mydf[,i])
#    for(i in 6:17)
 #       mydf[,i]=as.numeric( levels(mydf[,i]) )[mydf[,i]]
    mydf$tot=mydf$FN1+mydf$TP1+mydf$FP1+mydf$TN1
    mydf$tot2=mydf$FN2+mydf$TP2+mydf$FP2+mydf$TN2
#    mydf$tot3=mydf$FN3+mydf$TP3+mydf$FP3+mydf$TN3
    g1=ggplot(mydf,aes(x=cor1-gosel,y=param3,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g2=ggplot(mydf,aes(x=cor1-gosel,y=gosel,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g3=ggplot(mydf,aes(x=cor1-gosel,y=param2,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    g4=ggplot(mydf,aes(x=cor2-gosel,y=param3,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE POSITIVES")
    g5=ggplot(mydf,aes(x=cor2-gosel,y=gosel,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("FALSE POSITIVES")
    g6=ggplot(mydf,aes(x=cor2-gosel,y=param2,col=TP1))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+ggtitle("TRUE NEGATIVES")
    grid.arrange(g1,g2,g3,g4,g5,g6,ncol=3)
}




plotTFNP(tot)
plotTFNP2(tot)
plotTFNP3(tot)
plotTFNP4(tot)
plotTFN(tot)
plotTFNparam2(tot)
plotTFNparam3(tot)
plotTop(tot)

pdf("plot221217_yeast1.pdf")
plotTFNP(tot2)
dev.off()

pdf("plot221217_yeast2.pdf")
plotTFNP2(tot2)
dev.off()

pdf("plot221217_yeast3.pdf")
plotTFNP3(tot2)
dev.off()

pdf("plot221217_yeast4.pdf")
plotTFNP4(tot2)
dev.off()

pdf("plot221217_yeast5.pdf")
plotTFN(tot2)
dev.off()

pdf("plot221217_yeast6.pdf")
plotTFNparam2(tot2)
dev.off()

pdf("plot221217_yeast7.pdf")
plotTFNparam3(tot2)
dev.off()

pdf("plot221217_yeast8.pdf")
plotTop(tot2)
dev.off()
