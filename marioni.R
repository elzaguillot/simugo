## Script that study the empirical dataset from the 2008 Marioni paper

#data=read.table('marionitab2.txt')


require(gridExtra)

require(goseq)
require(ggplot)
require(dplyr)
require(GenomicFeatures)
require(TxDb)

data=read.table('marionitab3.txt',h=T)



head(data)

summary(data)


## convert the gene ID from version 48 to hg19
## use convert history ID tool on ensembl
## this link http://www.ensembl.org/Homo_sapiens/Tools/IDMapper?db=core
## with file moriani_ID48.txt
write.table(data$EnsemblGeneID,file="moriani_ID48.txt",row.names=F,col.names=F,quote=F)

convert_table=read.table("moriani_IDconverted.csv",sep=",",h=T)
newID=convert_table$Matched.ID.s
names(newID)=convert_table$Requested.ID

head(data)

transformID <- function(x)
{
    if(x %in% newID)
        return(newID[x])
    else
        return(NULL)
    }


data$gl=log(data$GeneEnd-data$GeneStart)
plot(data$gl,data$Affy.qvalue)
points(data$gl,data$Illumina.qvalue,col='red')

newdata <- mutate(data,EnsemblGeneID=transformID(EnsemblGeneID))  %>% filter(!is.null(EnsemblGeneID ))%>% distinct(EnsemblGeneID,.keep_all=T)
head(newdata)


cor(data$gl,data$Affy.qvalue,method='spearman')
cor(data$gl,data$Illumina.qvalue,use='pairwise.complete',method='spearman')
cor(data$gl,data$Affy.pvalue)
cor(data$gl,data$Illumina.pvalue,use='pairwise.complete')

library(ggplot2)
ggplot(data,aes(x=Illumina.qvalue,y=Affy.qvalue,col=data$gl))+geom_point()#

ggplot(data,aes(x=log(Illumina.pvalue),y=log(Affy.pvalue),col=data$gl))+geom_point()#



#pdf("Moriani2008_gl.pdf")
plot1 <- ggplot(data,aes(y=log(Illumina.qvalue),x=gl,col=data$gl))+geom_point()+geom_smooth(col="red")+labs(colour="log(gene length)")+xlab("log(gene length)")+ylab("Illumina pvalue")+ggtitle("Illumina")#
plot2 <- ggplot(data,aes(y=log(Affy.qvalue),x=gl,col=data$gl))+geom_point()+geom_smooth(col="red")+ylab("Affymetrix pvalue")+ggtitle("Affymetrix")#
grid.arrange(plot1,plot2,ncol=2)
#dev.off()

ggplot(data,aes(x=Illumina.qvalue,y=Affy.qvalue,col=data$gl))+geom_quantile()#


assayed.genes <- unique(newdata$EnsemblGeneID)
de.genes <- newdata %>% filter(Affy.pvalue<0.05) %>% dplyr::select(EnsemblGeneID) %>% unique
gene.vector=as.integer(assayed.genes%in%de.genes[[1]])
names(gene.vector)=assayed.genes

head(gene.vector)
pwf <- nullp(gene.vector,"hg19","ensGene")
go_map=getgo(names(gene.vector),"hg19","ensGene")
go=goseq(pwf,gene2cat=go_map)
go %>% top_n(n=20,1-over_represented_pvalue)


#cfamiliaris_gene_ensembl <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#txsByGene=transcriptsBy(cfamiliaris_gene_ensembl,"gene")
#lengthData=median(width(txsByGene))
#select_genes <- as.vector(names(lengthData)%in%names(de.genes))
#new_lengthData <- lengthData[select_genes]

## affymetrix
newdata <- mutate(data,EnsemblGeneID=transformID(EnsemblGeneID))  %>% filter(!is.null(EnsemblGeneID ))%>% distinct(EnsemblGeneID,.keep_all=T)
assayed.genes <- unique(newdata$EnsemblGeneID)
de.genes <- newdata %>% filter(Affy.pvalue<0.05) %>% dplyr::select(EnsemblGeneID) %>% unique
gene.vector=as.integer(assayed.genes%in%de.genes[[1]])
names(gene.vector)=assayed.genes
pwf1 <- nullp(gene.vector,bias.data=newdata$gl)#new_lengthData)
go_map=getgo(names(gene.vector),"hg19","ensGene")
affy_go=goseq(pwf1,gene2cat=go_map)
topaffy <- affy_go %>% top_n(n=20,1-over_represented_pvalue)


newdata <- mutate(data,EnsemblGeneID=transformID(EnsemblGeneID))  %>% filter(!is.null(EnsemblGeneID ))%>% distinct(EnsemblGeneID,.keep_all=T)
assayed.genes <- unique(newdata$EnsemblGeneID)
de.genes <- newdata %>% filter(Illumina.pvalue<0.05) %>% dplyr::select(EnsemblGeneID) %>% unique
gene.vector=as.integer(assayed.genes%in%de.genes[[1]])
names(gene.vector)=assayed.genes
pwf2 <- nullp(gene.vector,bias.data=newdata$gl)#new_lengthData)
go_map=getgo(names(gene.vector),"hg19","ensGene")
illumina_go=goseq(pwf2,gene2cat=go_map)
topillu <- illumina_go %>% top_n(n=20,1-over_represented_pvalue)


plot(illumina_go$over_represented_pvalue,affy_go$over_represented_pvalue)
allgo <- illumina_go%>% left_join(affy_go,by="category") 

## compare affymetrix and illumina
## different pvalue for go (quite different !)
## but almost no different call of in out DE and very slight difference in pwf
## why ? how ?
ggplot(allgo,aes(x=under_represented_pvalue.x,y=under_represented_pvalue.y,col=log(numDEInCat.x)))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")
ggplot(allgo,aes(x=numDEInCat.x,y=numDEInCat.y,col=log(numDEInCat.x)))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")


flatpwf=pwf1
flatpwf$pwf=1
flat_go=goseq(flatpwf,gene2cat=go_map)
topflat <- flat_go %>% top_n(n=20,1-over_represented_pvalue)

allgo <- illumina_go%>% left_join(flat_go,by="category") 
p1 <- ggplot(allgo,aes(x=under_represented_pvalue.x,y=under_represented_pvalue.y,col=log(numDEInCat.x)))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")

allgo <- affy_go%>% left_join(flat_go,by="category") 
p2 <- ggplot(allgo,aes(x=under_represented_pvalue.x,y=under_represented_pvalue.y,col=log(numDEInCat.x)))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("Affymetrix pvalue")+ylab("no bias correction pvalue")+title("Affymetrix")

#pdf("Moriani_bias_vs_nocorrection.pdf")
grid.arrange(p1,p2,ncol=2)
#dev.off()


## in the end, illumina closer to no bias than affymetrix


simugenes=pwf
names(simugenes)=c("selected","genelength","pwf")
simugenes$gene_id=rownames(simugenes)
#simugenes$go_id=go_map


goals=as.data.frame(unlist(go_map))
names(goals)="go_id"
goals$gene_id=rownames(goals)


unlisted=as.data.frame(unlist(go_map))
unlisted$gene_id=substr(rownames(unlisted),1,15)
#rownames(unlisted)=unlisted$gene_id

simugenes_final <- inner_join(unlisted,simugenes,by="gene_id")
names(simugenes_final)[1]="go_id"

lengthbygo=tapply(simugenes_final$genelength,simugenes_final$go_id,mean,na.rm=T)
lengthbygo2=tapply(simugenes_final$genelength,simugenes_final$go_id,median,na.rm=T)
lengthbygosd=tapply(simugenes_final$genelength,simugenes_final$go_id,sd,na.rm=T)
lbgo=as.data.frame(lengthbygo)
names(lbgo)[1]="genelength"
lbgo$category=rownames(lbgo)
lbgo$medgenelength=lengthbygo2
lbgo$sdgenelength=lengthbygosd



source("dotest.R")
simugo=data.frame("go_id"=unique(simugenes_final$go_id))
parameters=list()
simugenes=simugenes_final[!is.na(simugenes_final$genelength),]
results=dotestnewdf(simugo,simugenes)
results2=dotestPermdf(simugo,simugenes)
gofinal=results[[1]]
simugo2=results[[2]]
#parameters=

names(gofinal)[1]="category"
plot(gofinal[,2],gofinal[,3])
head(simugo2)
head(simugo)

head(illumina_go)

nbgenesgo=table(simugenes_final$go_id)

nbggo=data.frame("nbgenes"=as.vector(nbgenesgo),"category"=names(nbgenesgo))

allgo <- illumina_go%>% left_join(gofinal,by="category") %>% left_join(lbgo,by="category") %>% left_join(affy_go,by="category") %>% left_join(nbggo,by="category")

#simugenes%>%
save(allgo,file="allgomarioni.Rdata")
load("allgomarioni.Rdata")

p1 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval1,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+xlim(0,0.05)+ylim(0,0.05)+scale_color_gradientn(colours=rainbow(5))
p2 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval2,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+xlim(0,0.05)+ylim(0,0.05)+scale_color_gradientn(colours=rainbow(5))
grid.arrange(p1,p2,ncol=2)


p1 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval1,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))
p2 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval2,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))
p3 <- ggplot(allgo,aes(x=pval1,y=pval2,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))
grid.arrange(p1,p2,p3,ncol=3)## 
### zoom in the significant corner

p1 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval1,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("elsa pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))+xlim(0,0.10)+ylim(0,0.10)+geom_abline(slope=0,intercept=0.05)+geom_vline(xintercept=0.05)
p2 <- ggplot(allgo,aes(x=over_represented_pvalue.x,y=pval2,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))+xlim(0,0.10)+ylim(0,0.10)+geom_abline(slope=0,intercept=0.05)+geom_vline(xintercept=0.05)
p3 <- ggplot(allgo,aes(x=pval1,y=pval2,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("elsa pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))+xlim(0,0.10)+ylim(0,0.10)+geom_abline(slope=0,intercept=0.05)+geom_vline(xintercept=0.05)
grid.arrange(p1,p2,p3,ncol=3)


#pdf("marioni_newmethod.pdf")
grid.arrange(p1,p2,p3,ncol=3)
#dev.off()

p1 <- ggplot(allgo,aes(x=over_represented_pvalue.y,y=over_represented_pvalue.x,col=sdgenelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))
p2 <- ggplot(allgo,aes(x=over_represented_pvalue.y,y=over_represented_pvalue.x,col=genelength))+geom_point()+geom_smooth(col="green")+geom_abline(slope=1,intercept=0,col="red")+xlab("illumina pvalue")+ylab("no bias correction pvalue")+title("illumina")+scale_color_gradientn(colours=rainbow(5))
grid.arrange(p1,p2,ncol=2)


## differences

pdf("pvaldiff_genelength_marioni.pdf")
p1 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval1,x=log(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="black")+xlab("log(gene length)")+ylab(" pvalue difference")+ggtitle("GOSEQ - New method ")+scale_color_gradientn(colours=rainbow(5))+theme_minimal()
p2 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval2,x=log(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="black")+xlab("log(gene length)")+ylab("pvalue difference")+ggtitle("GOSEQ - no correction")+scale_color_gradientn(colours=rainbow(5))+theme_minimal()
p3 <- ggplot(allgo,aes(y=pval1-pval2,x=log(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="black")+xlab("log(gene length)")+ylab("pvalue difference")+ggtitle("New method - no correction")+scale_color_gradientn(colours=rainbow(5))+theme_minimal()
grid.arrange(p1,p2,p3,ncol=3)
dev.off()

## pval vs genelengvth

pdf("pval_genelength_marioni.pdf",height=5,width=12)

p1 <- ggplot(allgo %>% filter(nbgenes>50) ,aes(y=over_represented_pvalue.x,x=log(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="red")+xlab("log(gene length)")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
p2 <- ggplot(allgo%>% filter(nbgenes>50),aes(y=pval2,x=log(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="red")+ylab("pvalue  ")+xlab("log(gene length)")+ggtitle(" No correction")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
p3 <- ggplot(allgo %>% filter(nbgenes>50),aes(y=pval1,x=(genelength),col=log(nbgenes)))+geom_point()+geom_smooth(col="red")+ylab("pvalue")+xlab("log(gene length)")+ggtitle("new method")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
grid.arrange(p2,p1,p3,ncol=3)

dev.off()

pdf("pval_nbgenes_marioni.pdf",height=5,width=12)
p4 <- ggplot(allgo,aes(y=over_represented_pvalue.x,x=log(nbgenes),col=log(genelength)))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+xlab("log(nb genes)")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
p5 <- ggplot(allgo,aes(y=pval2,x=log(nbgenes),col=log(genelength)))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+ylab("pvalue  ")+xlab("log(nb genes)")+ggtitle(" No correction")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
p6 <- ggplot(allgo,aes(y=pval1,,x=log(nbgenes),col=log(genelength)))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+ylab("pvalue")+xlab("log(nb genes)")+ggtitle("new method")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
grid.arrange(p5,p4,p6,ncol=3)
dev.off()

pdf("pval_nbgenes_marioni2.pdf",height=5,width=12)

p5 <- ggplot(allgo,aes(y=pval2,x=log(nbgenes),col=log(genelength)))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+ylab("pvalue  ")+xlab("log(nb genes)")+ggtitle(" No correction")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()
p6 <- ggplot(allgo,aes(x=log(genelength),y=log(nbgenes),col=pval2))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+ylab("log(nbgenes)")+xlab("log(mean(gene length))")+ggtitle("Counding effect of gene length and GO size")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+theme_minimal()+labs(col="pvalue")
grid.arrange(p5,p6,ncol=2)

dev.off()

## with facet wrap on nb genes
p1 <- ggplot(allgo,aes(y=over_represented_pvalue.x,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+xlab("log(gene length)")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))+theme_minimal()+facet_wrap(~log(nbgenes)%/%1,ncol=5)
p2 <- ggplot(allgo,aes(y=pval2,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue  ")+xlab("log(gene length)")+ggtitle(" No correction")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))+theme_minimal()+facet_wrap(~log(nbgenes)%/%1,ncol=5)
p3 <- ggplot(allgo,aes(y=pval1,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue")+xlab("log(gene length)")+ggtitle("new method")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))+theme_minimal()+facet_wrap(~log(nbgenes)%/%1,ncol=5)

#ggplot(allgo,aes(y=over_represented_pvalue.x,x=log(nbgenes),col=log(genelength)))+geom_point(pch="+",cex=2)+geom_smooth(col="red")+xlab("log(nb genes)")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=colorRampPalette(c("black", "green"))(10))+ylim(c(0,1))+theme_minimal()

allgo %>% filter(nbgenes<10) %>% dim
allgo %>% dim
    
minigos <- allgo %>% filter(nbgenes<10)
p1 <- ggplot(minigos,aes(y=over_represented_pvalue.x,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+xlab("gene length")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
p2 <- ggplot(minigos,aes(y=pval2,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue  ")+xlab("gene length")+ggtitle(" No correction")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
p3 <- ggplot(minigos,aes(y=pval1,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue")+xlab("gene length")+ggtitle("Elsa")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
grid.arrange(p1,p2,p3,ncol=3)

pdf("pval_genelength_bigset_marioni.pdf")
biggos <- allgo %>% filter(nbgenes>20)
p1 <- ggplot(biggos,aes(y=over_represented_pvalue.x,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+xlab("gene length")+ylab("pvalue ")+ggtitle("Goseq")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
p2 <- ggplot(biggos,aes(y=pval2,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue  ")+xlab("gene length")+ggtitle(" No correction")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
p3 <- ggplot(biggos,aes(y=pval1,x=log(genelength),col=log(nbgenes)))+geom_point(pch="+",cex=2)+geom_smooth(col="black")+ylab("pvalue")+xlab("gene length")+ggtitle("Elsa")+scale_color_gradientn(colours=rainbow(5))+ylim(c(0,1))
grid.arrange(p1,p2,p3,ncol=3)
dev.off()

p1 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval1,x=log(genelength),col=medgenelength/genelength))+geom_point(pch="+",cex=2)+geom_smooth(col="green")+xlab("gene length")+ylab("pvalue difference")+ggtitle("Goseq - Elsa")+scale_color_gradientn(colours=rainbow(5))
p2 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval2,x=log(genelength),col=medgenelength/genelength))+geom_point(pch="+",cex=2)+geom_smooth(col="green")+ylab("pvalue difference ")+xlab("gene length")+ggtitle("Goseq - No correction")+scale_color_gradientn(colours=rainbow(5))
p3 <- ggplot(allgo,aes(y=pval1-pval2,x=log(genelength),col=medgenelength/genelength))+geom_point(pch="+",cex=2)+geom_smooth(col="green")+ylab("pvalue difference")+xlab("gene length")+ggtitle("Elsa - no correction")+scale_color_gradientn(colours=rainbow(5))
grid.arrange(p1,p2,p3,ncol=3)



p1 <- ggplot(allgo,aes(y=log(nbgenes),x=log(genelength),col=(pval1<0.05)))+geom_point(pch="+",cex=4)+geom_smooth(col="green")+xlab("nb genes")+ylab("log(genelenght)")+ggtitle("new - Elsa")
p2 <- ggplot(allgo,aes(y=log(nbgenes),x=log(genelength),col=(pval2<0.05)))+geom_point(pch="+",cex=4)+geom_smooth(col="green")+ylab("pvalue difference ")+xlab("gene length")+ggtitle("Goseq - No correction")
p3 <- ggplot(allgo,aes(y=log(nbgenes),x=log(genelength),col=(over_represented_pvalue.x<0.05)))+geom_point(pch="+",cex=4)+geom_smooth(col="green")+ylab("pvalue difference")+xlab("gene length")+ggtitle("Elsa - no correction")
grid.arrange(p1,p2,p3,ncol=3)


p1 <- ggplot(allgo,aes(y=log(nbgenes),x=(pval1<0.05)))+geom_boxplot()
p2 <- ggplot(allgo,aes(y=log(nbgenes),x=(pval2<0.05)))+geom_boxplot()
p3 <- ggplot(allgo,aes(y=log(nbgenes),x=(over_represented_pvalue.x<0.05)))+geom_boxplot()
p4 <- ggplot(allgo,aes(y=log(genelength),x=(pval1<0.05)))+geom_boxplot()
p5 <- ggplot(allgo,aes(y=log(genelength),x=(pval2<0.05)))+geom_boxplot()
p6 <- ggplot(allgo,aes(y=log(genelength),x=(over_represented_pvalue.x<0.05)))+geom_boxplot()
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)

lm1=

pdf("diff_pvalue.pdf")
p1 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval1,x=log(medgenelength),col=log(nbgenes))+geom_point(pch="+",cex=4)+geom_smooth(col="black")+xlab("gene length")+ylab("pvalue difference")+ggtitle("Goseq - Elsa")+scale_color_gradientn(colours=rainbow(5))+ylim(c(-1,1))
p2 <- ggplot(allgo,aes(y=over_represented_pvalue.x-pval2,x=log(medgenelength),col=log(nbgenes)))+geom_point(pch="+",cex=4)+geom_smooth(col="black")+ylab("pvalue difference ")+xlab("gene length")+ggtitle("Goseq - No correction")+scale_color_gradientn(colours=rainbow(5))+ylim(c(-1,1))
p3 <- ggplot(allgo,aes(y=pval1-pval2,x=log(medgenelength),col=log(nbgenes)))+geom_point(pch="+",cex=4)+geom_smooth(col="black")+ylab("pvalue difference")+xlab("gene length")+ggtitle("Elsa - no correction")+scale_color_gradientn(colours=rainbow(5))+ylim(c(-1,1))
grid.arrange(p1,p2,p3,ncol=3)
dev.off()

ggplot(allgo,aes(x=medgenelength,y=genelength,col=log(nbgenes)))+geom_point()+geom_smooth()+geom_abline(slope=1,intercept=0)

ggplot(allgo,aes(x=medgenelength,y=sdgenelength,col=log(nbgenes)))+geom_point()+geom_smooth()+geom_abline(slope=1,intercept=0)


pdf("genelength_bias_marioni.pdf")
ggplot(data,aes(gl,Illumina.pvalue))+geom_point()+geom_smooth(col="red")+ylab("pvalue")+xlab("log(gene length)")
dev.off()



ggplot(biggos,aes(col=over_represented_pvalue.x,x=log(genelength),y=log(nbgenes)))+geom_point(pch="+",cex=10)+geom_smooth(col="black")+xlab("gene length")+ylab("log(nb genes)")+ggtitle("Goseq")+scale_color_gradientn(colours=rainbow(5))+facet_wrap(~((over_represented_pvalue.x*10)%/%1))#+ylim(c(0,1))

library(dplyr)
library(reshape2)
library(plotly)
data_xyz <- biggos %>% dplyr::select(genelength,nbgenes,over_represented_pvalue.x)
#data_z <- acast(data_xyz, genelength~nbgenes,value.var="over_represented_pvalue.x")
plot_ly(data_xyz,x=~log(nbgenes),y=~genelength,z=~over_represented_pvalue.x,color=~over_represented_pvalue.x)


#data_xyz <- biggos %>% dplyr::select(genelength,nbgenes,over_represented_pvalue.x)
plot_ly(biggos,x=~log(nbgenes),y=~genelength,z=~pval1,color=~pval1,mode="lines",group=~pval1)
