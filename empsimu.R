## what is the distribution of genelength depending on GO terms

##source("https://bioconductor.org/biocLite.R")
##biocLite("biomaRt")


load(file="go_human_160617.txt")

library(tidyverse)

#filter by the number of genes per GO
Q3=list('go_id'=names(which(table(Q[,2])>10)))


## only keep GO with more than 100 terms
Q3=list('go_id'=names(which(table(Q[,2])>100)))
maxnbgenes=max(table(Q[,2]))
Qclean=Q[Q$go_id %in% unlist(Q3),]




Qtable=as_tibble(Qclean)

Qclean %>% group_by(ensembl_gene_id, go_id) %>% filter(row_number()==1) %>% filter(go_id!="")


gotable=data.frame(tapply(Qclean$transcript_length,Qclean$go_id,mean))
gotable$nbgenes=tapply(Qclean$transcript_length,Qclean$go_id,length)
colnames(gotable)[1]="meangl"
colnames(gotable)[2]="nbgenes"
## poisson log distributed
gotable$rank=rank(gotable$meangl)
    ##sort(gotable$meangl,index.return=T)$ix
gotable$go_id=rownames(gotable)


## remove the outliers with too many genes
gotable=gotable[-which(log(gotable$nbgenes)>8),]
#gotable=gotable[-which(gotable$go_id==" "),]

Qtot=merge(Qclean,gotable,by="go_id")

library(dplyr)
library(MASS)






simugenes=Qtot[,c(1,2,3)]
names(simugenes)[3]="genelength"
simugo=gotable[,c(1,2,4)]

for(model in 0:3)
{
    for(k in 1:1000)
    {
        doempsimunew(1,simugenes,simugo,model,100000)
    }
}

dotest<-function(filename,method=1)
{## method 1 is global resampling, method 2 is local resampling
    library(MASS)
    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters
    goterms=simugo$go_id
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    if(sum(simugenes$genelength<0)<=0)
    {        
        ## set appart the distro of gene length
        allGeneLL=simugenes$genelength
#        if(sum(is.na(positiveGeneLL))==0)
#        {
            ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            #print(i)
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                ctable=table(genego$selected)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                ## first co;pute the distribubtion
                normalfitGO <- fitdistr(geneLL,"normal",na.rm=T)$estimate
                gofinal[i,4]=mean(geneLL)
                ##        print(paste(selected,notselected,selected*1.0/notselected,go,sep=" ",collapse=" "))
                ##        print(paste(normalfitGO[1],normalfitGO[2]))
                ##           testsample=sample(1:length(genes[,1]),length(genes[,1])*10,replace=T,prob=dnorm(genes$genelength,normalfitGO[1],normalfitGO[2])*1/d)
                testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=T,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)   
                nbpositive=sum(simugenes$selected[testsample])
                nbnegative=length(simugenes$selected[testsample])-nbpositive        
            ##        print("3")
                notGOpos=nbpositive
                notGOneg=nbnegative
                ## if(notGOneg>0)
                ##   if(notGOneg>0)
            ##     contingencytable <- (matrix(c(selected[[1]],notselected[[1]],notGOpos[[1]],notGOneg[[1]]),nrow=2))
            ##  contingencytable <- (c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenes))
                ##  f=fisher.test(contingencytable,alternative="greater")
                f=phyper(selected[[1]]-1,notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
                ##  print(f$p.value)            
                ##  Q1$pvalue[i]=f$p.value
                gofinal[i,2]=f
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2
            }
            else{
                                        #                print(gofinal[i,])                
            }
        }
        return(list(gofinal,simugo,parameters))
    }
    ##    else
    return(list(gofinal=list(),simugo,parameters))    
    ## then without
    ##    success[k,5]=cor(gofinal[i,2],gofinal[i,4])
    ##    success[k,6]=cor(gofinal[i,2],gofinal[i,4]) 
   
}

tot=matrix(0,nrow=4000,ncol=17)
tot2=matrix(0,nrow=4000,ncol=17)
tot3=matrix(0,nrow=4000,ncol=17)
#tot2=matrix(0,nrow=50000,ncol=17)
i=0
for( filename in list.files(path="simuemp/"))
{
    print(paste("processing",paste("simuemp/",filename,sep="")))
    print(i)
    i=i+1
    t=dotest(paste("simuemp/",filename,sep=""))
    tot[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    t=dotestJensen(paste("simuemp/",filename,sep=""))
    tot2[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    t=dotestPerm(paste("simuemp/",filename,sep=""))
    tot3[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    save(tot,tot2,tot3,file="emptot.Rdata")
}
#save(tot,file="emptot.Rdata")


filename2=filename
t=dotest(filename2)
