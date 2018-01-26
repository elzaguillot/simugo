## script to simulate GO and test with new method


## empirical values from humans
meanGLLT=7.06
sdGLLT=0.866
nbsimu=100000
nbgenes=10000
nbgo=100
goID=paste("GO",1:nbgo,sep=":")
geneID=paste("GENE",1:nbgenes,sep=":")



dosimu<-function(k,meanGLLT,sdGLLT,nbsimu,nbgenes,nbgo,goID,geneID,model)
{
    library(MASS)
    saveOutput=rep(0,15)
    simugenes=data.frame("gene_id"=geneID,"genelength"=rep(0,nbgenes ))
    simugenes$go_id=sample(goID,nbgenes,replace=T)
    goterms=unique(simugenes$go_id)
    rownames(simugenes)=simugenes$gene_id
    simugenes$selected=rep(0,nbgenes)
    ## Random choice of paramters
    ## proportion of go selected vs non selected
    nbgoselected=runif(1,0,0.1)
    ## linear factor between probability of being selected and gene length
    advantage=runif(1,0.5,1)
    goforce=runif(1,0,0.5)
    goselected=rbinom(length(goterms),1,nbgoselected)
    ## save for later
    saveOutput[3]=nbgoselected
    saveOutput[4]=advantage
    saveOutput[5]=goforce
    ## create a fake distribution of mean of each go
    gogldistribution=rnorm(length(goterms),meanGLLT,sdGLLT)
    ## simulated go datasetm, ID , SELECTED, GENE LENGTH
    simugo=data.frame("go_id"=goterms,"selected"=goselected,"gl"=gogldistribution)
    ##
    ## Simulate a pvalue for each gene 
    dividor1=(max(simugenes$genelength)*advantage)/(1-goforce)
    dividor2=(max(exp(simugenes$genelength))*advantage)/(1-goforce)/2
    for(i in 1:nbgenes)
    {
        ## find the go
        thisgo=simugenes$go_id[i][1]
        tgo=which(simugo$go_id==thisgo)
        ## simulate gene length of the gene depending on GO mean GL
        simugenes$genelength[i]=rnorm(1,simugo[tgo,3],sd=sdGLLT)
        ## is that go selected ?
        sel=simugo$selected[tgo]        
        ##
        if(sel==1) ## if yes            
        {
            if(floor(model/2.0)==0)
            {
                simugenes$selected[i]=rbinom(1,1,prob=min(goforce+simugenes$genelength[i]/dividor1*advantage,1))
            }
            if(floor(model/2.0)==1)
                simugenes$selected[i]=rbinom(1,1,prob=min(goforce+exp(simugenes$genelength[i])/dividor2*advantage,1))

        }
        else
        {
            if((model%%2)==0)
                simugenes$selected[i]=rbinom(1,1,prob=min(0+simugenes$genelength[i]/dividor1*advantage/20,1))
            if((model%%2)==1)
                simugenes$selected[i]=rbinom(1,1,prob=0.05)
        }
    }
    ## save for each go terms: ID, pval1, pval2, meanGL    
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    ## set appart the distro of gene length
    allGeneLL=simugenes$genelength
    ## count nb of selected and not selected in this simulated dataset
    TOTnbpositive=sum(simugenes$selected)
    TOTnbnegative=length(simugenes$selected)-TOTnbpositive
    ## start test pval
    for( i in 1:length(gofinal[,1]))
    {
        go=gofinal[i,1]
        genego=simugenes[simugenes$go_id==go,]
        geneLL=genego$genelength
        nbgenesGO=length(geneLL)
        ctable=table(genego$selected)
        selected=max(ctable["1"],0,na.rm=T)
        notselected=max(ctable["0"],0,na.rm=T)
        ## first co;pute the distribubtion
        normalfitGO <- fitdistr(geneLL,"normal")$estimate
        gofinal[i,4]=mean(geneLL)
        normalfitTot <- fitdistr(allGeneLL,"normal")$estimate
        ##        print(paste(selected,notselected,selected*1.0/notselected,go,sep=" ",collapse=" "))
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])            
        ##        print("2")
        ##        print(paste(normalfitGO[1],normalfitGO[2]))
        ##            testsample=sample(1:length(genes[,1]),length(genes[,1])*10,replace=T,prob=dnorm(genes$genelength,normalfitGO[1],normalfitGO[2])*1/d)
        testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=T,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)   
        nbpositive=sum(simugenes$selected[testsample])
        nbnegative=length(simugenes$selected[testsample])-nbpositive        
        ##        print("3")
        notGOpos=nbpositive
        notGOneg=nbnegative
        ##            if(notGOneg>0)
        ##                if(notGOneg>0)
        ##                    contingencytable <- (matrix(c(selected[[1]],notselected[[1]],notGOpos[[1]],notGOneg[[1]]),nrow=2))
        ##     contingencytable <- (c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenes))
        ##        print(paste("treating go",go))
        ##        print(c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO))
        ##                    f=fisher.test(contingencytable,alternative="greater")
        f=phyper(selected[[1]]-1,notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
        ##            print(f$p.value)            
        ##            Q1$pvalue[i]=f$p.value
        gofinal[i,2]=f
        f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
        gofinal[i,3]=f2
    }
    ## then without
    ##    success[k,5]=cor(gofinal[i,2],gofinal[i,4])
    ##    success[k,6]=cor(gofinal[i,2],gofinal[i,4]) 
    ##    print(1)
    saveOutput[1]=table((gofinal[,2]<0.05)==simugo[,2])["TRUE"] ## accuracy 1
    saveOutput[2]=table((gofinal[,3]<0.05)==simugo[,2])["TRUE"] ## accuracy 2
    saveOutput[6]=cor(gofinal[,2],gofinal[,4]) # corellation GL/pvalue
    saveOutput[7]=cor(gofinal[,3],gofinal[,4]) #corellation GL/pvalue
    saveOutput[8]=max(table((gofinal[,2]<0.05)-simugo[,2])["-1"],0,na.rm=T) #FN 1
    saveOutput[9]=sum(((gofinal[,2]<0.05)==F)&(simugo[,2]=="0")) #TN1
    saveOutput[10]=sum(((gofinal[,2]<0.05)==T)&(simugo[,2]=="1")) #TP1
    saveOutput[11]=max(table((gofinal[,2]<0.05)-simugo[,2])["1"],0,na.rm=T)# FP2
    saveOutput[12]=max(table((gofinal[,3]<0.05)-simugo[,2])["-1"],0,na.rm=T) #FN2
    saveOutput[13]=sum(((gofinal[,3]<0.05)==F)&(simugo[,2]=="0")) #TN2
    saveOutput[14]=sum(((gofinal[,3]<0.05)==T)&(simugo[,2]=="1")) #TP2
    saveOutput[15]=max(table((gofinal[,3]<0.05)-simugo[,2])["1"],0,na.rm=T) #FP2
    ##    print(paste("simulation number",k))
    ##    print(table((gofinal[,2]<0.05)-simugo[,2]))
    ##    print(table((gofinal[,3]<0.05)-simugo[,2]))
    return(saveOutput)
}

#dosimu(1,meanGLLT,sdGLLT,nbsimu,nbgenes,nbgo,goID,geneID,1)


#library(parallel)
#makeCluster(2)


library(foreach)
library(doMC)
#registerDoMC(cores=2)


#total=foreach(k = 1:20,.combine=rbind) %dopar% dosimu(1,meanGLLT,sdGLLT,nbsimu,nbgenes,nbgo,goID,geneID)



#library(doMC)
registerDoMC(cores=10)
meanGLLT=7.06
sdGLLT=0.866
nbsimu=1000
nbgenes=10000
for(nbgo in c(500,1000))
{        
    for (nbgenes in c(2000,10000,50000))
        for(model in 0:4)            
        {
            goID=paste("GO",1:nbgo,sep=":")
            geneID=paste("GENE",1:nbgenes,sep=":")
            total=foreach(k = 1:nbsimu,.combine=rbind) %dopar% dosimu(1,meanGLLT,sdGLLT,nbsimu,nbgenes,nbgo,goID,geneID,model)
            names(total)=c("TPTN1","TPTN2","nbgoselected","advantage","goforce","cor1","cor2","FN1","TP1","TN1","FP1","FN2","TP2","TN2","FP2")
            save(total,file=paste("/scratch/cluster/monthly/eguillot/simu7_ge",nbgenes,"_go",nbgo,"_M",model,".Rdata",sep=""))
            write.table(total,file=paste("/scratch/cluster/monthly/eguillot/simu7_ge",nbgenes,"_go",nbgo,"_M",model,".txt",sep=""))
        }
    }

