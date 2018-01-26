#script that define function to test for go enrichment according to diffreent methods

dotestdf<-function(simugo,simugenes,method=1)
{## method 1 is global resampling, method 2 is local resampling
    library(MASS)
    ##    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters=simugo$go_id
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    if(sum(simugenes$genelength<0)<=0)
    {
        TOTnbpositive=nbpos_tot(simugenes)
        TOTnbnegative=nbneg_tot(simugenes)
        normalfitTot=fitall(simugenes)
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            genego=simugenes[simugenes$go_id==go,]
            gofinal[i,2]=pgo_sample
            gofinal[i,3]=pgo_basic            
        }
        gofinal[,2]=gofinal[,2]
        gofinal[,3]=gofinal[,3]
        return(list(gofinal,simugo,parameters))
    }
    return(list(gofinal=list(),simugo,parameters))    
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
        ##        if(sum(is.na(positiveGeneLL))==0)
        ##        {
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            ##print(i)
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
                testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)   
                nbpositive=sum(simugenes$selected[testsample])
                nbnegative=length(simugenes$selected[testsample])-nbpositive        
                notGOpos=nbpositive
                notGOneg=nbnegative
                contingencytable <- (matrix(c(selected[[1]],notselected[[1]],notGOpos[[1]],notGOneg[[1]]),nrow=2))
                  f=fisher.test(contingencytable,alternative="greater")
                gofinal[i,2]=f
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2
            }
        }
        return(list(gofinal,simugo,parameters))
    }
    ##    else
    return(list(gofinal=list(),simugo,parameters))    
    ## then without
    ##    success[k,5]=cor(gofinal[i,2],gofinal[i,4])
    ##    success[k,6]=cor(gofinal[i,2],gofinal[i,4]) 
    ##    print(1)
    ##    print(paste("simulation number",k))
    ##    print(table((gofinal[,2]<0.05)-simugo[,2]))    
    ##    print(table((gofinal[,3]<0.05)-simugo[,2]))
}


dotestnew<-function(filename,method=1)
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
        ##        if(sum(is.na(positiveGeneLL))==0)
        ##        {
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:nbgo)
        {
            go=gofinal[i,1]
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                upowered=rep(genego$selected,genego$genelength)
                upoweredgl=rep(genego$genelength,genego$genelength)
                ctable=table(genego$selected)
                ctable2=table(upowered)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                selected2=max(ctable2["1"],0,na.rm=T)*1/mean(genego$genelength)
                notselected2=max(ctable2["0"],0,na.rm=T)*1/mean(genego$genelength)
                ## first compute the distribubtion
                normalfitGO <- fitdistr(geneLL,"normal",na.rm=T)$estimate
                gofinal[i,4]=mean(geneLL)
                ## option1 
                d2=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])
                normalizationFactor=dim(simugenes)[1]/sum(d2) ## total nb of genes / sum of d2  
                notGOpos=sum(simugenes$selected*d2)*normalizationFactor
                notGOneg=sum((!simugenes$selected)*d2)*normalizationFactor
                ## option 2
#                testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)                                
#                nbpositive=sum(simugenes$selected[testsample])
#                nbnegative=length(simugenes$selected[testsample])-nbpositive     
#                notGOpos=nbpositive
#                notGOneg=nbnegative
                f=phyper(selected2[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
                gofinal[i,2]=f
                f2=phyper(selected[[1]],TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2
            }
            gofinal[,2]=p.adjust(gofinal[,2],method = "bonferroni")
            gofinal[,3]=p.adjust(gofinal[,3],method = "bonferroni")
        }
        return(list(gofinal,simugo,parameters))
    }
    ##    else
    return(list(gofinal=list(),simugo,parameters))    
}


dotestnewdf<-function(simugo,simugenes,method=1)
{## method 1 is global resampling, method 2 is local resampling
    library(MASS)
#    load(filename)
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
        ##        if(sum(is.na(positiveGeneLL))==0)
        ##        {
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            ##print(i)
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                upowered=rep(genego$selected,genego$genelength)
                upoweredgl=rep(genego$genelength,genego$genelength)
                ctable=table(genego$selected)
                ctable2=table(upowered)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                selected2=max(ctable2["1"],0,na.rm=T)*1/mean(genego$genelength)
                notselected2=max(ctable2["0"],0,na.rm=T)*1/mean(genego$genelength)
                ## first compute the distribubtion
                normalfitGO <- fitdistr(geneLL,"normal",na.rm=T)$estimate
                gofinal[i,4]=mean(geneLL)
                ## option1 
                d2=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])
                normalizationFactor=dim(simugenes)[1]/sum(d2) ## total nb of genes / sum of d2  
                notGOpos=sum(simugenes$selected*d2)*normalizationFactor
                notGOneg=sum((!simugenes$selected)*d2)*normalizationFactor
                ## option 2
                testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)                                
                nbpositive=sum(simugenes$selected[testsample])
                nbnegative=length(simugenes$selected[testsample])-nbpositive     
                notGOpos=nbpositive
                notGOneg=nbnegative
                f=phyper(selected2[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
                gofinal[i,2]=f*nbgenesGO
                f2=phyper(selected[[1]],TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2*nbgenesGO
            }
        }
        return(list(gofinal,simugo,parameters))
    }
    ##    else
    return(list(gofinal=list(),simugo,parameters))    
}


dotestgoseq<-function(filename,method=1)
{## method 1 is global resampling, method 2 is local resampling
    library(MASS)
    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters
    goterms=simugo$go_id
    gofinal=data.frame("go_id"=goterms,"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    simugenes$go_id=gsub(":","",simugenes$go_id)
    #simugenes$gene_id=gsub(":","",simugenes$gene_id)
    simugo$go_id=gsub(":","",simugo$go_id)
    rownames(simugo)=gsub(":","",simugo$go_id)
    lengthbygo=tapply(simugenes$genelength,simugenes$go_id,mean,na.rm=T)
    lengthbygo2=tapply(simugenes$genelength,simugenes$go_id,median,na.rm=T)
    lengthbygosd=tapply(simugenes$genelength,simugenes$go_id,sd,na.rm=T)
    lbgo=as.data.frame(lengthbygo)
    names(lbgo)[1]="genelength"
    lbgo$go_id=rownames(lbgo)
    nbgo=length(lbgo$go_id)
    lbgo$medgenelength=lengthbygo2
    lbgo$sdgenelength=lengthbygosd    
    assayed.genes <- (simugenes$gene_id)
    de.genes <- simugenes$gene_id[(simugenes$selected==1)]
    gene.vector=as.integer(assayed.genes%in%de.genes)
    pwf1 <- nullp(gene.vector,bias.data=simugenes$genelength) # new_lengthData)
    rownames(pwf1)=simugenes$gene_id
    go_map=data.frame("genes"=simugenes$gene_id,"category"=simugenes$go_id)
    go_map=as.list(tapply(simugenes$go_id,simugenes$gene_id,as.vector))
    goseq_go=goseq(pwf1,gene2cat=go_map)
    names(goseq_go)[1]="go_id"
    mysimu <- gofinal %>% inner_join(goseq_go,by="go_id")
    names(mysimu)[which(names(mysimu)=='over_represented_pvalue')]="pval1"
    if(sum(simugenes$genelength<0)<=0)
    {        
        ## set appart the distro of gene length
        allGeneLL=simugenes$genelength
        ##        if(sum(is.na(positiveGeneLL))==0)
        ##        {
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:nbgo)
        {
            go=gofinal[i,1]
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                ctable=table(genego$selected)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                f2=phyper(selected[[1]],TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2*nbgenesGO
            }
            gofinal[,2]=p.adjust(gofinal[,2],method = "bonferroni")
            gofinal[,3]=p.adjust(gofinal[,3],method = "bonferroni")
        }
        return(list(gofinal,simugo,parameters))
    }
    ##    else
    return(list(gofinal=list(),simugo,parameters))    
}




dotestJensen<-function(filename,method=1)
{## method 1 is global resampling, method 2 is local resampling
    library(MASS)
    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters
    goterms=simugo$go_id
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pval3"=rep(1,length(goterms)))
    if(sum(simugenes$genelength<0)<=0)
    {        
        allGeneLL=simugenes$genelength
        ctable=table(simugenes$selected)
        TOTnbpositive =max(ctable["1"],0,na.rm=T)#
        TOTnbnegative=max(ctable["0"],0,na.rm=T)#
        ## Fit the gene distributiion to a normal curve
        normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
        ## compute the probably of each gene under this model
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.000000000000000000000000000000000000000001
        ## compute the nb of selected / unselected er table
        gotable=table(simugenes$go_id,simugenes$selected)
        nbgo=length(gofinal[,1])
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                ctable=table(genego$selected)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                ## first compute the local go GL distribubtion
                normalfitGO <- fitdistr(geneLL,"normal",na.rm=T)$estimate
                gofinal[i,4]=mean(geneLL) ## save the mean gl of this go
                ## fi the probability under this new distribution
                d2=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])
                normalizationFactor=dim(simugenes)[1]/sum(d2) ## total nb of genes / sum of d2  
                notGOpos=sum(simugenes$selected*d2)*normalizationFactor
                notGOneg=sum((!simugenes$selected)*d2)*normalizationFactor        
                f=phyper(selected[[1]],notGOpos,notGOneg,nbgenesGO,lower.tail=F)
                gofinal[i,2]=f*nbgenesGO
                f2=phyper(selected[[1]],TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2*nbgenesGO
            }
        }
        return(list(gofinal,simugo,parameters))
    }
    return(list(gofinal=list(),simugo,parameters))    
}

dotest2 <- function(filename,method=1)
{
    library(MASS)
    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters
    goterms=simugo[,1]
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    if(sum(simugenes$genelength<0)<=0)
    {
        allGeneLL=simugenes$genelength
        positiveGeneLL=simugenes$genelength[simugenes$selected==1]
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval 
        normalfitPositive <- fitdistr(positiveGeneLL,"normal",na.rm=T)$estimate
        ##    gofinal[i,4]=mean(geneLL)
        normalfitTot <- fitdistr(allGeneLL,"normal",na.rm=T)$estimate
        ##        print(paste(selected,notselected,selected*1.0/notselected,go,sep=" ",collapse=" "))
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.0000000000000000000000000001            
        ##  print("2")
        ##  print(paste(normalfitGO[1],normalfitGO[2]))
        ##  testsample=sample(1:length(genes[,1]),length(genes[,1])*10,replace=T,prob=dnorm(genes$genelength,normalfitGO[1],normalfitGO[2])*1/d)
        testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitPositive[1],normalfitPositive[2])*1/d)
        nbpositive=sum(simugenes$selected[testsample])
        nbnegative=length(simugenes$selected[testsample])-nbpositive
        ##        print("3")
        notGOpos=nbpositive
        notGOneg=nbnegative
        testsamplesimu=simugenes[testsample,]
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            genego=testsamplesimu[testsamplesimu$go_id==go,]
            geneLL=genego$genelength
            nbgenesGO=length(geneLL)
            ctable=table(genego$selected)
            selected=max(ctable["1"],0,na.rm=T)
            notselected=max(ctable["0"],0,na.rm=T)
            f=phyper(selected[[1]]-1,notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
            gofinal[i,2]=f
            f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
            gofinal[i,3]=f2
        }
        p.adjust(gofinal[,2],n)
        p.adjust(gofinal[,3],n)
        return(list(gofinal,simugo,parameters))
    }
    return(list(gofinal=list(),simugo,parameters))
}




dotest3 <- function(filename,method=1)
{
    library(MASS)
    load(filename)
    ## load the following variables
    ## simugenes
    ## simugo
    ## parameters
    goterms=simugo[,1]
    gofinal=data.frame("ID"=goterms,"pval1"=rep(1,length(goterms)),"pval2"=rep(1,length(goterms)),"pva;l3"=rep(1,length(goterms)))
    if(sum(simugenes$genelength<0)<=0)
    {
        ##        if(sum(is.na(positiveGeneLL))==0)
        ##            {
        allGeneLL=simugenes$genelength
        positiveGeneLL=simugenes$genelength[simugenes$selected==1]
        ## count nb of selected and not selected in this simulated dataset
        TOTnbpositive=sum(simugenes$selected)
        TOTnbnegative=length(simugenes$selected)-TOTnbpositive
        ## start test pval 
        normalfitPositive <- fitdistr(positiveGeneLL,"normal",na.rm=T)$estimate
        ##    gofinal[i,4]=mean(geneLL)
        normalfitTot <- fitdistr(allGeneLL,"normal",na.rm=T)$estimate
        ##        print(paste(selected,notselected,selected*1.0/notselected,go,sep=" ",collapse=" "))
        d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])+0.0000000000000000000000000001            
        ##  print("2")
        ##  print(paste(normalfitGO[1],normalfitGO[2]))
        ##  testsample=sample(1:length(genes[,1]),length(genes[,1])*10,replace=T,prob=dnorm(genes$genelength,normalfitGO[1],normalfitGO[2])*1/d)
        testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitPositive[1],normalfitPositive[2])*1/d)
        nbpositive=sum(simugenes$selected[testsample])
        nbnegative=length(simugenes$selected[testsample])-nbpositive
        ##        print("3")
        notGOpos=nbpositive
        notGOneg=nbnegative
        testsamplesimu=simugenes[testsample,]
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            genego=testsamplesimu[testsamplesimu$go_id==go,]
            geneLL=genego$genelength
            nbgenesGO=length(geneLL)
            ctable=table(genego$selected)
            selected=max(ctable["1"],0,na.rm=T)
            notselected=max(ctable["0"],0,na.rm=T)
            f=clogit()
            ##phyper(selected[[1]]-1,notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
            gofinal[i,2]=f
            f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
            gofinal[i,3]=f2
        }
        return(list(gofinal,simugo,parameters))
    }
    return(list(gofinal=list(),simugo,parameters))
}


docompare<-function(gofinal,simugo,parameters)
{
    saveOutput=rep(0,15)
    saveOutput[1]=table((gofinal$pval1<0.05)==simugo$selected)["TRUE"] ## accuracy 1
    saveOutput[2]=table((gofinal$pval2<0.05)==simugo$selected)["TRUE"] ## accuracy 2
    saveOutput[3]=parameters[1]
    saveOutput[4]=parameters[2]
    saveOutput[5]=parameters[3]
    topgoindex= comparetopgos(gofinal,simugo)
    saveOutput[6]=topgoindex[1]  #cor(gofinal$pval1,gofinal[,4]) # corellation GL/pvalue
    saveOutput[7]=topgoindex[2] #corellation GL/pvalue
    saveOutput[8]=max(table((gofinal$pval1<0.05)-simugo$selected)["-1"],0,na.rm=T) #FN 1
    saveOutput[9]=sum(((gofinal$pval1<0.05)==F)&(simugo$selected=="0")) #TN1
    saveOutput[10]=sum(((gofinal$pval1<0.05)==T)&(simugo$selected=="1")) #TP1
    saveOutput[11]=max(table((gofinal$pval1<0.05)-simugo$selected)["1"],0,na.rm=T)# FP2
    saveOutput[12]=max(table((gofinal$pval2<0.05)-simugo$selected)["-1"],0,na.rm=T) #FN2
    saveOutput[13]=sum(((gofinal$pval2<0.05)==F)&(simugo$selected=="0")) #TN2
    saveOutput[14]=sum(((gofinal$pval2<0.05)==T)&(simugo$selected=="1")) #TP2
    saveOutput[15]=max(table((gofinal$pval2<0.05)-simugo$selected)["1"],0,na.rm=T) #FP2
    saveOutput[16]=parameters[4]
    saveOutput[17]=parameters[5]
    return(saveOutput)
}

accuracyofrank <- function(i,simugo,gofinal,j)
{
    top=sum(simugo$selected[sort(gofinal[,j],index.return=T)$ix[1:i]])/(i*1.0)
    }

comparetopgos<-function(gofinal,simugo)
{
                                        #    saveOutput=rep(0,15)
    top1= sapply(1:1000,accuracyofrank,simugo,gofinal,2)## accuracy 1
    top2= sapply(1:1000,accuracyofrank,simugo,gofinal,3)## accuracy 1
    plot(top1,type="l",col="blue",lwd=3)
    lines(top2,type="l",col="red",lwd=3)
    nbgoselected=sum(simugo$selected)    
    if(nbgoselected==0)
        return(0)
    nbgo=length(simugo$selected)
    lines(c(0,1000),c(nbgoselected*1.0/nbgo,nbgoselected*1.0/nbgo),col="grey",lty=2,lwd=2)
    lines(c(nbgoselected,nbgoselected),c(0,1),col="grey",lty=2,lwd=2)
    print(top1[nbgoselected])
    print(top2[nbgoselected])
    return(c(sum(top1-nbgoselected*1.0/nbgo,na.rm=T),sum(top2-nbgoselected*1.0/nbgo,na.rm=T)))
#    return(saveOutput)
}


docompareDF<-function(gotable,simugo)
{
    require(dplyr)
    saveOutput=rep(0,5)
    bigtable=simugo %>%left_join(gotable,by="go_id")
    saveOutput[1]=table((bigtable$over_represented_pvalue<0.05)==bigtable$selected)["TRUE"]    
    ##    saveOutput[1]=table((gotable[,2]<0.05)==simugo$selected)["TRUE"] ## accuracy 1
    saveOutput[2]=max(table((bigtable$over_represented_pvalue<0.05)-bigtable$selected)["-1"],0,na.rm=T) #FN 1
    saveOutput[3]=sum(((bigtable$over_represented_pvalue<0.05)==F)&(bigtable$selected=="0")) #TN1
    saveOutput[4]=sum(((bigtable$over_represented_pvalue<0.05)==T)&(bigtable$selected=="1")) #TP1
    saveOutput[5]=max(table((bigtable$over_represented_pvalue<0.05)-bigtable$selected)["1"],0,na.rm=T)# FP2
    return(saveOutput)
}


docompare2<-function(mysimu,parameters)
{
    saveOutput=rep(0,23)
    saveOutput[1]=table((mysimu$pval1<0.05)==mysimu$selected)["TRUE"] ## accuracy 1
    saveOutput[2]=table((mysimu$pval2<0.05)==mysimu$selected)["TRUE"] ## accuracy 2
    saveOutput[3]=table((mysimu$over_represented_pvalue<0.05)==mysimu$selected)["TRUE"] ## accuracy 3
    ##    saveOutput[4]=cor(mysimu$pval1,mysimu$genelength,use="complete") # corellation GL/pvalue
    ##    saveOutput[5]=cor(mysimu$pval2,mysimu$genelength,use="complete") # corellation GL/pvalue
    ##    saveOutput[6]=cor(mysimu$over_represented_pvalue,mysimu$genelength,use="complete") # corellation GL/pvalue
    saveOutput[7]=max(table((mysimu$pval1<0.05)-mysimu$selected)["-1"],0,na.rm=T) #FN 1
    saveOutput[8]=sum(((mysimu$pval1<0.05)==F)&(mysimu$selected=="0")) #TN1
    saveOutput[9]=sum(((mysimu$pval1<0.05)==T)&(mysimu$selected=="1")) #TP1
    saveOutput[10]=max(table((mysimu$pval1<0.05)-mysimu$selected)["1"],0,na.rm=T)# FP2
    saveOutput[11]=max(table((mysimu$pval2<0.05)-mysimu$selected)["-1"],0,na.rm=T) #FN2
    saveOutput[12]=sum(((mysimu$pval2<0.05)==F)&(mysimu$selected=="0")) #TN2
    saveOutput[13]=sum(((mysimu$pval2<0.05)==T)&(mysimu$selected=="1")) #TP2
    saveOutput[14]=max(table((mysimu$pval2<0.05)-mysimu$selected)["1"],0,na.rm=T) #FP2
    saveOutput[15]=max(table((mysimu$over_represented_pvalue<0.05)-mysimu$selected)["-1"],0,na.rm=T) #FN2
    saveOutput[16]=sum(((mysimu$over_represented_pvalue<0.05)==F)&(mysimu$selected=="0")) #TN2
    saveOutput[17]=sum(((mysimu$over_represented_pvalue<0.05)==T)&(mysimu$selected=="1")) #TP2
    saveOutput[18]=max(table((mysimu$over_represented_pvalue<0.05)-mysimu$selected)["1"],0,na.rm=T) #FP2
    saveOutput[19]=parameters[1]
    saveOutput[20]=parameters[2]
    saveOutput[21]=parameters[3]
    topgoindex=comparetopgos(gofinal,simugo)
    saveOutput[22]= topgoindex[1]#cor(gofinal[,2],gofinal[,4]) # corellation GL/pval
 #   saveOutput[22]=parameters[4]
    saveOutput[23]=topgoindex[2]
    return(saveOutput)
}
