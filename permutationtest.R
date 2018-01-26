## scrip that define a function to test for go enrichment using permutation

## a new testing method by permutation
library(MASS)
dotestPerm <- function(filename)
{
    ## parameters
    load(filename)
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
        nbgenes=length(simugenes$genelength)
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            ## print(i)
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                ctable=table(genego$selected)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                normalfitGO <- fitdistr(geneLL,"normal",rm.na=T)$estimate
                ##                plot(sort(geneLL))
                bootPos=numeric(1000)
                d2=dnorm(simugenes$genelength,normalfitGO[1],normalfitGO[2])
                nbgenesInGo=length(geneLL)
                for(j in 1:1000)
                {
                    ##                       lines(sort(simugenes$genelength[sample(1:length(simugenes$genelength),length(geneLL),prob=dnorm(simugenes$genelength,normalfitGO[1],normalfitGO[2]))]))
                    simulatedGO=simugenes$selected[sample(1:nbgenes,nbgenesInGo,prob=d2)]
                    bootPos[j]=sum(simulatedGO)
                }
                ##                hist(bootPos)
                if(selected>quantile(bootPos,0.95))
                    gofinal[i,2]=0
                else
                    gofinal[i,2]=1                
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2
            }
        }
    }    
    return(list(gofinal,simugo,parameters))
}

dotestPermdf <- function(simugenes,simugo)
{
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
        nbgenes=length(simugenes$genelength)
        for( i in 1:length(gofinal[,1]))
        {
            go=gofinal[i,1]
            print(i)
            genego=simugenes[simugenes$go_id==go,]
            if(dim(genego)[1]>2)
            {
                geneLL=genego$genelength
                nbgenesGO=length(geneLL)
                ctable=table(genego$selected)
                selected=max(ctable["1"],0,na.rm=T)
                notselected=max(ctable["0"],0,na.rm=T)
                normalfitGO <- fitdistr(geneLL,"normal",rm.na=T)$estimate
                ##                plot(sort(geneLL))
                bootPos=numeric(1000)
                d2=dnorm(simugenes$genelength,normalfitGO[1],normalfitGO[2])
                nbgenesInGo=length(geneLL)
                for(j in 1:1000)
                {
                    ##                       lines(sort(simugenes$genelength[sample(1:length(simugenes$genelength),length(geneLL),prob=dnorm(simugenes$genelength,normalfitGO[1],normalfitGO[2]))]))
                    simulatedGO=simugenes$selected[sample(1:nbgenes,nbgenesInGo,prob=d2)]
                    bootPos[j]=sum(simulatedGO)
                }
                ##                hist(bootPos)
                if(selected>quantile(bootPos,0.95))
                    gofinal[i,2]=0
                else
                    gofinal[i,2]=1                
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2
            }
        }
    }    
    return(list(gofinal,simugo,parameters))
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
        nbgo=length(gofinal[,1])
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
                ##  print(paste("treating go",go))
                ##  print(c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO))
                ##  f=fisher.test(contingencytable,alternative="greater")
                f=phyper(selected[[1]]-1,notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
                ##  print(f$p.value)            
                ##  Q1$pvalue[i]=f$p.value
                gofinal[i,2]=f*nbgo
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2*nbgo
            }
            else{
                print(gofinal[i,])
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



dotestJensen<-function(filename,method=1)
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
        nbgo=length(gofinal[,1])
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
                d2=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])
                ##        print(paste(selected,notselected,selected*1.0/notselected,go,sep=" ",collapse=" "))
                ##        print(paste(normalfitGO[1],normalfitGO[2]))
                ##           testsample=sample(1:length(genes[,1]),length(genes[,1])*10,replace=T,prob=dnorm(genes$genelength,normalfitGO[1],normalfitGO[2])*1/d)
                                        #                testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=T,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)
                normalizationFactor=dim(simugenes)[1]/sum(d2)
                nbpositive=sum(simugenes$selected*d2)*normalizationFactor
                nbnegative=sum((!simugenes$selected)*d2)*normalizationFactor        
                ##        print("3")
                notGOpos=nbpositive
                notGOneg=nbnegative
                ## if(notGOneg>0)
                ##   if(notGOneg>0)
                ##     contingencytable <- (matrix(c(selected[[1]],notselected[[1]],notGOpos[[1]],notGOneg[[1]]),nrow=2))
                ##  contingencytable <- (c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenes))
                ##  print(paste("treating go",go))
                ##  print(c(selected[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO))
                ##  f=fisher.test(contingencytable,alternative="greater")
                f=phyper(selected[[1]]-1,notGOpos,notGOneg,nbgenesGO,lower.tail=F)
                ##  print(f$p.value)            
                ##  Q1$pvalue[i]=f$p.value
                gofinal[i,2]=f*nbgo
                f2=phyper(selected[[1]]-1,TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
                gofinal[i,3]=f2*nbgo
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


docompare<-function(gofinal,simugo,parameters)
{
    saveOutput=rep(0,15)
    saveOutput[1]=table((gofinal[,2]<0.05)==simugo$selected)["TRUE"] ## accuracy 1
    saveOutput[2]=table((gofinal[,3]<0.05)==simugo$selected)["TRUE"] ## accuracy 2
    saveOutput[3]=parameters[1]
    saveOutput[4]=parameters[2]
    saveOutput[5]=parameters[3]
    saveOutput[6]=cor(gofinal[,2],gofinal[,4]) # corellation GL/pvalue
    saveOutput[7]=cor(gofinal[,3],gofinal[,4]) #corellation GL/pvalue
    saveOutput[8]=max(table((gofinal[,2]<0.05)-simugo$selected)["-1"],0,na.rm=T) #FN 1
    saveOutput[9]=sum(((gofinal[,2]<0.05)==F)&(simugo$selected=="0")) #TN1
    saveOutput[10]=sum(((gofinal[,2]<0.05)==T)&(simugo$selected=="1")) #TP1
    saveOutput[11]=max(table((gofinal[,2]<0.05)-simugo$selected)["1"],0,na.rm=T)# FP2
    saveOutput[12]=max(table((gofinal[,3]<0.05)-simugo$selected)["-1"],0,na.rm=T) #FN2
    saveOutput[13]=sum(((gofinal[,3]<0.05)==F)&(simugo$selected=="0")) #TN2
    saveOutput[14]=sum(((gofinal[,3]<0.05)==T)&(simugo$selected=="1")) #TP2
    saveOutput[15]=max(table((gofinal[,3]<0.05)-simugo$selected)["1"],0,na.rm=T) #FP2
    saveOutput[16]=parameters[4]
    saveOutput[17]=parameters[5]
    return(saveOutput)
}

tot=matrix(0,nrow=50000,ncol=17)
tot2=matrix(0,nrow=50000,ncol=17)
tot3=matrix(0,nrow=50000,ncol=17)
i=0
for( filename in list.files(path="simuemp/"))
{
    f=paste("simuemp/",filename,sep="")
    print(paste("processing",f))
    print(i)
    i=i+1
    t=dotest(f)
    tot[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    t=dotestJensen(f)
    tot2[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
    t=dotestPerm(f)
    tot3[i,]=as.vector(docompare(t[[1]],t[[2]],t[[3]]))
}


save(tot,tot2,tot3,file="simuperm>Rdata")

n1=1
n2=167
df1=data.frame(tot[n1:n2,])
names(df1)=c("accuracy_new","accuracy_old","go_p_sel","advantage","goforce","cor_new","cor_old","FN_new","TN_new","TP_new","FP_new","FN_old","TN_old","TP_old","FP_old","distr","model")
for( i in c(1:15))
    df1[,i]=as.numeric(tot[n1:n2,i])
df1$FPR_old=df1$FP_old/(df1$FP_old+df1$TN_old)
df1$FPR_new=df1$FP_new/(df1$FP_new+df1$TN_new)
df1$TPR_new=df1$TP_new/(df1$TP_new+df1$FN_new)
df1$TPR_old=df1$TP_old/(df1$TP_old+df1$FN_old)

df2=data.frame(tot2[n1:n2,])
names(df2)=c("accuracy_new","accuracy_old","go_p_sel","advantage","goforce","cor_new","cor_old","FN_new","TN_new","TP_new","FP_new","FN_old","TN_old","TP_old","FP_old","distr","model")
for( i in c(1:15))    
    df2[,i]=as.numeric(tot2[n1:n2,i])
df2$FPR_old=df2$FP_old/(df2$FP_old+df2$TN_old)
df2$FPR_new=df2$FP_new/(df2$FP_new+df2$TN_new)
df2$TPR_new=df2$TP_new/(df2$TP_new+df2$FN_new)
df2$TPR_old=df2$TP_old/(df2$TP_old+df2$FN_old)
df3=data.frame(tot3[n1:n2,])
names(df3)=c("accuracy_new","accuracy_old","go_p_sel","advantage","goforce","cor_new","cor_old","FN_new","TN_new","TP_new","FP_new","FN_old","TN_old","TP_old","FP_old","distr","model")
for( i in c(1:15))    
    df3[,i]=as.numeric(tot3[n1:n2,i])
df3$FPR_old=df3$FP_old/(df3$FP_old+df3$TN_old)
df3$FPR_new=df3$FP_new/(df3$FP_new+df3$TN_new)
df3$TPR_new=df3$TP_new/(df3$TP_new+df3$FN_new)
df3$TPR_old=df3$TP_old/(df3$TP_old+df3$FN_old)
df1$tot=df1$TP_new+df1$TN_new+df1$FP_new+df1$FN_new
df2$tot=df2$TP_new+df2$TN_new+df2$FP_new+df2$FN_new
df3$tot=df3$TP_new+df3$TN_new+df3$FP_new+df3$FN_new

summary(df1)
summary(df3)
summary(df2)

 
plot(df3$accuracy_new*1.0/df3$tot,df1$accuracy_new*1.0/df1$tot,col=rainbow(5)[df1$go_p_sel*16+1],main="difference in UNIVERSE methods - accuracy",xlab="Marc's Universe",ylab="Elsa's Universe")
abline(0,1)
text(800,500,"color represent percentage of selected go in simulations")

plot(df2$FP_new*1.0/df2$tot,df1$FP_new*1.0/df1$tot,col=rainbow(5)[df1$go_p_sel*16+1],main="difference in UNIVERSE methods - false positive",xlab="Marc's Universe",ylab="Elsa's Universe")
abline(0,1)
text(800,500,"color represent percentage of selected go in simulations")

plot(df2$FN_new*1.0/df2$tot,df1$FN_new*1.0/df2$tot,col=rainbow(5)[df1$go_p_sel*16+1],main="difference in UNIVERSE methods - false negative",xlab="Marc's Universe",ylab="Elsa's Universe")
abline(0,1)
text(800,500,"color represent percentage of selected go in simulations")


plot(df3$accuracy_new*1.0/df3$tot,df1$accuracy_new*1.0/df1$tot,col=rainbow(5)[df1$go_p_sel*16+1],main="difference in UNIVERSE methods - accuracy",xlab="Marc's Universe",ylab="Elsa's Universe")
abline(0,1)
text(800,500,"color represent percentage of selected go in simulations")

library(dplyr)
nbgenesgo=table(simugenes$go_id)
nbggo=data.frame("nbgenes"=as.vector(nbgenesgo),"category"=names(nbgenesgo))
allgo <- df1 %>% left_join(df2,by="go_id") %>% left_join(df3,by="go_id") %>%  left_join(nbggo,by="category")
