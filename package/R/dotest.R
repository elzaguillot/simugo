#' dotest
#'
#' `dotest` - run enrichment test for gene sets
#'
#' This function runs enrichment on simulated or obsevered gene sets, based on binary data (activity or no activity).
#' It uses a hypergeometric test, as well as a corrected hypergeometric test for bias such as gene length.
#'
#' @param adjustPvalMethod correction method, see `p.adjust`
#' @param experiment

## script that define function to test for go enrichment according to diffrent methods

dotest<-function(experiment,adjustPvalMethod="fdr")
{## method 1 is global resampling, method 2 is local resampling
  simugo=experiment[[1]]
  simugenes=experiment[[2]]
  parameters=experiment[[3]]
  library(MASS)
  #    load(filename)
  ## load the following variables
  ## simugenes
  ## simugo
  ## parameters
  goterms=simugo$go_id
  gofinal=data.frame("go_id"=goterms,"pval_classic"=rep(1,length(goterms)),"pval_corrected"=rep(1,length(goterms)),"meanBias"=rep(1,length(goterms)))
  ## set appart the distro of gene length
  allGeneLL=simugenes$bias
  #        print(simugenes)
  #        print(allGeneLL)
  ##        if(sum(is.na(positiveGeneLL))==0)
  ##        {
  ## count nb of active and not active in this simulated dataset
  TOTnbpositive=sum(simugenes$active)
  TOTnbnegative=length(simugenes$active)-TOTnbpositive
  ## start test pval
  normalfitTot <- fitdistr(allGeneLL,"normal",rm.na=T)$estimate
  d=dnorm(allGeneLL,normalfitTot[1],normalfitTot[2])#+0.000000000000000000000000000000000000000001
  # because we divide by this value later, can not be ==0
  gotable=table(simugenes$go_id,simugenes$active)
  nbgo=length(gofinal[,1])
  for( i in seq_along((gofinal[,1])))
      {
      go=gofinal[i,1]
      genego=simugenes[simugenes$go_id==go,]
      if(dim(genego)[1]>2)
        {
        bias=genego$bias
        nbgenesGO=length(bias)
        upowered=rep(genego$active,genego$bias)
        upoweredgl=rep(genego$bias,genego$bias)
        ctable=table(genego$active)
        ctable2=table(upowered)
        active=max(ctable["1"],0,na.rm=T)
        notactive=max(ctable["0"],0,na.rm=T)
        active2=max(ctable2["1"],0,na.rm=T)*1/mean(genego$bias)
        notactive2=max(ctable2["0"],0,na.rm=T)*1/mean(genego$bias)
        ## first compute the distribubtion
        normalfitGO <- fitdistr(bias,"normal",na.rm=T)$estimate
        gofinal[i,4]=mean(bias)
        ## option1
        d2=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])
        normalizationFactor=dim(simugenes)[1]/sum(d2) ## total nb of genes / sum of d2
        notGOpos=sum(simugenes$active*d2)*normalizationFactor
        notGOneg=sum((!simugenes$active)*d2)*normalizationFactor
        ## option 2
        testsample=sample(1:length(simugenes[,1]),length(simugenes[,1]),replace=TRUE,prob=dnorm(allGeneLL,normalfitGO[1],normalfitTot[2])*1/d)
        nbpositive=sum(simugenes$active[testsample])
        nbnegative=length(simugenes$active[testsample])-nbpositive
        notGOpos=nbpositive
        notGOneg=nbnegative
        f=phyper(active2[[1]],notGOpos[[1]],notGOneg[[1]],nbgenesGO,lower.tail=F)
        gofinal[i,2]=p.adjust(f,method=adjustPvalMethod)
        f2=phyper(active[[1]],TOTnbpositive,TOTnbnegative,nbgenesGO,lower.tail=F)
        gofinal[i,3]=p.adjust(f2,method=adjustPvalMethod)
        }
    }
    return(merge(gofinal,simugo,by="go_id"))
  #    }
  ##    else
  #   return(list(gofinal=list(),simugo,parameters))
}

