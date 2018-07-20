simuExp<-function(geneset,meanBias,sdBias, distribution,popGoExp,advantage,goforce,nbgoactive)
{
  library(MASS)
  parameters = rep(0, 5)
  geneId = unlist(geneIds(geneset))
  goId = goIds(geneset)
  nbgo=length(goId)
  nbgenes=  length(geneId)
#  simugenes = data.frame("gene_id" = geneset$gene_id, "bias" = rep(0, nbgenes ),go_id=geneset$go_id)
 # simugenes$exp = rep(0, nbgenes)
  parameters[1] = nbgoactive
  parameters[2] = advantage
  parameters[3] = goforce
  parameters[4] = distribution
  ## create a fake distribution of mean of each go    
  goactive=rbinom(n = nbgo,size = 1,prob = nbgoactive)
  tmp=simubias(geneset,meanBias,sdBias,goactive,distribution)
  simugenes=tmp[[2]]
  simugo=tmp[[1]]
  ##    gobiasdistribution = rgamma(length(goterms), meanBias/5, scale = 5)
  ##   gobiasdistribution = rcauchy(length(goterms), meanBias, sdBias/10)
  ##    print(distribution)
  ## simulated go datasetm,  ID ,  ACTIVE,  GENE LENGTH
  ## Simulate a pvalue for each gene
  divisor1<-(max(simugenes$bias)*advantage)/(1-goforce)
  #    for(i in 1:nbgenes)
  ## loop over all goes in data set
  ## check if the go is active and then decide whether the gene is active based on the go
  for(g in 1:nbgo)
  {
    ## find the go
    thisgo = simugo$go_id[g]	
    tgenes = which(as.character(simugenes$go_id) == as.character(thisgo))
    ## is that go active ?
    sel = simugo$active[g]
    if(sel == 1) ## if go is active, then binomial probability that the gene is active, ax+b, goforce+bias*advantage            
    {
      simugenes$active[tgenes] = rbinom(length(tgenes), 1, prob = min(goforce+simugenes$bias[tgenes]/divisor1*advantage, 1))
    }
    else ## if go is NOS active, then binomial probability that the gene is active, ax+b, bias*advantage
    {
      simugenes$active[tgenes] = rbinom(length(tgenes), 1, prob = min(0+simugenes$bias[tgenes]/divisor1*advantage/2, 1))
    }
    ## difference between the two cases is goforce term added
  }
  uniquegenes = unique(simugenes$gene_id)
  for(t in 1:length(uniquegenes))
  {
    thesegenes = which(simugenes$gene_id == uniquegenes[t])
    if(length(thesegenes)>1)
    {
      simugenes$active[thesegenes] = sample(simugenes$active[thesegenes], 1)  ## randomly sample one of the go output
      simugenes$active[thesegenes] = floor(median(simugenes$active[thesegenes]))  ## majority rule = take the median, if even number and 0.5 is median, take 1
    }
  }
  simugenes
  #    simugenes <- simugenes %>% distinct(gene_id) #remove the doublons
  activeset = which(simugenes$go_id %in% simugo$go_id[(simugo$active == 1)]) ## save the genes of go that are sleected
#  save(simugo, simugenes, parameters, file = paste("simusept/simu", paste(parameters[c(4, 5)], collapse = "-"), nbgo, nbgenes, format(Sys.time(),  "%H:%M:%S"), ".Rdata", sep = "-"))
  #return(c(cor(simugenes[selectset, c(2, 4)])[1, 2], parameters))
  return(list(simugo,simugenes,parameters))
}

simubias<-function(geneset,meanBias,sdBias,goactive,distribution="normal")
{
  geneId = unlist(geneIds(geneset))
  goId = goIds(geneset)
  nbgo=length(goId)
  nbgenes=  length(geneId)
  simugenes = data.frame("gene_id" = geneset$gene_id, "bias" = rep(0, nbgenes ),"go_id"=geneset$go_id,'active'=rep(0, nbgenes ))
  simugenes$active = rep(0, nbgenes)
  gobiasdistribution = rlnorm(length(goId), meanBias, sdBias/10)
  simugo = data.frame("go_id" = goId, "active" = goactive, "biasteo" = gobiasdistribution, "biasemp" = 0, "nbgenes" = 0)
  rownames(simugo) = goId
  for(i in 1:nbgenes)
  {
      ## find the go
      thisgo = simugenes$go_id[i][1]
      tgo = which(as.character(simugo$go_id) == as.character(thisgo))
      ## simulate gene length of the gene depending on GO mean bias
      if (distribution == "normal")
          simugenes$bias[i] = max(0.1, rnorm(1, simugo[tgo, 3], sd = sdBias))
      if (distribution == "gamma")
      {
          ##          print("fail")
          simugenes$bias[i] = max(0.1, rgamma(1, simugo[tgo, 3], scale = sdBias))
          ##            gobiasdistribution = rgamma(length(goterms), meanBias/5, scale = 5)
      }
      if (distribution == "cauchy")
      {
          simugenes$bias[i] = max(0.1, rcauchy(1, simugo[tgo, 3], scale = sdBias))
          ##            print(paste(simugenes$bias[i], simugo[tgo, 3]))
      }
      simugo[tgo, ]$biasemp = simugo[tgo, ]$biasemp+simugenes$bias[i]
      simugo[tgo, ]$nbgenes = simugo[tgo, ]$nbgenes+1
  }
  ## compute the average bias per go
  simugo$biasemp = simugo$biasemp*1.0/simugo$nbgenes+1
  return(list(simugo,simugenes))
}


analyzesimu<-function(experiment)
{
    simugo=experiment[[1]]
    simugenes=experiment[[2]]
    parameters=experiment[[3]]
    ## proportion of active genes in a active go
    ## proportrion active genes in a negative go
    posgo=which(simugo$active==T)
    posgog=simugenes$go_id %in% simugo[posgo,]$go_id
    goactive=length(posgo)/length(simugo$active)
    if((sum(posgog>0))) # warning, if no go active, the value is NA, but that is what we want
    {
        popposgene=sum(simugenes[posgog,]$active)
    }
  else
      popposgene=0
  if((sum((!posgog)>0))) # warning, if all go active, the value is NA, but that is what we want
  {popneggene=sum(simugenes[!posgog,]$active)
    }
  else
      popneggene=0
    tmp=c(popposgene,popneggene,goactive)
    ##  simugenes[posgog,]
    names(tmp) <- c("truePositiveSignal","falsePositiveSignal","proportionActiveGo")
  return(tmp)
}

geneIds<-function(geneset)
{
  return(mockgeneset$gene_id)
}


geneIdsU<-function(geneset)
{
  return(unique(geneIds(geneset)))
}

goIds<-function(geneset)
{
  return(unique(mockgeneset$go_id))
}
