simuExp<-function(geneset,meanBias,sdBias, distribution,popGoExp,advantage,goforce,nbgoselected)
{
  library(MASS)
  parameters = rep(0, 5)
  geneId = unlist(geneIds(geneset))
  goId = goIds(geneset)
  nbgo=length(goId)
  nbgenes=  length(geneId)
#  simugenes = data.frame("gene_id" = geneset$gene_id, "bias" = rep(0, nbgenes ),go_id=geneset$go_id)
 # simugenes$exp = rep(0, nbgenes)
  parameters[1] = nbgoselected
  parameters[2] = advantage
  parameters[3] = goforce
  parameters[4] = distribution
  ## create a fake distribution of mean of each go    
  goselected=rbinom(n = nbgo,size = 1,prob = nbgoselected)
  tmp=simubias(geneset,meanBias,sdBias,goselected,distribution)
  simugenes=tmp[[2]]
  simugo=tmp[[1]]
  ##    gobiasdistribution = rgamma(length(goterms), meanBias/5, scale = 5)
  ##   gobiasdistribution = rcauchy(length(goterms), meanBias, sdBias/10)
  ##    print(distribution)
  ## simulated go datasetm,  ID ,  SELECTED,  GENE LENGTH
  ## Simulate a pvalue for each gene
  divisor1<-(max(simugenes$bias)*advantage)/(1-goforce)
  #    for(i in 1:nbgenes)
  ## loop over all goes in data set
  ## check if the go is selected and then decide whether the gene is selected based on the go
  for(g in 1:nbgo)
  {
    ## find the go
    thisgo = simugo$go_id[g]	
    tgenes = which(as.character(simugenes$go_id) == as.character(thisgo))
    ## is that go selected ?
    sel = simugo$selected[g]
    if(sel == 1) ## if go is selected, then binomial probability that the gene is selected, ax+b, goforce+bias*advantage            
    {
      simugenes$selected[tgenes] = rbinom(length(tgenes), 1, prob = min(goforce+simugenes$bias[tgenes]/divisor1*advantage, 1))
    }
    else ## if go is NOS selected, then binomial probability that the gene is selected, ax+b, bias*advantage
    {
      simugenes$selected[tgenes] = rbinom(length(tgenes), 1, prob = min(0+simugenes$bias[tgenes]/divisor1*advantage/2, 1))
    }
    ## difference between the two cases is goforce term added
  }
  uniquegenes = unique(simugenes$gene_id)
  for(t in 1:length(uniquegenes))
  {
    thesegenes = which(simugenes$gene_id == uniquegenes[t])
    if(length(thesegenes)>1)
    {
      simugenes$selected[thesegenes] = sample(simugenes$selected[thesegenes], 1)  ## randomly sample one of the go output
      simugenes$selected[thesegenes] = floor(median(simugenes$selected[thesegenes]))  ## majority rule = take the median, if even number and 0.5 is median, take 1
    }
  }
  simugenes
  #    simugenes <- simugenes %>% distinct(gene_id) #remove the doublons
  selectset = which(simugenes$go_id %in% simugo$go_id[(simugo$selected == 1)]) ## save the genes of go that are sleected
#  save(simugo, simugenes, parameters, file = paste("simusept/simu", paste(parameters[c(4, 5)], collapse = "-"), nbgo, nbgenes, format(Sys.time(),  "%H:%M:%S"), ".Rdata", sep = "-"))
  #return(c(cor(simugenes[selectset, c(2, 4)])[1, 2], parameters))
  return(list(simugo,simugenes,parameters))
}

simubias<-function(geneset,meanBias,sdBias,goselected,distribution="normal")
{
  geneId = unlist(geneIds(geneset))
  goId = goIds(geneset)
  nbgo=length(goId)
  nbgenes=  length(geneId)
  simugenes = data.frame("gene_id" = geneset$gene_id, "bias" = rep(0, nbgenes ),"go_id"=geneset$go_id,'selected'=rep(0, nbgenes ))
  simugenes$selected = rep(0, nbgenes)
  gobiasdistribution = rlnorm(length(goId), meanBias, sdBias/10)
  simugo = data.frame("go_id" = goId, "selected" = goselected, "biasteo" = gobiasdistribution, "biasemp" = 0, "nbgenes" = 0)
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


analyzesimu<-function(simugo,simugenes,parameters)
{## proportion of positive genes in a positive go
    ## proportrion positive genes in a negative go
  posgo=which(simugo$selected==T)
  posgog=simugenes$go_id %in% simugo[posgo,]$go_id
  goselected=length(posgo)/length(simugo$selected)
  if((sum(posgog>0))) # warning, if no go selected, the value is NA, but that is what we want
      {popposgene=mean(simugenes[posgog,]$selected)}
  else
      popposgene=0
  if((sum((!posgog)>0))) # warning, if all go selected, the value is NA, but that is what we want
      popneggene=mean(simugenes[!posgog,]$selected)
  else
      popneggene=0
  tmp=c(popposgene,popneggene,goselected)
#  simugenes[posgog,]
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
