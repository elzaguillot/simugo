#' simuExp
#'
#' simuExp simulate the activity of go and genes
#'
#' This function randomly simulates a certain proportion of go as active,
#' and another proportion as inactive. It then randomly simulates genes as
#' active or inactive, with a probability depending on whether they belong to
#' an active or inactive gene set
#'
#' P(gene active| geneset active)=advantage + goforce*biaselement
#'
#' P(gene active| geneset NOT active)= goforce*biaselement
#'
#'
#' @param geneset (dataframe) contains geneid and goid
#' @param meanBias (double)  mean of the element (e.g. genelength) that create a bias
#' @param sdBiasGo the standard deviation of the covariant between gene sets (if =0 all go have similar distribution)
#' @param sigma paramater of variation in the distributin of the covariant foe each gene within a go
#' @param advantage (double) proportional to the probability that a gene is active if a go is active
#' @param goforce (double) the probability that a gene is active is proportional to the bias element and go force
#' @param nbgoactive (double) proportion of the go that will randomly selected as active following a binomial distribution
#'
#' @return Returns a list of 3 tables: simugs, simugenes and parameters

simuExp<-function(geneset,meanBias,sdBiasGo,sigma, distribution,advantage,goforce,nbgoactive)
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
  # randomly assign each geneset to active or not active
  goactive=rbinom(n = nbgo,size = 1,prob = nbgoactive)
  # simulate the covariant in genes
  tmp=simubias(geneset,meanBias,sdBiasGo,sigma,goactive,distribution)
  simugenes=tmp[[2]]
  simugo=tmp[[1]]
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
  for(t in seq_along(uniquegenes))
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
