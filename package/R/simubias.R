#' simubias
#'
#' `simubias` - simulate a biased covariant in the genes
#'
#' This function simulate a biased covariant in the genes that will affect the simulated activity. This bias can follow several distributions.
#'
#'@param geneset table with geneid and goid
#'@param meanBias the mean of the covariant
#'@param sdBias the standard deviation of the covariant
#'@param goactive the percentage of go that is simulated as active
#'@param distribution the distribution of the covariant

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
