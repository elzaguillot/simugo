#' analyzesimu
#'
#' `analyzesimu` - compare the reel simulated activity to the on inferred on the simulated dataset
#'
#' This function compare the actual simulated activity in gene set, to the ones inferred from the hypergemotric test on the genes
#'
#' @param experiement tba

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
