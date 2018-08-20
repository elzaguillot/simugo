# function that compare the infered enrichment and the simulated ones


docompareDF<-function(gotable,simugo)
{
  require(dplyr)
  saveOutput=rep(0,5)
  bigtable=simugo %>%left_join(gotable,by="go_id")
  saveOutput[1]=table((bigtable$over_represented_pvalue<0.05)==bigtable$active)["TRUE"]
  ##    saveOutput[1]=table((gotable[,2]<0.05)==simugo$active)["TRUE"] ## accuracy 1
  saveOutput[2]=max(table((bigtable$over_represented_pvalue<0.05)-bigtable$active)["-1"],0,na.rm=T) #FN 1
  saveOutput[3]=sum(((bigtable$over_represented_pvalue<0.05)==F)&(bigtable$active=="0")) #TN1
  saveOutput[4]=sum(((bigtable$over_represented_pvalue<0.05)==T)&(bigtable$active=="1")) #TP1
  saveOutput[5]=max(table((bigtable$over_represented_pvalue<0.05)-bigtable$active)["1"],0,na.rm=T)# FP2
  return(saveOutput)
}

