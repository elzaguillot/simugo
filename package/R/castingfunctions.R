
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
