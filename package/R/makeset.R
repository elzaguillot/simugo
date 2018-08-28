#' makeset
#'
#' Makeset- Create a SIMULATED set of go / genes
#'
#' This function creates a random bipartite graph simulated gene sets
#'
#' @param nbgo integer - total number of unique go
#' @param nbgenes integer - total number of unique genes
#' @param nbgopergenes float - average nb of go that a genes belong to

makeset <- function(nbgo,nbgenes,nbgopergenes)
{
    if(nbgopergenes==1)
    {
        go_id=paste0("GO",1:nbgo)
        gene_id=paste0("gene_id",1:nbgenes)
        return(data.frame("go_id"=sample(go_id,nbgenes,replace=T),"gene_id"=gene_id))
    }
    else
    {
        go_id=paste0("GO",1:nbgo)
                nbgoforeachgene=sapply(rpois(nbgenes,nbgopergenes),min,nbgo)
        gene_id=paste0("gene_id",rep(1:nbgenes,nbgoforeachgene))
        go_id2=unlist(sapply(nbgoforeachgene,sample,x=go_id,replace=F))
        return(data.frame("go_id"=go_id2,"gene_id"=gene_id))
    }
}


