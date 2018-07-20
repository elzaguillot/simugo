## create a fake set of go / genes
## 3 parameters:
## nbgo: total number of unique go
## nbgenes: total number of unique genes
## nbgopergenes: average nb of go that a genes belong to

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
        nbgoforeachgene=min(rpois(1,2),nbgo)
        gene_id=paste0("gene_id",rep(1:nbgenes,nbgoforeachgene))
        go_id2=unlist(lapply(1:nbgenes,y <- function(x)sample(go_id,nbgoforeachgene,replace=F)))
        return(data.frame("go_id"=go_id2,"gene_id"=gene_id))
    }
}


