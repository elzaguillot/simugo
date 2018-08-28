#' yeastnetwork
#'
#' `yeastnetwork` return a table containing the yeast go network
#'
#' The function return the Saccharomyces cerevisiae go network and its associated genes by
#' by querying the ensemble database fungi.ensembl.org
#'
#' @import biomaRt
#'
#' @return Return a table of two columns of genes and go ensembl id. Each line represent one association between a gene and a geneset



yeastnetwork <- function()
{
if(!require(biomaRt))
{
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
}
## load biomart
mart = useMart("fungi_mart",host="fungi.ensembl.org")
ensembl = useDataset("scerevisiae_eg_gene", mart = mart)
Q = getBM(attributes=c("ensembl_gene_id","go_id"), mart=ensembl)
colnames(Q)=c("gene_id","go_id")
return(Q)
}

