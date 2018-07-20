## this method retrieve a genes and geneset data that represents the real netowrk found in yeast

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

