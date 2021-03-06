% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yeast.R
\name{yeastnetwork}
\alias{yeastnetwork}
\title{yeastnetwork}
\usage{
yeastnetwork()
}
\value{
Return a table of two columns of genes and go ensembl id. Each line represent one association between a gene and a geneset
}
\description{
yeastnetwork() return a table containing the yeast go network
}
\details{
The function return the Saccharomyces cerevisiae go network and its associated genes by
by querying the ensemble database fungi.ensembl.org
}
