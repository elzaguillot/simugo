% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simuExp.R
\name{simuExp}
\alias{simuExp}
\title{simuExp}
\usage{
simuExp(geneset, meanBias, sdBias, distribution, advantage, goforce,
  nbgoactive)
}
\arguments{
\item{geneset}{(dataframe) contains geneid and goid}

\item{meanBias}{(double)  mean of the element (e.g. genelength) that create a bias}

\item{sdBias}{(double) standard deviation of the element (e.g. genelength) that create a bias}

\item{advantage}{(double) proportional to the probability that a gene is active if a go is active}

\item{goforce}{(double) the probability that a gene is active is proportional to the bias element and go force}

\item{nbgoactive}{(double) proportion of the go that will randomly selected as active following a binomial distribution}
}
\value{
Returns a list of 3 tables: simugs, simugenes and parameters
}
\description{
simuExp simulate the activity of go and genes
}
\details{
This function randomly simulates a certain proportion of go as active,
and another proportion as inactive. It then randomly simulates genes as
active or inactive, with a probability depending on whether they belong to
an active or inactive gene set

P(gene active| geneset active)=advantage + goforce*biaselement

P(gene active| geneset NOT active)= goforce*biaselement
}
