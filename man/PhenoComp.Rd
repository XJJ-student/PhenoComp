\name{PhenoComp}
\alias{PhenoComp}
\title{Identification of population-level differentially expressed genes in one-phenotype data}
\usage{
PhenoComp (stable_pair,expE,individualfdr,method,populationfdr)
}
\arguments{
  \item{stable_pair}{Stable gene pairs identified by the function of StablePairs.}
  \item{expE}{A (non-empty) numeric matrix of disease samples. The first columns is the Entrez gene IDs and the remaining columns are the gene expression values of the disease samples.}
  \item{individualfdr}{The threshold of FDR for identifying individual-level differentially expressed genes}
  \item{method}{Method determines how to estimate the p_up and p_down. Method=1: the p_up and p_down were estimated as the median values of the frequencies of up-regulated and down-regulated genes for individual disease samples.Method=2: the p_up and p_down were estimated as the mean values of the frequencies of up-regulated and down-regulated genes for individual disease samples.
}
  \item{populationfdr}{The threshold of FDR for identifying population-level differentially expressed genes.}
}
\description{
PhenoComp is an algorithm to identify population-level differential genes in one-phenotype data. This algorithm is based on RankComp, an algorithm used to identify individual-level differentially expressed genes in each sample.
}
\examples{
expdata_list=list(normalexp1,normalexp2)
stable_pair=StablePairs(expdata_list,0.99)
PhenoComp(stable_pair,expE,0.05,1,0.05)
}
