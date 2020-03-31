\name{StablePairs}
\alias{StablePairs}
\title{Identification of population-level differentially expressed genes in one-phenotype data}
\usage{
StablePairs(expdata_list,freq)
}
\arguments{
  \item{expdata_list}{The gene expression profile list of the normal samples of one human tissue-type which is measured by one or more platforms. The first columns is the Entrez gene IDs and the remaining columns are the gene expression values of the normal samples.}  
  \item{freq}{The criteria for identifying stable gene pairs in normal samples. The default setting of freq is 0.99.}
  }
\description{
The StablePairs is to identify stable gene pairs by using the gene expression profile list of the normal samples of one human tissue-type.
}
\examples{
expdata_list=list(normalexp1,normalexp2)
stable_pair=StablePairs(expdata_list,0.99)
}
