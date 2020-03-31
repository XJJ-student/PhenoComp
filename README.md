# PhenoComp
PhenoComp is an algorithm to identify population-level differential genes in one-phenotype data. This algorithm is based on RankComp, an algorithm used to identify individual-level differentially expressed genes in each sample.
# Install
To install the PhenoComp, install from github using devtools
```
library(devtools)
install_github("XJJ-student/PhenoComp")
```
Or you can download the .ZIP file and unzip it.
```
install.packages("PhenoComp",repos = NULL,type="source")
#The "PhenoComp" should be combined with the absolute path.
```
# Usage
StablePairs(expdata_list,freq)
Arguments|Description
:--|:---
expdata_list|The gene expression profile list of the normal samples of one human tissue-type which is measured by one or more platforms. The first columns is the Entrez gene IDs and the remaining columns are the gene expression values of the normal samples.
freq|The criteria for identifying stable gene pairs in normal samples. The default setting of freq is 0.99.

PhenoComp (stable_pair,expE,individualfdr,method,populationfdr)
Arguments|Description
:--|:---
stable_pair|Stable gene pairs identified by the function of StablePairs.
expE|A (non-empty) numeric matrix of disease samples. The first columns is the Entrez gene IDs and the remaining columns are the gene expression values of the disease samples.
individualfdr|The threshold of FDR for identifying individual-level differentially expressed genes.
method|Method determines how to estimate the p_up and p_down. Method=1: the p_up and p_down were estimated as the median values of the frequencies of up-regulated and down-regulated genes for individual disease samples.Method=2: the p_up and p_down were estimated as the mean values of the frequencies of up-regulated and down-regulated genes for individual disease samples.
populationfdr|The threshold of FDR for identifying population-level differentially expressed genes.



# Example
```
library(PhenoComp)
data(example)
expdata_list=list(normalexp1,normalexp2)
stable_pair=StablePairs(expdata_list,0.99)
PhenoComp(stable_pair,expE,0.05,1,0.05)
```
The gene expression data are from the database Gene Expression Omnibus (GEO).        
The **normalexp1** and **normalexp2** are the normal sample expression profiles of the first 1000 genes of GSE26887 and GSE46224, respectively.          
**expE** is the disease sample expression profile of the first 1000 genes of GSE26887.          
**Recommendation**: if the number of genes exceeds 10000, it is recommended to use the server.


# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to:
Jiajing Xie <xiejiajing_fjmu@163.com>; Haidan Yan <Joyan168@126.com>
