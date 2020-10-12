#' Identification of population-level differentially expressed genes in one-phenotype data
#'
#' PhenoComp is an algorithm to identify population-level differential genes
#' in one-phenotype data. This algorithm is based on RankComp,  an algorithm
#' used to identify individual-level differentially expressed genes in each
#' sample.
#'
#' @param stable_pair Stable gene pairs identified by the function of StablePairs.
#' @param expE A (non-empty) numeric matrix of disease samples.
#'  The first columns is the Entrez gene IDs and the remaining columns
#'  are the gene expression values of the disease samples.
#' @param individualfdr The threshold of FDR for identifying individual-level
#'  differentially expressed genes
#' @param method Method determines how to estimate the p_up and p_down.
#'  Method=1: the p_up and p_down were estimated as the median values of
#'  the frequencies of up-regulated and down-regulated genes for individual
#'  disease samples.Method=2: the p_up and p_down were estimated as the mean
#'  values of the frequencies of up-regulated and down-regulated genes for
#'  individual disease samples.
#' @param populationfdr The threshold of FDR for identifying population-level
#'  differentially expressed genes.
#'
#' @return List with upregulated and downregulated gene IDs.
#' @export
#'
#' @examples
#'
#' data(example)
#' expdata_list=list(normalexp1, normalexp2)
#' stable_pair=StablePairs(expdata_list, 0.99)
#' PhenoComp(stable_pair, expE, 0.05, 1, 0.05)
#'
PhenoComp <- function(stable_pair, expE, individualfdr, method, populationfdr){
  gene <- expE[, 1]
  con_gid <- intersect(unique(c(stable_pair[, 1], stable_pair[, 2])), gene)
  expE <- expE[match(con_gid, gene), -1]
  gene <- con_gid
  index1 <- stable_pair[, 1] %in% gene
  index2 <- stable_pair[, 2] %in% gene
  stable_pair <- stable_pair[index1+index2==2, ]
  case_exp <- expE
  cat("differential expressed analysis... \n")
  outlier_dir <- NULL
  outlier_pvalue <- NULL
  Lgene <- length(gene)
  for (k in 1:Lgene){
    cat(k, "\n")
    colC <- dim(case_exp)[2]
    Nloc_up <- stable_pair[stable_pair[, 1] %in% gene[k], 2]
    Nloc_down <- stable_pair[stable_pair[, 2] %in% gene[k], 1]
    reverse <- matrix(0, colC, 4)
    reverse[, 1] <- rep(length(Nloc_up), colC)
    reverse[, 2] <- rep(length(Nloc_down), colC)
    Tcanc=case_exp[k, ]
    if (length(Nloc_up)>0){
      Nloc_up1 <- match(Nloc_up, gene)
      N_tmp <- matrix(rep(Tcanc, length(Nloc_up)), ncol=colC, byrow=T)-case_exp[Nloc_up1, ]
      case_p <- colSums(N_tmp<0)
      reverse[, 3] <- case_p
    }
    if (length(Nloc_down)>0){
      Nloc_down1 <- match(Nloc_down, gene)
      N_tmpp <- matrix(rep(Tcanc, length(Nloc_down)), ncol=colC, byrow=T)-case_exp[Nloc_down1, ]
      case_pp <- colSums(N_tmpp>0)
      reverse[, 4] <- case_pp;
    }
    GenePair_sig <- NULL
    GenePair <- rep(0, colC)
    GenePair[which(reverse[, 3]>reverse[, 4])] <- -1
    GenePair[which(reverse[, 3]<reverse[, 4])] <-  1
    tmp <- matrix(c(reverse[, 1], reverse[, 2], reverse[, 1]-reverse[, 3]+reverse[, 4],  reverse[, 2]-reverse[, 4]+reverse[, 3]), ncol=4)
    GenePair_sig <- apply(tmp, 1, function(x) fisher.test(matrix(x, ncol=2, byrow=T))$p.value)
    outlier_dir <- rbind(outlier_dir, GenePair)
    outlier_pvalue <- rbind(outlier_pvalue, GenePair_sig)
  }
  fdr <- apply(outlier_pvalue, 2, function(x) p.adjust(x, method="fdr", length(x)))

  individual_level_DEGs <- matrix(0, nrow(outlier_dir), ncol(outlier_dir))
  for(i in 1:ncol(outlier_dir)){
    a <- which(fdr[, i]<individualfdr)
    individual_level_DEGs[a, i] <- outlier_dir[a, i]
  }
  p_diff <- c()
  for(i in 1:ncol(individual_level_DEGs)){
    want_col_DEG <- individual_level_DEGs[, i]
    up_num <- sum(want_col_DEG==1)
    down_num <- sum(want_col_DEG==-1)
    all_num <- nrow(individual_level_DEGs)
    p_diff1<-(up_num+down_num)/all_num
    p_diff <- c(p_diff, p_diff1)
  }
  p_up_result <- c()
  for(i in 1:ncol(individual_level_DEGs)){
    want_col_DEG <- individual_level_DEGs[, i]
    up_num <- sum(want_col_DEG==1)
    down_num <- sum(want_col_DEG==-1)
    p_up_ratio <- up_num/(up_num+down_num)
    p_up_result <- c(p_up_result, p_up_ratio)
  }
  if(method==1){
    p0_up <- (median(p_diff))*(median(p_up_result))
  }
  if(method==2){
    p0_up <- (mean(p_diff))*(mean(p_up_result))
  }
  up_result <- matrix(0, dim(individual_level_DEGs)[1], 4)
  z <- 1
  for(i in 1:dim(individual_level_DEGs)[1]){
    k <- sum(individual_level_DEGs[i, ]==1)
    s <- dim(individual_level_DEGs)[2]
    up_result[z, 2] <- k/s
    up_result[z, 3] <- 1-pbinom(k-1, s, p0_up)
    z <- z+1
  }
  up_result[, 1] <- gene
  up_result[, 4] <- p.adjust(up_result[, 3], method="BH")
  colnames(up_result) <- c("geneid", "up_ratio", "p.value", "FDR")
  up_DEG <- up_result[which(up_result[, 4]<populationfdr), 1]
  cat(length(up_DEG))
  p_down_result <- c()
  for(i in 1:ncol(individual_level_DEGs)){
    want_col_DEG <- individual_level_DEGs[, i]
    up_num <- sum(want_col_DEG==1)
    down_num <- sum(want_col_DEG==-1)
    p_down_ratio <- down_num/(up_num+down_num)
    p_down_result <- c(p_down_result, p_down_ratio)
  }
  if(method==1){
    p0_down <- (median(p_diff))*(median(p_down_result))
  }
  if(method==2){
    p0_down <- (mean(p_diff))*(mean(p_down_result))
  }
  down_result <- matrix(0, dim(individual_level_DEGs)[1], 4)
  z <- 1
  for(i in 1:dim(individual_level_DEGs)[1]){
    k <- sum(individual_level_DEGs[i, ]==-1)
    s <- dim(individual_level_DEGs)[2]
    down_result[z, 2] <- k/s
    down_result[z, 3] <- 1-pbinom(k-1, s, p0_down)
    z <- z+1
  }
  down_result[, 1] <- gene
  down_result[, 4] <- p.adjust(down_result[, 3], method="BH")
  colnames(down_result) <- c("geneid", "down_ratio", "p.value", "FDR")
  down_DEG <- down_result[which(down_result[, 4]<populationfdr), 1]
  cat(length(down_DEG))
  overlap_gene <- intersect(up_DEG, down_DEG)
  overlap_gene_num <- length(overlap_gene)
  cat(overlap_gene_num)
  if(overlap_gene_num==0){
    cat(length(up_DEG)+length(down_DEG))
    population_DEGs <- list(up=up_DEG, down=down_DEG)
    return(population_DEGs)
  }
  if(overlap_gene_num>0){
    up_gene <- up_DEG[-match(overlap_gene, up_DEG)]
    down_gene <- down_DEG[-match(overlap_gene, down_DEG)]
    cat(length(up_gene))
    cat(length(down_gene))
    cat(length(up_gene)+length(down_gene))
    population_DEGs <- list(up=up_gene, down=down_gene)
    return(population_DEGs)
  }
}



