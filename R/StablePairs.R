
#' Identification of population-level differentially expressed genes in one-phenotype data
#'
#' The StablePairs is to identify stable gene pairs by using the gene expression profile list of the normal samples of one human tissue-type.
#'
#' @param expdata_list The gene expression profile list of the normal samples of
#'  one human tissue-type which is measured by one or more platforms. The first
#'  columns is the gene IDs and the remaining columns are the gene expression
#'  values of the normal samples.
#' @param freq The criteria for identifying stable gene pairs in normal samples.
#'  The default setting of freq is 0.99.
#'
#' @return matrix with each row indicating a stable pair of
#'  gene IDs across the datasets in \code{expdata_list}.
#' @export
#'
#' @examples
#'
#' data(example)
#' expdata_list <- list(normalexp1,normalexp2)
#' stable_pair <- StablePairs(expdata_list)
#'
StablePairs  <-  function(expdata_list, freq = 0.99){
  n <- length(expdata_list)
  data1 <- expdata_list[[1]]
  gid1 <- as.matrix(data1[,1])

  # row names not needed and double memory usage of freqs
  expdata_list <- lapply(expdata_list, function(x) {row.names(x) <- NULL; x})


  # one-time generation of pairs if all in same order
  gids <- lapply(expdata_list, function(x) x[,1])
  gids_eq <- all(sapply(gids, identical, gids[[1]]))

  pairs <-  t(combn(gid1,2))
  data1 <- as.matrix(data1[,-1])

  freqs=list()
  for(i in 1:(length(gid1)-1)){
    cat(i,"\n")
    gid11 <- gid1[-c(1:i)]
    pair1 <- cbind(gid1[i],gid11)
    coms <- data1[match(pair1[,1],gid1),,drop=F]-data1[match(pair1[,2],gid1),,drop=F]
    freq1 <- rowMeans(coms>0)
    freqs[[i]] <- freq1
  }
  freqs <- unlist(freqs)
  stable_pair <- rbind(pairs[freqs>freq,],pairs[freqs<(1-freq),c(2,1)])
  if(n==1){
    return(stable_pair)
  }
  if(n>1){
    # speed up merges
    stable_pair <- data.table::data.table(stable_pair)

    for (k in 2:n){
      data1 <- as.matrix(expdata_list[[k]])
      gid1 <- as.matrix(data1[,1])
      data1 <- as.matrix(data1[,-1])
      freqs <- list()
      for(i in 1:(length(gid1)-1)){
        cat(i,"\n")
        gid11 <- gid1[-c(1:i)]
        pair1 <- cbind(gid1[i],gid11)
        coms <- data1[match(pair1[,1],gid1),,drop=F]-data1[match(pair1[,2],gid1),,drop=F]
        freq1 <- rowMeans(coms>0)
        freqs[[i]] <- freq1
      }
      if (!gids_eq) pairs <- t(combn(gid1,2))
      freqs <- unlist(freqs)
      stablepair <- rbind(pairs[freqs>freq,],pairs[freqs<(1-freq),c(2,1)])

      stablepair <- data.table::data.table(stablepair)
      stable_pair <- merge(stable_pair,stablepair)

    }
    return(as.matrix(stable_pair))
  }
}
