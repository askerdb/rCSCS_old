library(foreach)

distance_cscs <- function(features, css, cosine_threshold = 0.6, weighted = T, verbose = F){
  stopifnot(nrow(features) == nrow(css))
  diag(css) <- 1
  css[css < cosine_threshold] <- 0
  css <- as.matrix(css)
  sample_names <- colnames(features)[-1]
  distm <- matrix(0, nrow = length(sample_names), ncol = length(sample_names),
                  dimnames = list(sample_names, sample_names))
  spn <- combn(sample_names, 2, simplify=FALSE)
  distlist <- foreach( pair = spn) %dopar% {
    i <- pair[1]
    j <- pair[2]
    feature_union <- which(rowSums(features[,c(i,j)]) > 0)
    css_tmp <- css[feature_union,feature_union]
    a <-  features[feature_union,i]/sum(as.numeric(features[feature_union,i]))
    b <- features[feature_union,j]/sum(as.numeric(features[feature_union,j]))
    if (weighted == FALSE){
      a[a > 0] <- 1
      b[b > 0] <- 1
    }
    abt <- a %*% t(b)
    aat <- a %*% t(a)
    bbt <- tcrossprod(b, b)
    cssAB <- css_tmp * abt
    d <- sum(cssAB)/max(sum(css_tmp * aat), sum(css_tmp * bbt))
  }
  if (verbose) print(paste(i,j, d))
  distm[lower.tri(distm)] <- as.numeric(distlist)
  distmt <- t(distm)
  distmt[lower.tri(distmt)] <- as.numeric(distlist)
  colnames(distmt) <- sample_names
  row.names(distmt) <- sample_names

  # Make dissimilarity
  ones = matrix(data = rep(1, nrow(distmt) * ncol(distmt)),nrow = nrow(distmt), ncol = ncol(distmt))
  distmt = ones - distmt
  diag(distmt) <- 0

  return(distmt)

}
