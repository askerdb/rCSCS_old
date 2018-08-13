library(foreach)

#' Calculate CSCS distance
#'
#' @param features Matrix of feature intensities in each sample (pxn)
#' @param css Square matrix of cosine similarities of features (pxp)
#' @param dissimilarity Output dissimilarity matrix
#' @param cosine_threshold Only include features below this threshold
#' @param weighted Weight features by intensity (TRUE) or absence/presence (FALSE)
#' @param verbose Print extra details
#'
#' @details The value of cosine_threshold is 0.6 in Sedio et al. but for certain applications other values might be better.
#' @references Sedio et al, 2016
#' @return A pxp matrix of CSCS distances
#' @export
#'
#' @examples
#'
#' #Get arbitrary GNPS project
#' gnps <- prepare_GNPS("0310e20491314ddbbf12d56b592548b4")
#' dist <- distance_cscs(gnps$features, gnps$css)
#'
#' #GlobalEuphorbia data
#' distance_cscs(GEfeatures, GEcss)
distance_cscs <- function(features, css, dissimilarity = F, cosine_threshold = 0.6, weighted = T, verbose = F){
  stopifnot(nrow(features) == nrow(css))
  diag(css) <- 1
  css[css < cosine_threshold] <- 0
  css <- as.matrix(css)
  sample_names <- colnames(features)
  distm <- matrix(0, nrow = length(sample_names), ncol = length(sample_names),
                  dimnames = list(sample_names, sample_names))
  spn <- combn(sample_names, 2, simplify=FALSE)
  distlist <- foreach::foreach( pair = spn) %dopar% {
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
  if (dissimilarity == T){
  distmt <- 1 - distmt
  }

  diag(distmt) <- 1
  return(as.dist(distmt))

}

distance_cscs_mh <- function(features, css, dissimilarity = F, cosine_threshold = 0.6, weighted = T, verbose = F){
  stopifnot(nrow(features) == nrow(css))
  diag(css) <- 1
  css[css < cosine_threshold] <- 0
  css <- as.matrix(css)
  sample_names <- colnames(features)
  distm <- matrix(0, nrow = length(sample_names), ncol = length(sample_names),
                  dimnames = list(sample_names, sample_names))
  spn <- combn(sample_names, 2, simplify=FALSE)
  distlist <- foreach::foreach( pair = spn) %dopar% {
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
  diag(distmt) <- 1
  # Make dissimilarity
  if (dissimilarity == T){
    distmt <- 1 - distmt
  }


  return(as.dist(distmt))

}
