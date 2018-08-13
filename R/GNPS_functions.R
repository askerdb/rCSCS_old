library(reshape2)
library(igraph)

#' download_GNPS
#'
#' @param ID The GNPS id to retrieve
#' @param dir Location of the output
#'
#' @return Downloads and extracts a file of the GNPS content, returns a list with paths to the buckettable and network files
#' @export
#'
#' @examples  download_GNPS("0310e20491314ddbbf12d56b592548b4", ".")
download_GNPS <- function(ID, dir){
  #TODO: handle several datasets extracted in the same dir
  url <- paste0('https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=', ID, "&view=download_cluster_buckettable")
  print(paste(ID, dir, url))
  temp_gnps <- paste0(dir, "/", ID,"_GNPS.zip")
  #f = CFILE(file = temp_gnps , mode = "wb")
  #results <- postForm(url, dummy="", binary=T, .opts = curlOptions(writedata = f@ref))
  #print(paste("Download status", results))
  #RCurl::close(f)
  system(paste0('curl -d "" \"', url, '\" -o ', temp_gnps), show.output.on.console = F)
  print(paste0('curl -d "" ', url, ' -o ', temp_gnps))
  unzip(temp_gnps, exdir = paste(dir,"/",ID,sep=""))
  files <- Sys.glob(file.path(paste(dir,"/",ID,sep=""), "*"))
  netattr <- list.files(grep("clusterinfosummarygroup_attributes_withIDs_withcomponentID", files, value = T), full.names = T)
  buckettable <- grep("download_cluster_buckettable-main.tsv", files, value = T)
  edges <- list.files(grep("networkedges_selfloop", files, value = T), full.names = T)
  ids <- list.files(paste(dir,"/",ID,"/","clusterinfosummarygroup_attributes_withIDs",sep=""), full.names = T)
  return(list(buckettable = buckettable, edges = edges, attr = netattr, ids = ids))
}



#' prepare_GNPS
#'
#' @param ID GNPS ID to prepare
#' @param dir Path to put GNPS files
#' @param select Only output features in the feature table found in the CSS matrix
#' @return A list of a feature matrix and a chemical similarity matrix (CSS)
#' @export
#'
#' @examples prepare_GNPS("0310e20491314ddbbf12d56b592548b4")
prepare_GNPS <- function(ID, dir = ".", select = T){
  gnps <- download_GNPS(ID, dir)
  features <- read.table(gnps$buckettable, sep = "\t", header=T, row.names=1, comment.char="")
  csslong <- read.table(gnps$edges, header=T)
  
  g <- igraph::graph.data.frame(csslong, directed=FALSE)
  css <- get.adjacency(g, attr="Cosine", sparse=FALSE)
  diag(css) <- 1

  return(list(features = features, css = css))
}