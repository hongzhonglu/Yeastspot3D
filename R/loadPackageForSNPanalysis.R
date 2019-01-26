#' This function is used to load the package used for SNP analysis
#'
#' @return
#' @export
#'
#' @examples
loadPackageForSNPanalysis <- function() {
  # During the analysis we need the followed package, if your computer don't have them, please install it firstly before
  # the analysis
  library(Yeastspot3D)
  library(tidyverse)
  library(stringr)
  library(readxl)
  library(Biostrings)
  library(filesstrings) # move the files
  library(hongR)
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("Biostrings")
  library(centiserve) # this package is used to calculate the closeness centrality
  library(igraph) # form the unique clust based on Floyd-Warshall shortest-paths algorithm
  library(seqinr)
  library(readr)
}
