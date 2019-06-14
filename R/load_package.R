#' Load essential packages
#'
#' Load the package used for SNP analysis
#'
#' @return Display whether all essential packages are installed.
#' @export
#'
#' @examples
#' library(Yeastspot3D)
#' loadPackageForSNPanalysis()
loadPackageForSNPanalysis <- function() {
  # During the analysis we need the followed package, if your computer don't have them, please install it firstly before
  # the analysis
    essential_package <- c('Yeastspot3D','tidyverse','stringr','readxl','Biostrings','filesstrings','centiserve','igraph','readr')

  if (!all(essential_package %in% rownames(installed.packages()))){
    print('Please firstly install the followed packages')
    print(essential_package)} else{
      print('All essential packages have been installed')
      print('Load essential packages')
      library(Yeastspot3D)
      library(tidyverse)
      library(stringr)
      library(readxl)
      library(Biostrings)
      library(filesstrings) # move the files
      #library(hongR) # have moved three common functions from hongR into Yeastspot3D, so this package was not used.
      #source("https://bioconductor.org/biocLite.R")
      #biocLite("Biostrings")
      library(centiserve) # this package is used to calculate the closeness centrality
      library(igraph) # form the unique clust based on Floyd-Warshall shortest-paths algorithm
      library(seqinr)
      library(readr)
    }
}


