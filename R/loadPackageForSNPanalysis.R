#' This function is used to load the package used for SNP analysis
#'
#' @return
#' @export
#'
#' @examples
loadPackageForSNPanalysis <- function() {
  # During the analysis we need the followed package, if your computer don't have them, please install it firstly before
  # the analysis
  print('Please firstly install the followed packages')
  essential_package <- c('Yeastspot3D','tidyverse','stringr','readxl','Biostrings','filesstrings','centiserve','igraph','readr')
  print(essential_package)
  #library(Yeastspot3D)
  #library(tidyverse)
  #library(stringr)
  #library(readxl)
  #library(Biostrings)
  #library(filesstrings) # move the files
  #library(hongR) # have moved three common functions from hongR into Yeastspot3D, so this package was not used.
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("Biostrings")
  #library(centiserve) # this package is used to calculate the closeness centrality
  #library(igraph) # form the unique clust based on Floyd-Warshall shortest-paths algorithm
  #library(seqinr)
  #library(readr)
  for (i in essential_package){
    checkInstallation(i)
  }
}



checkInstallation <- function(package0){
  if(package0 %in% rownames(installed.packages())){
    print(paste(package0,'is installed'))
    library(package0)
  } else{
    print(paste('Please install',package0))
  }
}