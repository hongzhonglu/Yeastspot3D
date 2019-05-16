# This script contains obtain all the parameters related each PDB file.
# Based on the parameters obtained, the quality analysis can be conducted.


# this code is used to calculate the distance between each two alpha carbons
#' Title
#'
#' @param pdbdir dir for a pdb file
#' @param chainID chainID for a chain
#'
#' @return
#' @export k a matrix
#'
#' @examples
pdb.ResidueDistance <- function (pdbdir,chainID) {


  pdb <- read.pdb(pdbdir)
  sele <- atom.select(pdb, type = "ATOM", "calpha", chain = chainID)
  pdb0 <- trim.pdb(pdb, sele)
  newPDB <- pdb0$atom
  k <- dm(pdb0, inds = "calpha")
  # for the element on the symmetry of matrix, we set it at 0
  diag(k) <- 0
  # obtain the element on the lower triangle
  ss <- dim(k)
  for (i in 1:ss[1]){
    for (j in 1:i){
      if(i ==j){
        k[i,j] <- 0
      } else{
        k[i,j] <- k[j,i]
      }
    }
  }

  return(k)
}

# example
# library(bio3d)
# library(seqinr)
# infile <- "data/"
# pdbid <- '6cp6.pdb'
# pdbdir <- paste(infile, pdbid, sep = "")
# pdb.ResidueDistance(pdbdir, chainID = 'K')
