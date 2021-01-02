# -----------------Note
# This script contains obtain all the parameters related each PDB file.
# Based on the parameters obtained, the quality analysis can be conducted.


#' Obtain the resolution for one experimetal PDB file
#'
#' Get the resolution of an experimental PDB file
#'
#' @param pdbdir Directory of a pdb file
#'
#' @return Resolution for the pdb file
#' @export
#'
#' @examples
#' pdbdir0 <- "xx/data/6cp6.pdb"
#' pdb.ResolutionEX(pdbdir=pdbdir0)
pdb.ResolutionEX <- function(pdbdir) {
  # input a pdb file
  # output the resolution for the pdb file
  # experiment pdb file
  pdb <- scan(pdbdir, sep = "\n", what = "complex")
  ss1 <- which(str_detect(pdb, "REMARK   2") == TRUE)
  pdb1 <- pdb[ss1]
  ss2 <- which(str_detect(pdb1, "RESOLUTION") == TRUE)
  pdb2 <- pdb1[ss2]
  pdb3 <- str_split(pdb2, "RESOLUTION")
  # obtain the resolution
  resolution <- pdb3[[1]][2] %>%
    str_extract(.,"\\d+\\.*\\d*")
  return(resolution)
}


#' Obtain the resolution through parsing the pdb web
#'
#' Parse the pdb web https://www.rcsb.org/ to obtain the resolution for each pdb file
#' Please firstly library(rvest) when using this function
#'
#'
#' @param pdb0 A pdb id
#'
#' @return Resolution of a pdb file
#' @export
#'
#' @examples
#' pdb.ResolutionAll(pdb0=6cp6)
pdb.ResolutionAll <- function(pdb0) {
  url0 <- "https://www.rcsb.org/structure/"
  url <- paste(url0, pdb0, sep = "")
  webpage <- read_html(url)
  resolution_file <- webpage %>%
    html_nodes("li") %>%
    html_text()
  resolution_file <- resolution_file[str_detect(resolution_file, "Resolution")]
  resolution <- resolution_file[1] %>%
    str_extract(., "\\d+\\.*\\d*")
  print(resolution)
  return(as.numeric(resolution))
}


#' Get PDB residue sequence
#'
#' Extract the residue sequence from both experimental and homology pdb files
#' Please library(bio3d) firstly when using this function
#'
#'
#' @param pdbdir Directory of a pdb file
#' @param pdbid Name of the pdb file, like '6cp6.pdb'
#'
#' @return A list contains the residues sequence of each chain in a PDB file
#' @export
#'
#' @examples
#' pdbdir0 <-  "xx/data/6cp6.pdb"
#' pdb.Sequence(pdbdir=pdbdir0, pdbid='6cp6.pdb')
pdb.Sequence <- function(pdbdir, pdbid) {
  pdb <- read.pdb(pdbdir)
  atom1 <- pdb$atom
  chainAll <- unique(atom1$chain)
  seq_list <- list()
  for (j in chainAll) {
    cat("chainID")
    print(j)
    pdbid0 <- str_replace_all(pdbid, "pdb", j)
    # should remove the legend not belong to amino acids
    atom_choose <- atom1[atom1$chain == j & atom1$type == "ATOM", ]
    atom_choose_o <- select(atom_choose, resid, chain, resno)
    for (i in seq_along(atom_choose$resid)) {
      atom_choose$resid[i] <- str_to_title(atom_choose$resid[i])
      if (is.na(a(atom_choose$resid[i])) & atom_choose$resid[i] == "Fme") {
        # sometimes a molecular can not be a amino acid residue
        atom_choose_o$resid[i] <- "M"
      } else if (is.na(a(atom_choose$resid[i])) & atom_choose$resid[i] != "Fme") {
        atom_choose_o$resid[i] <- "*"
      } else {
        # change three letter into one
        atom_choose_o$resid[i] <- a(atom_choose$resid[i])
      }
    }
    atom_choose_o$not_duplicated <- !duplicated(atom_choose_o$resno)
    seq0 <- atom_choose_o[atom_choose_o$not_duplicated == TRUE, ] %>%
      select(., resid)
    seq0 <- paste(seq0$resid, collapse = "")
    print(seq0)
    seq_list[[pdbid0]] <- seq0
  }

  return(seq_list)
}



#' Modeled PDB files filtration
#'
#' Remove modeled PDB file of low quality
#'
#' @param pdb_homo A dataframe contains paramters for each homology pdb file
#' @param qmean0  A paramter to evaluate the pdb quality defined in Swiss model database
#' @param identity0 Seq-identity, identity between target protein sequence and the residues sequence of templated PDB files
#' @param similarity0 Seq_similarity, similarity between target protein sequence and the residues sequence of templated PDB files
#' @param resolution0 Resolution of the template pdb file
#'
#' @return A dataframe contain the modeled PDB files of high quality
#' @export
#'
#' @examples
#' s1 = data.frame(qmean=c(-6,1,2),Seq_identity=c(10,30,40),Seq_similarity=c(0.15,20,45), Resolution=c(6,2,3), stringsAsFactors = FALSE)
#' pdbHomoFilter(pdb_homo=s1)
pdbHomoFilter <- function(pdb_homo, qmean0=-4, identity0=25, similarity0=0.31, resolution0=3.4 ){
  colname0 <- colnames(pdb_homo)
  #check whether the four parameters as the column names
  key_para_filter <- c('qmean','Seq_identity','Seq_similarity','Resolution')
  if (all(key_para_filter %in% colname0)){
    print('all four key parameters exist')
    print('Conduct the filteration for homology PDB files')
    pdb_homo$qmean <- as.numeric(pdb_homo$qmean)
    pdb_homo$Seq_identity <- as.numeric(pdb_homo$Seq_identity)
    pdb_homo$Seq_similarity <- as.numeric(pdb_homo$Seq_similarity)
    pdb_homo$Resolution <- as.numeric(pdb_homo$Resolution)
    pdb_homo_filter <- filter(pdb_homo, pdb_homo$qmean >= qmean0 &
                                pdb_homo$Seq_identity >= identity0 &
                                pdb_homo$Seq_similarity >=similarity0 &
                                pdb_homo$Resolution <= resolution0 )
  } else{
    print('Please prepare the detailed parameters for the quality analysis')
  }

  # this function missing return() 

}



#'  Experimental PDB files filtration
#'
#'  Remove experimental PDB files of low quality
#'
#' @param pdb_EX A dataframe contains paramters for each experimental pdb file
#' @param pident0 Identity between target protein sequence and the residues sequence included in experimental PDB files
#' @param mismatch0 Mismatch between target protein sequence and the residues sequence included in experimental PDB files
#' @param resolution0 Resolution of PBB file
#'
#' @return A dataframe contain the experimental PDB files of high quality
#' @export
#'
#' @examples
#' s1 = data.frame(pident=c(99,98,100),mismatch=c(2,3,0), Resolution=c(6,2,3), stringsAsFactors = FALSE)
#' pdbExFilter(pdb_EX=s1)
pdbExFilter <- function(pdb_EX, pident0=100, mismatch0=0, resolution0=3.4 ){
  colname0 <- colnames(pdb_EX)
  #check whether the four parameters as the column names
  key_para_filter <- c('pident','mismatch','Resolution')
  if (all(key_para_filter %in% colname0)){
    print('all three key parameters exist')
    print('Conduct the filteration for experimental PDB files')
    pdb_EX$pident <- as.numeric(pdb_EX$pident)
    pdb_EX$mismatch <- as.numeric(pdb_EX$mismatch)
    pdb_EX$Resolution <- as.numeric(pdb_EX$Resolution)
    pdb_EX_filter <- filter(pdb_EX, pdb_EX$pident >= pident0 &
                              pdb_EX$mismatch <= mismatch0 &
                              pdb_EX$Resolution <= resolution0 )
  } else{
    print('Please prepare the detailed parameters for the quality analysis')
  }
  
  # this function missing return() 
}




#' Calculate the distance between each two alpha carbons
#'
#'
#' @param pdbdir Dir for a pdb file
#' @param chainID ChainID for a chain
#'
#' @return A matrix contains the distance between each two alpha carbons
#' @export
#'
#' @examples
#' library(bio3d)
#' library(seqinr)
#' infile <- "xx/data/"
#' pdbid <- '6cp6.pdb'
#' pdbdir <- paste(infile, pdbid, sep = "")
#' pdb.ResidueDistance(pdbdir, chainID = 'K')
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
