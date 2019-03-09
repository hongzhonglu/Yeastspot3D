# This script contains obtain all the parameters related each PDB file.
# Based on the parameters obtained, the quality analysis can be conducted.


#' homology pdb file filter function
#'
#' @param pdb_homo a dataframe contains paramters for each homology pdb file
#' @param qmean0  paramter to evaluate the pdb quality defined in Swiss model database
#' @param identity0 Seq-identity, identity between target protein sequence and the residues sequence of templated PDB files
#' @param similarity0 Seq_similarity, similarity between target protein sequence and the residues sequence of templated PDB files
#' @param resolution0 Resolution
#'
#' @return
#' @export
#'
#' @examples
pdbHomoFilter <- function(pdb_homo, qmean0=-4, identity0=25, similarity0=0.31, resolution0=3.4 ){
  colname0 <- colnames(pdb_homo)
  #check whether the four parameters as the column names
  key_para_filter <- c('qmean','Seq_Identity','Seq_similarity','Resolution')
  if (all(key_para_filter %in% colname0)){
    print('all four key parameters exist')
    print('Conduct the filteration for homology PDB files')
    pdb_homo$qmean <- as.numeric(pdb_homo$qmean)
    pdb_homo$Seq_Identity <- as.numeric(pdb_homo$Seq_Identity)
    pdb_homo$Seq_similarity <- as.numeric(pdb_homo$Seq_similarity)
    pdb_homo$Resolution <- as.numeric(pdb_homo$Resolution)
    pdb_homo_filter <- filter(pdb_homo, pdb_homo$qmean >= qmean0 &
                                pdb_homo$Seq_Identity >= identity0 &
                                pdb_homo$Seq_similarity >=similarity0 &
                                pdb_homo$Resolution <= resolution0 )
  } else{
    print('Please prepare the detailed parameters for the quality analysis')
  }
}



#' experimental pdb file filter function
#'
#' @param pdb_EX a dataframe contains paramters for each experimental pdb file
#' @param pident0 identity between target protein sequence and the residues sequence included in experimental PDB files
#' @param mismatch0 mismatch between target protein sequence and the residues sequence included in experimental PDB files
#' @param resolution0
#'
#' @return
#' @export
#'
#' @examples
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
}
