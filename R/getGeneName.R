#' function get the gene name based on the chr name and mutation position, which can be used for the SNP annotation
#'
#' @param chr
#' @param mutated_positions
#' @param gene_annotation
#'
#' @return
#' @export
#'
#' @examples
getGeneName <- function(chr,mutated_positions,gene_annotation = gene_feature0){
  #input:
  #1. chr: chromsome
  #2. mutated_positiion
  #3. gene_featured0: contains the gene sequence information from chromsome of sec-s288c ,like the start and end
  #output:
  # the gene name contained this mutation
  ss <- filter(gene_feature0,
               chromosome == chr &
                 start <= mutated_positions &
                 end >= mutated_positions)
  if(length(ss$locus_tag)){
    ss0 <- ss$locus_tag
  } else{
    ss0 <- "INTERGENIC"
  }
  return(ss0)
}
