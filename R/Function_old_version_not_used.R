
# Note: this function is not used in the new version of Yeastspot3D
# this function is used to filter fasta file of genes which belong to metabolic
#
# @return
# @export
#
# @examples
#filterMetabolicGene <- function() {
#  geneName0 <- list.files("1011_project")
#  geneMetabolic <- paste(str_trim(gene_feature_GEM$locus_tag, side = "both"), ".fasta", sep = "")
#  geneMetabolic0 <- intersect(geneName0, geneMetabolic)
#  dir.create("target_gene")
#
#  for (i in seq_along(geneMetabolic0)) {
#    geneX <- geneMetabolic0[i]
#    file0 <- paste("1011_project/", geneX, sep = "")
#    file.copy(file0, "target_gene")
#  }
#}



#' this function is used to just preprocess the fasta file without filteration, used for the early version
#'
#' @param gene_test
#'
#' @return
#' @export
#'
#' @examples
processFasta <- function(gene_test) {
  # read each fasta file and change it into a dataframe
  #gene_test <- "YIL160C.fasta"
  gene_name_test <- str_replace_all(gene_test, ".fasta", "")
  fastaFile <- readDNAStringSet(paste("target_gene/", gene_test, sep = ""))
  # obtain the strain name information and sequence information
  seq_name <- names(fastaFile)
  sequence <- paste(fastaFile)

  # establish a dataframe contains the strain name and sequnece information
  df <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  df_list <- list()
  for (i in seq(length(df$sequence))) {
    df_list[i] <- str_split(df$sequence[i], "")
  }

  return(df_list)
}




#' this function is used to choose the gene based on strain phenotype information
#'
#' @param gene_test
#' @param strain_select
#'
#' @return
#' @export
#'
#' @examples
filterMutationStrainType <- function(gene_test, strain_select) {
  # read each fasta file and change it into a dataframe
  # then the filter can be used for each dataframe to obtain the strains we need
  # gene name
  # exampel:gene_test <- "YAL012W.fasta"
  gene_name_test <- str_replace_all(gene_test, ".fasta", "")
  fastaFile <- readDNAStringSet(paste("target_gene/", gene_test, sep = ""))
  # obtain the strain name information and sequence information
  seq_name <- names(fastaFile)
  sequence <- paste(fastaFile)

  # establish a dataframe contains the strain name and sequnece information
  df <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  strain_select['index_strain'] <- paste(strain_select$Standardized_name, "_", gene_name_test, "_", sep = "")


  for (j in seq_along(df$seq_name)) {
    exist_sign <- vector()
    for (i in seq_along(strain_select$index_strain)) {
      exist_sign[i] <- str_detect(df$seq_name[j], strain_select$index_strain[i])
    }

    if (any(exist_sign) == TRUE) {
      df$choosed[j] <- "YES"
    } else {
      df$choosed[j] <- "NO"
    }
  }

  df_refine <- filter(df, choosed == "YES")


  df_list <- list()
  for (i in seq(length(df_refine$sequence))) {
    df_list[i] <- str_split(df_refine$sequence[i], "")
  }

  dir.create("target_gene_processed")
  filename0 <- paste("target_gene_processed/",gene_name_test, ".RData", sep = "")
  save(df_list, file = filename0)
  return(df_list)
}




#' the followed two function  were used to estimate the mutation information based on input DNA fasta file
#' ---------------version 1
#' @param alted_seq
#' @param geneName
#'
#' @return
#' @export
#'
#' @examples
findPPosition0 <- function(alted_seq, geneName){
  #this function is used to find the postion of mutated amino acids based on genomics mutation
  #alted_seq <- df_list[[1]]
  #geneName <- 'YAL012W'
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature_GEM)
  gene_snp[['gene']] <- alted_seq

  #translation
  #using package seqinr
  realcds <- str_to_lower(paste(gene_snp[['gene']],collapse = ""))
  toycds <- s2c(realcds)
  gene_snp[['protein_mutated']] <- translate(seq = toycds)

  #find the relative postion of mutated amino acids
  aa_position <- which(gene_snp[['protein']] != gene_snp[['protein_mutated']] )

  #calculate the mutation number in the mutated postion (for specific strain -x)
  gene_snp[['mutation_position']] <- rep(0,length(gene_snp[['protein']])) #initialize the start value for each positions
  gene_snp[['mutation_position']][aa_position] <- 1
  result <- unlist(gene_snp[['mutation_position']])
  return(result)
}




#' function to obtain the mutation of amino acids residue in the 3d structure based on sequence blast anlysis
#' but it can be very dangeous using this method, as the insertion or deletion could lead to many mutation in a seq
#' --------------version1
#' @param geneName
#' @param mutated_gene_seq
#'
#' @return
#' @export
#'
#' @examples
countMutationProtein0 <- function (geneName, mutated_gene_seq, gene_annotation0){

  #mutated_gene_seq <- df_list
  #geneName = 'YAL012W'
  df_list <- mutated_gene_seq
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_annotation0)
  tt <- rep(0,length(gene_snp[['protein']]))
  for (i in seq(length(mutated_gene_seq))){
    if(length(gene_snp[['gene']]) != length(df_list[[i]])) {
      tt <- tt + rep(0,length(gene_snp[['protein']]))
    }
    ##to avoide the insertion or deletion in the seq
    else{
      tt <- tt + findPPosition0(df_list[[i]],geneName)
    }
  }

  return(tt)
}
