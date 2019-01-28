#' function used to give summary of SNP for gene based on a specific strain phenotype
#'
#' @param gene0 A gene systematic name
#' @param SNPlist0 A SNP list for the strains from specific phenotype
#' @param gene_annotation0 The gene annotation summary
#' @param pdbID0 A pdb id which has been mapped onto the gene name (the first parameter)
#' @param sstart0 The start residue coordinate for the resdiues in the PDB file
#' @param send0 The end residue coordinate for the residues in the PDB file
#'
#' @return
#' @export CLUMPS analysis result
#'
#' @examples
printSNPforGene <- function(gene0 = ss0,
                            SNPlist0 = mutated_gene1,
                            gene_annotation0 = gene_feature_GEM,
                            pdbID0 = pdbID,
                            sstart0 = p1,
                            send0 = p2) {
  # input the gene name
  ss <- gene0
  # input the protein coordinate ID
  p3 <- paste(sstart0, send0, sep = "-")
  seq_3D_origin <- sstart0:send0 # this is the coordinated of original protein sequence and should changed into 3D structure coordinates

  # produce the mutaed residue and its position
  pos_residue1 <- list()
  for (j in 1:nrow(SNPlist0)) {
    cat('Process all the SNPs from SNP list to obtain the SNP belong to the input gene:')
    print(j)
    pos_residue1[[j]] <- PositionResidueSNP(SNPlist0$Pos[j], SNPlist0$Alt[j], ss)
  }

  pos_residue_df <- ResidueSum(pos_residue1)

  # mapping the mutate residue onto the original protein sequence
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_annotation0)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  gene_snp[["residue"]] <- getMultipleReactionFormula(pos_residue_df$residue, pos_residue_df$pos, gene_snp[["pro_coordinate"]])
  residue_3D <- gene_snp[["residue"]][seq_3D_origin]

  # analyse the mutation frequence of each residue
  residue_3D0 <- residue_3D[!is.na(residue_3D)]
  residue_3D0 <- str_split(residue_3D0, ";")
  residue_3D0 <- unlist(residue_3D0)
  unique(residue_3D0)
  tmp <- table(residue_3D0)
  result <- as.data.frame(tmp)
  result <- result %>% separate(., residue_3D0, into = c("alt", "position"), sep = "@@")
  result$position <- as.numeric(result$position)
  result <- result %>% arrange(., position)
  # obtain the original ref residue
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_annotation0)
  # output the result
  result$ref <- getSingleReactionFormula(gene_snp[["protein"]], gene_snp[["protein_coordinate"]], result$position)
  result$orf <- ss
  result$pdbID <- pdbID0
  result <- select(result, orf, ref, position, alt, Freq, pdbID)
  cat('The SNP summary:')
  print(result)
  # save the results in format which can predict the protein function
  result0 <- select(result, orf, ref, position, alt)
  result0 <- result0 %>% unite(mutation, c("ref", "position", "alt"), sep = "")
  # merge the combined information with the result, the new column could be used for the annotation of SNP based on mutfunc database
  result$mutation <- result0$mutation
  return(result)
}
