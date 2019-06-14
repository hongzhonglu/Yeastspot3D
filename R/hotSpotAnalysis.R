#' Mutation hotspot analysis
#'
#' @param gene0 A gene systematic name
#' @param pdbID A string
#' @param SNPlist0 A SNP list for the strains from specific phenotype
#' @param gene_annotation0 The gene annotation summary
#' @param pdb The dir store the residue distance matrix, should be a txt file seperated by ',' or the residue distance matrix
#' @param sstart0 The start residue coordinate for the resdiues in the original protein
#' @param send0 The end residue coordinate for the residues in the the original protein
#' @param qstart0 The start residue coordinate for the resdiues in the PDB file
#' @param qend0 The end residue coordinate for the residues in the PDB file
#' @param result_dir The directory to save the hot spot analysis result
#' @param strain_type A string to represent the source of SNP
#' @param input_dir
#'
#'
#' @return
#' @export important_hot analysis result
#'
#' @examples
hotSpotAnalysis <- function (gene0 = ss0,
                             pdbID = NA,
                             SNPlist0 = mutated_gene1,
                             gene_annotation0 = gene_feature0,
                             pdb = distance,
                             sstart0 = p1,
                             send0 = p2,
                             qstart0 = q1,
                             qend0 = q2,
                             result_dir = outfile0,
                             strain_type = NA,
                             input_dir=TRUE) {

  # step 1
  # preprocess the SNP information
  gene_snp <- getGeneCoordinate(gene_name = gene0, genesum = gene_annotation0)
  gene_snp[["pro_mutation_count"]] <- countMutationProtein(gene_name = gene0, mutation_annotation = SNPlist0, gene_snp0 = gene_snp)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  pos_mutation <- which(gene_snp[["pro_mutation_count"]] != 0)


  # step 2 input the structure information
  # input the distance of all the pired residues
  r3 <- paste(qstart0, qend0, sep = "-")

  # the residue distance can be a directory or matrix
  if(input_dir){
    ResidueDistance0 <- read.table(pdb, sep = ",") # in the followed calculation, the matrix dosen't have the col and row names
  } else{
    ResidueDistance0 <- pdb
  }


  ResidueDistance0 <- as.matrix(ResidueDistance0)
  ResidueDistance <- ResidueDistance0 # [r1:r2,r1:r2]


  # the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
  # obtain the mutation information for the structure
  seq_3D_origin <- sstart0:send0 # seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
  p3 <- paste(sstart0,send0, sep = "-")
  amino_acid_3D <- gene_snp[["protein"]][seq_3D_origin]
  count_mutation_3D <- gene_snp[["pro_mutation_count"]][seq_3D_origin]

  # mutation position on structure and #mutation number on structure
  pos_mutation_3D <- which(count_mutation_3D != 0)
  seq_3D <- 1:length(count_mutation_3D) # seq0 is the coordinate of PDB structure
  mutation_count_3D <- count_mutation_3D[pos_mutation_3D]

  # there should be two postions which have mutations
  if (length(pos_mutation_3D) >= 2) {
    # wap calculation for each pair mutated residue
    # calculate the standardard sample number
    sample_standard1 <- sampleStand(count_mutation_3D)

    # step 3 hot spot analysis
    # this main function will be divided into different parts for easy understanding
    pos_residue1 <- list()
    for (i in seq_along(SNPlist0$Chr)) {
      pos_residue1[[i]] <- PositionResidueSNP(SNPlist0$Pos[i], SNPlist0$Alt[i], gene0, gene_feature = gene_annotation0)
    }

    pos_residue_df <- ResidueSum(pos_residue1)


    # mapping the mutate residue onto the original protein sequence
    gene_snp[["residue"]] <- getMultipleMatchParameter(pos_residue_df$residue, pos_residue_df$pos, gene_snp[["pro_coordinate"]])
    residue_3D <- gene_snp[["residue"]][seq_3D_origin]


    # obtain the paired residue
    # if the residue position is samller than 2, the folloed function will raise error and should be passed

    residue_pair <- getHotVertice(aa_3d = seq_3D, residue0 = residue_3D, aa_pro = seq_3D_origin, distance0 = ResidueDistance)
    # remove the *@@54
    residue_pair <- removeStopCoden(residue_pair)

    if (length(residue_pair$V1) >= 1) {
      # calculate closeness of each cluster
      important_hot <- clusterAnalysis(residue_pair)


      # further calculate the p value of each choosed cluster
      # calculate the pvalue for each clust
      important_hot$pvalue <- getHotPvalue(
        cluster0 = important_hot$cluster,
        sample_standard = sample_standard1,
        distance = ResidueDistance,
        seq = seq_3D
      )

      # add sample information for the result
      important_hot$gene <- gene0
      important_hot$seq_3D_origin <- p3
      important_hot$structure <- pdbID
      important_hot$seq_3D <- r3
      important_hot$strain_type <- strain_type
      outfile <- paste(result_dir, "/", pdbID, "_", gene0, ".txt", sep = "")
      # last step: get the mutate residue coordinate from protein seqence
      # coordinate mapping
      coordinate_mapping <- mappingCoordinateFrom3DtoProtein(aa_3d = seq_3D, residue0 = residue_3D)
      important_hot$cluster <- getOriginalCoordinateProtein(coordinate0 = important_hot$cluster, coordinate_mapping0 = coordinate_mapping)


      write.table(important_hot, outfile, row.names = FALSE, sep = "\t")
    } else {
      print("------NO sigificant pairs")
    }
  } else {
    print("------Not enough mutation")
  }

}
