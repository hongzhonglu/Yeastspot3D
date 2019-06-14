#' Mutation enrichment analysis
#'
#' Conduct the mutation mapping analysis based on protein 3D structures
#'
#' @param gene0 A gene systematic name
#' @param SNPlist0 A SNP list for the strains from specific phenotype
#' @param gene_annotation0 The gene annotation summary
#' @param pdb The dir store the residue distance matrix, should be a txt file seperated by ',' or it is a residue distance matrix
#' @param sstart0 The start residue coordinate for the resdiues in the PDB file
#' @param send0 The end residue coordinate for the residues in the PDB file
#' @param input_dir A logical vector
#'
#' @return p_value
#' @export
#'
#' @examples
#' # Load the gene feature gene and SNP dataset
#' data('gene_feature0')
#' data('snp_data')
#' # Annotate the SNP
#' mutated_gene <- annotateSNP(snp_input = snp_data, gene_feature = gene_feature0)
#' # Load the residue distance matrix
#' data('ResidueDistance_YPR184W')
#' # Get the SNP list for specific gene
#' mutated_gene1 <- filter(mutated_gene, Gene2 == 'YPR184W')
#' result0 <- clumpsAnalysis(gene0 = 'YPR184W',
#'                          SNPlist0 = mutated_gene1,
#'                          gene_annotation0 = gene_feature0,
#'                          pdb = ResidueDistance_YPR184W,
#'                          sstart0 = 2,
#'                          send0 = 1534,
#'                          input_dir= FALSE)
clumpsAnalysis <- function(gene0, SNPlist0, gene_annotation0, pdb, sstart0, send0, input_dir=TRUE) {
  # step 1
  # preprocess the SNP information
  gene_snp <- getGeneCoordinate(gene_name = gene0, genesum = gene_annotation0)
  gene_snp[["pro_mutation_count"]] <- countMutationProtein(gene_name = gene0, mutation_annotation = SNPlist0, gene_snp0=gene_snp)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  pos_mutation <- which(gene_snp[["pro_mutation_count"]] != 0)


  # step 2 input the structure information
  # input the distance of all the pired residues
  if(input_dir){
    ResidueDistance0 <- read.table(pdb, sep = ",") # in the followed calculation, the matrix dosen't have the col and row names
  } else{
    ResidueDistance0 <- pdb
  }

  ResidueDistance0 <- as.matrix(ResidueDistance0)
  ResidueDistance <- ResidueDistance0 # [r1:r2,r1:r2]


  # the amino acid sequence in structure is from 2:394 while  the original sequence is from 1:394
  # obtain the mutation information for the structure
  p1 <- sstart0
  p2 <- send0
  p3 <- paste(p1, p2, sep = "-")
  seq_3D_origin <- p1:p2 # seq_from_3D <- 2:394 #"YAL012W.fasta"#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
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

    # step 3
    # calculate the standardard sample number
    sample_standard1 <- sampleStand(count_mutation_3D)

    # calculate the wap for each pair of mutated residues based on mutation postion
    wap_original <- getTotalWAP(pos_mutation_3D, sample_standard1, ResidueDistance)

    # change the postion of mutation while keep the mutation number in each postion
    # only change the postion but not change the mutated number???
    wap_sample0 <- getSampleWAP(pos_mutation_3D, sample_standard1, ResidueDistance, seq = seq_3D, n = 10000)

    # analyze the result
    plotNullDistribution(wap_sample0)
    p_value <- getPvalue(wap_original, wap_sample0)
    print(paste("-------p_value=", p_value, sep = ""))
  } else {
    p_value <- "NA"
    print("------Not enough mutation")
  }

  return(p_value)
}
