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



#' this function was used to parse the genome annotation of s288c
#' the result of this script is the base for all other analysis and should be run firstly
#' firstly we need download gene annotation information from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64
#' which include:
#' 1) GCF_000146045.2_R64_genomic.gbff
#' 2) GCF_000146045.2_R64_cds_from_genomic.fna
#' secondly we download the gene annotation information from SGD
#' 3) s288_genome, which include the colnames:  [1] "DBID"               "systematic_name"    "gene_standard_name" "gene_name"
#' "locations_strand"   "locations_start"    "locations_end"      "chromosome"
#' "sequence"           "sequnece_length"    "protein_residue"    "protein_length"
#'
#' @return
#' @export
#'
#' @examples
mergeGeneAnnotationFeature <- function() {
  dir1 <- "data/GCF_000146045.2_R64_genomic.gbff"
  dir2 <- "data/GCF_000146045.2_R64_cds_from_genomic.fna"
  dir3 <- "data/s288_genome.tsv"
  ## parse the genomic gbff file
  cat("parse the genome annotation from GCF_000146045.2_R64_genomic.gbff")
  s288 <- scan(dir1, sep = "\n", what = "complex")
  s288_n <- s288[72:length(s288)]
  s288_n <- str_replace_all(s288_n, "     ", "")

  # gene feature summary
  # get the index for each gene
  gene <- which(str_detect(s288_n, "gene  ") == TRUE)
  gene_name <- gene + 1
  gene_orf <- gene + 2

  # produce the dataframe
  gene_annotation <- data.frame(gene_name = character(length = length(gene_name)), stringsAsFactors = FALSE)
  gene_annotation$gene_name <- s288_n[gene_name]
  gene_annotation$gene_orf <- s288_n[gene_orf]
  gene_annotation$complement <- s288_n[gene]
  gene_annotation$gene_name <- str_replace_all(gene_annotation$gene_name, "/", "")
  gene_annotation$gene_orf <- str_replace_all(gene_annotation$gene_orf, "/", "")

  for (i in seq(length(gene_annotation$gene_name))) { # seq in R and range in python
    if (str_detect(gene_annotation$gene_name[i], "locus_tag")) {
      gene_annotation$gene_orf[i] <- gene_annotation$gene_name[i]
      gene_annotation$gene_name[i] <- ""
    } else {
      gene_annotation$gene_orf[i] <- gene_annotation$gene_orf[i]
      gene_annotation$gene_name[i] <- gene_annotation$gene_name[i]
    }
  }

  gene_annotation0 <- select(gene_annotation, gene_orf, complement)
  colnames(gene_annotation0) <- c("locus_tag", "location")

  gene_annotation0$locus_tag <- str_replace_all(gene_annotation0$locus_tag, "locus_tag=", "") %>%
    str_replace_all(., "\"", "") %>%
    str_trim(., side = "both")

  gene_annotation0$location <- str_replace_all(gene_annotation0$location, "gene ", "")



  ## cds fna file analysis
  cat("parse the gene annotation from cds zone")
  s288_cds <- scan(dir2, sep = "\n", what = "complex")
  s288_cds[1:200]
  index1 <- which(str_detect(s288_cds, ">"))

  cds <- s288_cds[index1]
  cds <- str_replace_all(cds, ">", "")
  cds <- str_replace_all(cds, " \\[", "@")
  cds <- str_replace_all(cds, "\\]", "")
  cds0 <- str_split(cds, "@")

  locus <- vector()
  gene <- vector()
  location_index <- vector()
  location <- vector()

  for (i in seq(length(cds0))) {
    locus[i] <- which(str_detect(cds0[[i]], "locus_tag="))
    gene[i] <- cds0[[i]][locus[i]]
    location_index[i] <- which(str_detect(cds0[[i]], "location="))
    location[i] <- cds0[[i]][location_index[i]]
  }

  cds_fna <- data.frame(gene = gene, location = location, stringsAsFactors = FALSE)


  seq_cds <- list()
  for (i in seq(length(index1) - 1)) {
    seq_cds[[i]] <- s288_cds[(index1[i] + 1):(index1[i + 1] - 1)]
  }
  seq_cds[[6008]] <- s288_cds[(index1[length(index1)] + 1):length(s288_cds)]

  nchar(paste(seq_cds[[6008]], sep = "", collapse = ""))


  for (i in 1:6008) {
    cds_fna$cds[i] <- paste(seq_cds[[i]], sep = "", collapse = "")
    cds_fna$length_cds[i] <- nchar(cds_fna$cds[i])
  }

  cds_fna$gene <- str_replace_all(cds_fna$gene, "locus_tag=", "")


  # give the mRNA seq, mRNA length, amino acid sequece and amino acid length
  # here we merge the gene info from gene_annotation0, sds_fna info and SGD info
  cat("input the gene annotation from sgd database")
  s288_SGD <- read_tsv(dir3)
  gene_annotation0$cds_location <- getSingleMatchParameter(cds_fna$location, cds_fna$gene, gene_annotation0$locus_tag)
  gene_annotation0$cds_seq <- getSingleMatchParameter(cds_fna$cds, cds_fna$gene, gene_annotation0$locus_tag)
  gene_annotation0$cds_length <- getSingleMatchParameter(cds_fna$length_cds, cds_fna$gene, gene_annotation0$locus_tag)
  gene_annotation0 <- gene_annotation0[gene_annotation0$cds_location != "NA", ] # total 6008 genes coud can be translated into proteins

  gene_annotation0$aa_seq <- getSingleMatchParameter(s288_SGD$protein_residue, s288_SGD$systematic_name, gene_annotation0$locus_tag)
  gene_annotation0$aa_length <- getSingleMatchParameter(s288_SGD$protein_length, s288_SGD$systematic_name, gene_annotation0$locus_tag)
  gene_annotation0$chromosome <- getSingleMatchParameter(s288_SGD$chromosome, s288_SGD$systematic_name, gene_annotation0$locus_tag)
  gene_annotation0$start <- getSingleMatchParameter(s288_SGD$locations_start, s288_SGD$systematic_name, gene_annotation0$locus_tag)
  gene_annotation0$end <- getSingleMatchParameter(s288_SGD$locations_end, s288_SGD$systematic_name, gene_annotation0$locus_tag)
  # gene_annotation0$DNA_SGD <- getSingleMatchParameter(s288_SGD$sequence,s288_SGD$systematic_name,gene_annotation0$locus_tag) It should be note the gene sequence from SGD seems not right
  gene_annotation0$start <- as.numeric(gene_annotation0$start)
  gene_annotation0$end <- as.numeric(gene_annotation0$end)
  gene_annotation0 <- gene_annotation0[gene_annotation0$chromosome != "NA", ]
  gene_annotation0$complement_sign <- str_detect(gene_annotation0$cds_location, "complement")
  cat("The genome annotation summary is done!")
  return(gene_annotation0)
}




#' this functio was used to evaluate the quality of gene annotation from cds to protein sequence
#' three check step
#' check 1, ratio of gene and the translated protein----OK
#' check 2, length of gene----OK
#' check 3, translated protein seq is equal to that in SGD, index_need_check 570 1220 1221 1222 1223 1224 1225 1226
#'
#' @param gene_feature_original
#'
#' @return
#' @export
#'
#' @examples
qualityCheckFromCDStoProtein <- function(gene_feature_original) {
  print("start check the translation")
  translatedProtein <- checkTanslatedProtein(gene_feature_original)
  print("the whole genome coordinate information have been checked")
  gene_feature_need_check2 <- translatedProtein[[2]]
  for (i in seq_along(gene_feature_need_check2$locus_tag)) {
    gene_feature_need_check2$aa_seq[i] <- getProteinFromCDS(gene_feature_need_check2$cds_seq[i])
  }
  gene_feature_right2 <- translatedProtein[[1]]
  gene_feature_original <- rbind.data.frame(gene_feature_right2, gene_feature_need_check2)
  #check the relation between the protein length and the gene length
  gene_feature_original$check <- ((as.numeric(gene_feature_original$cds_length))/3-1) == as.numeric(gene_feature_original$aa_length)
  if(length(which(gene_feature_original$check==FALSE))>=1){
    print('The lenght of some proteins is not consistent with the related gene lenght')
  } else{
    print('The lenght of proteins is consistent with the related gene lenght')
  }
  return(gene_feature_original)
}


#' Function to parse the UniProt annotation with only single site(?)
#'
#' @param site_dataframe a dataframe to store the active site annotation from UniProt database
#' @param site_type0 the site type classification by UniProt database
#'
#' @return
#' @export
#'
#' @examples
parseSingleSite <- function(site_dataframe, site_type0){
  site_type <- site_type0
  colnames(site_dataframe) <- c('Entry',site_type)
  active_site <- site_dataframe
  active_site[[site_type]] <- str_replace_all(active_site[[site_type]], "\\{.*?\\}", "") %>%
    str_replace_all(.,"\\. \\.","") %>%
    str_replace_all(.,"\\.;",";")


  active_site1 <- splitAndCombine(active_site[[site_type]], active_site$Entry, sep0 = ";")

  ss <- active_site1$v1
  ss1 <- str_split(ss," ")
  singleSite <- vector()
  for (i in seq_along(ss1)) {
    if (site_type %in% ss1[[i]]) {
      ss <- which(ss1[[i]] == site_type)
      singleSite[i] <- ss1[[i]][ss + 1]
    } else {
      singleSite[i] <- "NA"
    }
  }

  active_site1$v3 <- singleSite
  active_site1$v4 <- site_type

  return(active_site1)
}




#' Function to parse the UniProt annotation with only single or several site(?)
#'
#' @param site_dataframe
#' @param site_type0
#' @param single_site
#'
#' @return
#' @export
#'
#' @examples
parseMutipleSite <- function(site_dataframe, site_type0, single_site = TRUE) {
  site_type <- site_type0
  colnames(site_dataframe) <- c("Entry", site_type)
  active_site <- site_dataframe
  active_site[[site_type]] <- str_replace_all(active_site[[site_type]], "\\{.*?\\}", "") %>%
    str_replace_all(., "\\. \\.", "") %>%
    str_replace_all(., "\\.;", ";")


  active_site1 <- splitAndCombine(active_site[[site_type]], active_site$Entry, sep0 = ";")

  ss <- active_site1$v1
  ss1 <- str_split(ss, " ")
  singleSite <- vector()
  for (i in seq_along(ss1)) {
    if (site_type %in% ss1[[i]]) {
      ss <- which(ss1[[i]] == site_type)
      if (single_site) {
        singleSite[i] <- ss1[[i]][ss + 1]
      } else {
        singleSite[i] <- paste(ss1[[i]][ss + 1], ss1[[i]][ss + 2], sep = "-")
      }
    } else {
      singleSite[i] <- "NA"
    }
  }

  active_site1$v3 <- singleSite
  active_site1$v4 <- site_type

  return(active_site1)
}



#' Function to summarize all the active site of potein for S. cerevisiae s288c based on
#' protein function annotation from uniprot database
#' @return
#' @export
#'
#' @examples
mergeActiveSiteInf <- function() {
  # download the sce-active site information from uniprot database for S. cerevisiae s288c
  # also prepare the id mapping between the id from uniprot and the gene systematic name

  dir1 <- "data/sce_active site.xlsx"
  dir2 <- "data/uniprotGeneID_mapping.xlsx"
  # input the data
  sce_site00 <- read_excel(dir1)

  # parse the active site
  active_site0 <- select(sce_site00, Entry, Active_site) %>%
    filter(., !is.na(Active_site))

  active_site1 <- parseSingleSite(active_site0, site_type0 = "ACT_SITE")

  # parse the metal binding
  metal_site <- select(sce_site00, Entry, Metal_binding) %>%
    filter(., !is.na(Metal_binding))

  metal_site1 <- parseSingleSite(metal_site, site_type0 = "METAL")

  # parse the binding site
  binding_site <- select(sce_site00, Entry, Binding_site) %>%
    filter(., !is.na(Binding_site))

  binding_site1 <- parseSingleSite(binding_site, site_type0 = "BINDING")

  # parse the modified residue
  modified_site <- select(sce_site00, Entry, Modified_residue) %>%
    filter(., !is.na(Modified_residue))

  modified_site1 <- parseSingleSite(modified_site, site_type0 = "MOD_RES")


  # parse the PTM
  # will do based on the need
  PTM_site <- select(sce_site00, Entry, Post_translational_modification) %>%
    filter(., !is.na(Post_translational_modification))


  # parse the DNA binding
  Nucleotide_binding_site <- select(sce_site00, Entry, Nucleotide_binding) %>%
    filter(., !is.na(Nucleotide_binding))

  Nucleotide_binding_site1 <- parseMutipleSite(Nucleotide_binding_site, site_type0 = "NP_BIND", single_site = FALSE)


  # merge the above information
  sce_site_refine <- rbind.data.frame(active_site1, metal_site1, binding_site1, modified_site1, Nucleotide_binding_site1, Nucleotide_binding_site1)
  colnames(sce_site_refine) <- c("description", "Entry", "coordinate", "site_type")
  sce_site_refine <- filter(sce_site_refine, coordinate != "NA")

  # find the gene orf name
  ID_mapping <- read_excel(dir2)
  sce_site_refine$orf <- getMultipleMatchParameter(ID_mapping$GeneName, ID_mapping$Entry, sce_site_refine$Entry)
  which(str_detect(sce_site_refine$orf, ";"))

  # refine the dataformat
  # one entry id from uniprot database could mapping onto gene orf ID
  ss <- sce_site_refine %>% unite(., inf, c("description", "Entry", "coordinate", "site_type"), sep = "@@")
  tt <- str_split(ss$orf, ";")
  ss0 <- list()
  for (i in seq_along(ss$inf)) {
    ss0[[i]] <- paste(ss$inf[i], tt[[i]], sep = "@@")
  }
  ss1 <- unlist(ss0)

  sce_site_refine0 <- data.frame(inf = ss1, stringsAsFactors = FALSE)
  sce_site_refine0 <- sce_site_refine0 %>% separate(., inf, into = c("description", "Entry", "coordinate", "site_type", "orf"), sep = "@@")
  sce_site_refine0$id_p <- paste(sce_site_refine0$orf, sce_site_refine0$coordinate, sep = "@")
  sce_site_refine0$description <- str_replace_all(sce_site_refine0$description, "\\.", "")
  sce_site_refine0$description <- str_trim(sce_site_refine0$description, side = "both")
  return(sce_site_refine0)
}



#' This a general function used to annotate the small SNP dataset
#'
#' @param snp_input a dataframe contains columns of "Chr","Pos","Gene","Ref","Alt" for each SNP
#' @param gene_feature a dataframe contains the detailed information of a gene
#'
#' @return
#' @export mutated_gene0 # a dataframe contains detailed annotation information of each SNP
#'
#' @examples
annotateSNP <- function(snp_input, gene_feature=gene_feature0) {
  mutated_test <- snp_input
  mutated_test$complement_sign <- getSingleMatchParameter(gene_feature$complement_sign, gene_feature$locus_tag, mutated_test$Gene2)
  mutated_test$Chr <- str_trim(mutated_test$Chr, side = "both")
  mutated_test$Pos <- as.numeric(mutated_test$Pos)
  mutated_test$Gene2 <- NA
  # get the gene name
  for (i in seq(length(mutated_test$Chr))){
    print(i)
    mutated_test$Gene2[i] <- getGeneName(mutated_test$Chr[i],mutated_test$Pos[i])
    print(getGeneName(mutated_test$Chr[i],mutated_test$Pos[i]))
  }

  mutated_test0 <- filter(mutated_test, Gene2 != "INTERGENIC") ##filter the mutated test
  #choose the metabolic gene
  #if the gene is type of "complement", then the complement_sign is "TRUE"
  #else the complement_sign is "FALSE"
  gene_feature0$complement_sign <- str_detect(gene_feature0$cds_location,"complement")
  index_m <- which(mutated_test0$Gene2 %in% gene_feature0$locus_tag ==TRUE)
  mutated_gene <- mutated_test0[index_m,]

  mutated_gene$Ref <- str_trim(mutated_gene$Ref, side = "both")
  mutated_gene$Alt <- str_trim(mutated_gene$Alt, side = "both")

  mutated_gene$complement_sign <- getSingleMatchParameter(gene_feature0$complement_sign,gene_feature0$locus_tag,mutated_gene$Gene2)
  mutated_gene0 <- mutated_gene

  for (i in seq(length(mutated_gene0$Chr))){
    if(mutated_gene0$complement_sign[i]){
      mutated_gene0$Ref[i] <- changeATCG(mutated_gene0$Ref[i])
      mutated_gene0$Alt[i] <- changeATCG(mutated_gene0$Alt[i])

    } else{
      mutated_gene0$Ref[i] <- mutated_gene0$Ref[i]
      mutated_gene0$Alt[i] <- mutated_gene0$Alt[i]
    }
  }

  return(mutated_gene0)
}



#' Parse the gene coordinate
#' Obtain a list which contains the coordinate information of each gene
#'
#' @param gene_name A string represent the gene name
#' @param genesum A dataframe contains the detailed gene annotation
#'
#' @return A list contains the coordinate of gene and its protein
#' @export
#'
#' @examples
#' data('gene_feature0')
#' getGeneCoordinate(gene_name ='YPR184W', genesum = gene_feature0)
getGeneCoordinate <- function(gene_name, genesum = gene_feature_GEM ){
  #genesum = gene_feature0
  #gene_name <- "YMR242C" # example
  ss <- filter(genesum, locus_tag==gene_name)
  gene_snp <- list()
  cds_seq <- vector()
  if(str_detect(ss$location[1], "complement")==FALSE){ #complement means the seq in the -1 strand
    ll <- ss$cds_location
    ll1 <- unlist(str_split(ll, ","))
    ll1 <- str_replace_all(ll1, "location=join\\(","") %>%
      str_replace_all(.,"\\)","") %>%
      str_replace_all(.,"location=", "")
    tt <- list()
    for (i in seq(length(ll1))){
      if(str_detect(ll1[i],"\\.\\." )){
        tt[[i]] <- unlist(str_split(ll1[i],"\\.\\."))
      } else{
        tt[[i]] <- ll1[i]
      }
    }

    tt0 <- list()
    for (i in seq(length(tt))){
      if(length(tt[[i]])==2) {
        tt0[[i]] <- seq(as.numeric(tt[[i]][1]), as.numeric(tt[[i]][2]),1)
      } else{
        tt0[[i]] <-tt[[i]][1]
      }
    }

    cds_seq <- unlist(tt0)

    gene_snp[['gene']] <- unlist(strsplit(ss$cds_seq[1], split = ""))
    # cds_location
    gene_snp[['gene_coordinate']] <- cds_seq
    gene_snp[['protein']] <- unlist(strsplit(ss$aa_seq[1], split = ""))
    gene_snp[['protein_coordinate']] <- seq(as.numeric(ss$aa_length[1]))
  } else{
    ss <- filter(genesum, locus_tag==gene_name)
    ll <- ss$cds_location
    ll1 <- unlist(str_split(ll, ",")) # coordinate of cds
    ll1 <- str_replace_all(ll1, "location=complement\\(join\\(","") %>%
      str_replace_all(.,"\\)","") %>%
      str_replace_all(.,"location=complement\\(", "")
    tt <- list()

    for (i in seq(length(ll1))){
      if(str_detect(ll1[i],"\\.\\." )){
        tt[[i]] <- unlist(str_split(ll1[i],"\\.\\."))
      } else{
        tt[[i]] <- ll1[i]
      }
    }

    tt0 <- list()

    for (i in seq(length(tt),1,-1)){
      #i=2
      if(length(tt[[i]])==2){ #code error here, original code:if(length(tt[[i]]==2)); new code:if(length(tt[[i]])==2)
        tt0[[i]] <- seq(as.numeric(tt[[i]][2]), as.numeric(tt[[i]][1]),-1)
        cds_seq <- c(cds_seq, unlist(tt0[[i]]))
      } else{
        tt0[[i]] <- as.numeric(tt[[i]][1])
        cds_seq <- c(cds_seq, unlist(tt0[[i]]))
      }

    }

    gene_snp[['gene']] <- unlist(strsplit(ss$cds_seq[1], split = ""))
    # cds_location
    gene_snp[['gene_coordinate']] <- cds_seq
    gene_snp[['protein']] <- unlist(strsplit(ss$aa_seq[1], split = ""))
    gene_snp[['protein_coordinate']] <- seq(as.numeric(ss$aa_length[1]))

  }

  return(gene_snp)

}











