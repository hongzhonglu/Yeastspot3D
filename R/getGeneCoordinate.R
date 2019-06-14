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
