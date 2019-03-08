#' fast estabolish the mapping relation between two element
#'
#' @param gene vetor
#' @param rxn vector
#' @param sep0 string
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' gene <- c('a&b','c')
#' rxn <- c('r1','r2')
#' splitAndCombine(gene, rxn, sep0="&")
#'
splitAndCombine <- function(gene, rxn,sep0) { ##one rxn has several genes, this function was used to splite the genes
  library(stringr)
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")

  }

  gene1 <- unique(unlist(gene0))
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }

  return(rxnGene)
}

