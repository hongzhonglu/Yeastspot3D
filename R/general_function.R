#' Get multiple mapping results for one element
#'
#' @param description A vector contains the description for each element from 'reaction'
#' @param reaction A vector contrains the element been mapped
#' @param ko A vector contains the element need mapping
#'
#' @return A vector, which contains several mapping results seperated by ";"
#' @export
#' @examples
#' w <- 1:6
#' v <- c('z','a','b','a','b', 'e')
#' testData <- c('a','b','g')
#' getMultipleMatchParameter(w,v,testData)
#'
getMultipleMatchParameter <- function(description, reaction, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length( which (reaction %in%  ko[i]))){
      index <- which (reaction %in%  ko[i])
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{

      result[i] <- NA
    }
  }
  return(result)
}



#' Get single mapping result for one element
#'
#' @param description A vector contains the description for each element from 'reaction'
#' @param reaction A vector contrains the element been mapped
#' @param ko A vector contains the element need mapping
#'
#' @return A vector
#' @export
#' @examples
#' w <- 1:6
#' v <- c('z','a','b','a','b', 'e')
#' testData <- c('a','b','g')
#' getSingleMatchParameter(w,v,testData)
#'
getSingleMatchParameter <- function(description, reaction, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction))){
      index <- match(ko[i],reaction)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{

      result[i] <- NA
    }
  }
  return(result)
}



#' Fast estabolish the mapping relation between two element
#'
#' @param gene A vetor contains elements with "&" or ","
#' @param rxn A vector contains unique identifier
#' @param sep0 A string, i.e, "&". ","
#'
#' @return A dataframe contains the single mapping relation for the element from gene and rxn
#' @export
#' @examples
#' gene <- c('a&b','c')
#' rxn <- c('r1','r2')
#' splitAndCombine(gene, rxn, sep0="&")
#'
splitAndCombine <- function(gene, rxn, sep0) {
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
