library(readxl)
library(tidyverse)
library(hongR)
# main function

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



