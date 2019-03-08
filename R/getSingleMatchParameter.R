#' get single mapping result for one element
#'
#' @param description vector contains the description for each element from 'reaction'
#' @param reaction vector contrains the element been mapped
#' @param ko vector contains the element need mapping
#'
#' @return one element
#' @export
#'
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
#devtools::document()
