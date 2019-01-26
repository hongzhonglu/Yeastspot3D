#-----------------------------------------------------------------
# hotspot: merege multiple file of hotspot analysis result into one
#-----------------------------------------------------------------

#' Title
#'
#' @param infile the directory contains the hotspot analysis result
#' @param outfile the directory contains the hotspot analysis result
#'
#' @return
#' @export
#'
#' @examples
mergeHotspot <- function(infile, outfile) {

  file_list <- list.files(infile)

  for (file in file_list) {
    file <- paste(infile, "/", file, sep = "")
    # if the merged dataset doesn't exist, create it
    if (!exists("hotspot_merge")) {
      hotspot_merge <- read.table(file, header = TRUE, sep = "\t")
    } else {
      temp_dataset <- read.table(file, header = TRUE, sep = "\t")
      hotspot_merge <- rbind(hotspot_merge, temp_dataset)
      rm(temp_dataset)
    }
  }
  # output a new file named hotspot_merge.txt which contains all the results from hotspots

  write.table(hotspot_merge, paste(outfile, "/", "hotspot_merge.txt", sep = ""), row.names = FALSE, sep = "\t")
}
