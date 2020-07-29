#' extract a table of the phylum of bacteria based on the taxize package
#'
#' to compute later
#'
#' @return a table with the one line for each phylum
#' @export
#'
extract.phylum <- function()
{
  library(taxize)
  library(dplyr)
  phylum <- itis_downstream(id = 50, downto="phylum")
  phylum <- tibble(phylum)
  phylum <- phylum %>% select(tsn,taxonname)
  return(phylum)
}


