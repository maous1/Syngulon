#' Telecharge les accessions list sur NCBI
#'
#' @param species
#' @param title
#' @param accessionDir
#'
#' @return
#' @export
#'
#' @examples
download.accession.bacteriocin <- function(accessionDir)
{
  library(Biostrings)
  library(reutils)
  library(seqinr)

  demo.search <- esearch(term = paste0("bacteriocin [title]"), db = 'nuccore', usehistory = TRUE) #search
  accessions <- efetch(demo.search, rettype = "acc",retmode = "text",outfile= paste0(accessionDir,'bacteriocin.csv'))#fetch accessions
  accession  <- read.csv(paste0(accessionDir,'bacteriocin.csv'),header=F,stringsAsFactors = F)
  accession$V1 <- unlist(lapply(strsplit(accession$V1,split='\\.'),function(x) x[[1]]))
  write.csv(accession,paste0(accessionDir,'bacteriocin.csv'),row.names = F)
  print(i)

}
