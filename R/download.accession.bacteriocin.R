#' Telecharge les accessions de bacteriocin sur NCBI
#'
#' @param outDir The output directory
#'
#' @return Ecrit dans le repertoire un fichier csv contenant les accessions listes
#' @export
#'
#' @examples
download.accession.bacteriocin <- function(outDir)
{
  library(Biostrings)
  library(reutils)
  library(seqinr)

  demo.search <- esearch(term = paste0("bacteriocin [title]"), db = 'nuccore', usehistory = TRUE) #search
  accessions <- efetch(demo.search, rettype = "acc",retmode = "text",outfile= paste0(outDir,'bacteriocin.csv'))#fetch accessions
  accession  <- read.csv(paste0(outDir,'bacteriocin.csv'),header=F,stringsAsFactors = F)
  accession$V1 <- unlist(lapply(strsplit(accession$V1,split='\\.'),function(x) x[[1]]))
  write.csv(accession,paste0(outDir,'bacteriocin.csv'),row.names = F)
  print(i)

}
