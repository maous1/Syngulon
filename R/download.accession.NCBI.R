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
download.accession.NCBI <- function(species,title,accessionDir,index)
{
  library(Biostrings)
  library(reutils)
  species2 = gsub("_"," ",species)
  for (i in index:length(species)) {

  demo.search <- esearch(term = paste0(species2[i],"[orgn] and ",title,"[title]"), db = 'nuccore', usehistory = TRUE) #search
  accessions <- efetch(demo.search, rettype = "acc",retmode = "text",outfile= paste0(accessionDir,species[i],".csv"))#fetch accessions
  accession  <- read.csv(paste0('01-accession-list/',species[i],".csv"),header=F,stringsAsFactors = F)
  accession$V1 <- unlist(lapply(strsplit(accession$V1,split='\\.'),function(x) x[[1]]))
  write.csv(accession,paste0(accessionDir,species[i],".csv"),row.names = F)
  print(i)

}
}
