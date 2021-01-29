#' Telecharge les accessions list des chromosomes 1 et 2 de vibrio cholerae sur NCBI
#'
#' @param outDir output directory
#'
#' @return Ecrit dans le repertoire un fichier csv contenant les accessions listes
#' @export
#'
#' @examples
download.accession.chromosome <- function(outDir)
{
  library(Biostrings)
  library(reutils)
  library(seqinr)
  for (i in 1:2) {

    demo.search <- esearch(term = paste0("Vibrio cholerae[orgn] and chromosome ",i," [title] and complete [title]"), db = 'nuccore', usehistory = TRUE) #search
    accessions <- efetch(demo.search, rettype = "acc",retmode = "text",outfile= paste0(outDir,'vcholerae_c',i,".csv"))#fetch accessions
    accession  <- read.csv(paste0(outDir,"vcholerae_c",i,".csv"),header=F,stringsAsFactors = F)
    accession$V1 <- unlist(lapply(strsplit(accession$V1,split='\\.'),function(x) x[[1]]))
    write.csv(accession,paste0(outDir,'vcholerae_c',i,".csv"),row.names = F)
    print(i)

  }
}
