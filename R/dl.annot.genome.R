#' telecharge les fichiers d'annotations et des genomes
#'
#' @param species
#' @param NmaxPlasmid
#' @param annotationDir
#' @param genomeDir
#' @param accessionDir
#'
#' @return
#' @export
#'
#' @examples
dl.annot.genome <- function(species,NmaxPlasmid=1000, annotationDir,genomeDir,accessionDir)
{
  for (j in 1:length(species)) {
    dir.create(paste0(annotationDir,species[j]))
    dir.create(paste0(genomeDir,species[j]))
    accession <- read.csv(paste0(accessionDir,species[j],".csv"),stringsAsFactors = F)
    accession  <- accession$V1
    if (length(accession)>NmaxPlasmid) {
      accession <- accession[sample(1:length(accession), NmaxPlasmid, replace=F)]
    }
    N.accession <- length(accession)
    for(i in 1:N.accession)
    {
      seq <- try(read.GenBank(accession[i]),silent = T)
      current.annotation <- try(getAnnotationsGenBank(access.nb=accession[i], quiet = TRUE))
      if (class(seq)=="DNAbin" & class(current.annotation)=="data.frame"& (all(is.element(c('start','end','type','product'),colnames(current.annotation))))) {
        write.FASTA(seq,paste0(genomeDir,species[j],"/",accession[i],'.fasta'))
        seq <- readDNAStringSet(paste0(genomeDir,species[j],"/",accession[i],'.fasta'))
        writeXStringSet(seq,paste0(genomeDir,species[j],"/",accession[i],'.fasta'))
        current.annotation <- tibble(current.annotation)
        current.annotation <- current.annotation %>% select(start,end,type,gene,product)
        write.csv(current.annotation,paste0(annotationDir,species[j],"/",accession[i],'.csv'),row.names = F)
      }

      print(i)
    }
  }
}
