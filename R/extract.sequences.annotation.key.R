#' extrait les sequences comportant colicin et ne comportant pas resistance et immunity dans le produit des annotations utilisé pour les key
#'
#' @param species
#' @param annotationDir
#' @param outDir
#' @param genomeDir
#'
#' @return
#' @export
#'
#' @examples
extract.sequences.annotation <- function(species,annotationDir,outDir,genomeDir)
{
  library(dplyr)
  library(Biostrings)

  boolean = 1
  for (i in 1:length(species)) {
    listannotation = list.files(paste0(annotationDir,species[i]),full.names =T)
    print(i)
    if (!isEmpty(listannotation)) {
      for (j in 1:length(listannotation)) {
        currentannotation = data.frame(read.csv(listannotation[j]))
        currentannotation$name = unlist(rep(paste0(species[i],"-",basename(listannotation[j])),dim(currentannotation)[1]))
        currentannotation = currentannotation%>% filter(.,grepl("colicin",product))%>%filter(.,!grepl("immunity",product))%>% filter(.,!grepl("resistance",product))
        ecrire = DNAStringSet()
        if (dim(currentannotation)[1]>1) {
          dir.create(paste0(outDir,species[i]))
          for (k in 1:dim(currentannotation)[1]) {
            sequences = readDNAStringSet(paste0(genomeDir,species[i],"/",gsub(".csv",".fasta",basename(listannotation[j]))))
            gene = DNAStringSet(sequences[[1]][currentannotation$start[k]:currentannotation$end[k]])
            if (currentannotation$start[k]>currentannotation$end[k]) {
              gene = reverseComplement(gene)
              gene = reverse(gene)
            }
            names(gene)= names(sequences)
            ecrire = c(ecrire,gene)
          }
          writeXStringSet(ecrire,paste0(outDir,species[i],"/",names(sequences),".fasta"))
        }

        if (boolean==0 & dim(currentannotation)[1]>1) {
          currentannotationfull = add_row(currentannotationfull,currentannotation)
        }
        if (boolean==1 & dim(currentannotation)[1]>1) {
          currentannotationfull = currentannotation
          boolean =0
        }

      }
    }
  }
  return(currentannotationfull)
}
