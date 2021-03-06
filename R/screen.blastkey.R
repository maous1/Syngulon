#' Screen a sequence using Blast. utiliser pour les keys
#'
#' You provide a reference and a querry and the function compute the percentage of the reference which is covered by the querry
#' @param reference the reference sequence that you want to screen. Fasta file in one or severa sequences
#' @param querry the querry sequence. Fasta file in one or severa sequences
#' @param min.pc.ident
#' @param min.pc.length
#' @param geneDir
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry and write dans le geneDir
#' la sequences des gènes similaires
#' @export
#'
#' @examples
screenBlastkey <- function (reference, querry,min.pc.ident,min.pc.length,geneDir)
{
  library(Biostrings)
  library(GenomicRanges)
  try(unlink("temp", recursive = TRUE))
  dir.create("temp")
  dir.create("temp/dbblast")
  myarg <- paste0("-in ", reference, " -out temp/dbblast/db -dbtype nucl")
  system2(command = "makeblastdb", args = myarg, stdout = F)

  myarg <- paste0("-query ", querry, " -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 1000 -outfmt \"7 qframe sacc bitscore pident length slen qstart qend sstart send\"")
  system2(command = "blastn", args = myarg)
  blast <- try(read.table("temp/blast.txt", comment.char = "#"), silent = T)

  if (class(blast) == "data.frame")
  {
    sequences <- readDNAStringSet(reference)
    genes <-names(sequences)
    gene.levels <- levels(factor(genes))


    colnames(blast) <- c("qframe","subj.access", "bitscore", "pident", "align.length", "subj.len","qstart","qend","sstart","send")
    genes <- blast$subj.access
    blast <- data.frame(genes,blast)
    blast$pc.length <- round(100*(blast$align.length/blast$subj.len))
    blast$pc.length[blast$pc.length>100] <- 100

    ngenes <- length(gene.levels)
    toreturn <- character(length =ngenes )
    boolean=1
    for(i in 1:ngenes)
    {
      blast.selected <- blast[blast$gene==gene.levels[i],]
      blast.selected <- blast.selected[(blast.selected$pident>min.pc.ident)&(blast.selected$pc.length>min.pc.length),]
      if(dim(blast.selected)[1]>0 |!(isEmpty(blast.selected$bitscore)))
      {
        if (boolean==1) {
          species = gsub("/","",querry)
          species = gsub("03-genome","",species)
          species = gsub(paste0(basename(querry)),"",species)
          dir.create(paste0(geneDir,species))
          querry = readDNAStringSet(querry)
          dir.create(paste0(geneDir,species,"/",names(querry)))
          boolean=0
        }

        geneclone = DNAStringSet(querry[[1]][blast.selected$qstart:blast.selected$qend])
        if (blast.selected$sstart>blast.selected$send) {
          geneclone = reverseComplement(geneclone)
        }
        names(geneclone)= paste0(species,":",gene.levels[i],":",names(querry))
        writeXStringSet(geneclone,paste0(geneDir,species,"/",names(querry),"/",gene.levels[i],".fasta"))
        toreturn1 <- unique(unlist(lapply(strsplit(blast.selected$subj.access,split=':'),function(x) x[[1]])))
        toreturn1 <- paste(toreturn1[1:min(3,length(toreturn1))],collapse='|')
        blast.selected <- blast.selected[sort.list(blast.selected$bitscore,decreasing=T),]
        blast.selected <- blast.selected[which.max(blast.selected$pc.length),]
        toreturn2  <- paste0('pc.id = ',round(blast.selected$pident),' ; pc.cov = ',blast.selected$pc.length)
        toreturn[i] <- paste(toreturn2,toreturn1,sep=';')
      }
      else{toreturn[i] <- ''}
    }

  }
  else{toreturn <- character(length =ngenes )}
  unlink('temp',recursive = T)
  return(toreturn)
}
