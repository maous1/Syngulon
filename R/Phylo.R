#' fait un arbre phylogenetic. Utiliser avec de grand nombre de sequences
#'
#' @param bacteria.table
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
phylo <- function(sequences)
{
  library(muscle)
  library(phangorn)
  library(ggtree)
  maxlength <- max(width(sequences))

  align.muscle <- muscle::muscle(sequences)
  dist1 <- stringDist(as(align.muscle,"DNAStringSet"), method="hamming")
  dist1 <- 100*dist1/maxlength
  mytree1 <- upgma(dist1)

  nseqeuences <- length(sequences)

  pdf(paste0('phylogenetic.otus.pdf'),width=10,height = 10*log10(nseqeuences))
  p <- ggtree(mytree1)
  p <- p  + geom_tiplab(offset=0) + xlim(NA, 150) +geom_treescale(0.05,-2.5,width=10,fontsize = 2,linesize = 0.5)
  plot(p)
  dev.off()
}
