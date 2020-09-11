#' fait un arbre phylogenetic. Utiliser avec de grand nombre de sequences
#'
#' @param bacteria.table
#' @param Dir
#'
#' @return
#' @export
#'
#' @examples
phylo.all <- function(bacteria.table,Dir)
  {
  correspondance.organism.subgroup <- bacteria.table%>%group_by(Organism) %>% dplyr::count(Organism, SubGroup) %>% dplyr::slice(which.max(n)) %>% dplyr::rename(species=Organism)
  filelist <- list.files(Dir,full.names = T,recursive = T)
  filelist <- filelist[grep(pattern = '.fasta',x =filelist )]
  genename <- gsub(basename(filelist),pattern = '.fasta',replacement = '')
  ngenes <- length(genename)

  sequences <- readDNAStringSet(filelist[i])
  maxlength <- max(width(sequences))

  align.muscle <- muscle::muscle(sequences)
  dist1 <- stringDist(as(align.muscle,"DNAStringSet"), method="hamming")
  dist1 <- 100*dist1/maxlength
  mytree1 <- upgma(dist1)

  nseqeuences <- length(sequences)

  pdf(paste0('99-results/phylogenetic.otus.',genename[i],'.pdf'),width=10,height = 10*log10(nseqeuences))
  p <- ggtree(mytree1)
  p <- p  + geom_tiplab(offset=0) + xlim(NA, 150) +geom_treescale(0.05,-2.5,width=10,fontsize = 2,linesize = 0.5)
  plot(p)
  dev.off()

}
