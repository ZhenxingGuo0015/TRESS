### get bin strand given bin sites and genome
getStrand <- function(bins, anno_TXDB){
  # anno_TXDB: an annotation file in form of txdb
  txdb = anno_TXDB
  allGenes = genes(txdb)
  allGenes.strand = as.character(strand(allGenes))

  tmp = findOverlaps(bins, allGenes)
  length(tmp) ## there are bins overlapping more than one gene.
  binIdx = queryHits(tmp)
  geneIdx = subjectHits(tmp)

  tt = table(binIdx)
  binNames = as.integer(names(tt))

  ### start to get strands
  binStrand = rep(".", length(bins))

  ## unique ones
  ix.unique = binNames[which(tt==1)]
  tmpidx = which(binIdx %in% ix.unique)
  geneidx.unique = geneIdx[tmpidx]
  binStrand[ix.unique] = allGenes.strand[geneidx.unique]

  ## deal with bins overlapping multiple genes
  ix.dup = binNames[which(tt>1)]
  tmpidx = which(binIdx %in% ix.dup)
  geneidx.dup = geneIdx[tmpidx]
  strands.dup = allGenes.strand[geneidx.dup]
  binIdx.dup = binIdx[tmpidx]
  res <- tapply(strands.dup, binIdx.dup,
                function(x) {
                  if(all(x == "+")) return("+")
                  else if(all(x == "-")) return("-")
                  else return("*")
                })
  binStrand[ix.dup] = res
  return(binStrand)
}
