### divide genes into bins
exonBins.byTXDB <- function(txdb,
                            binsize = 50,
                            IncludeIntron = FALSE){
  #### each gene will be divided into bins of length "binsize"
  ### when the annotation file are provided by
  ###    user instead of downloading from internet
  tmp = unique(genes(txdb))
  idx = (grepl("random", tmp@seqnames)|
           grepl("M", tmp@seqnames)|
           grepl("hap", tmp@seqnames)|
           grepl("_", tmp@seqnames))
  allgenes = tmp[!idx]
  unique.gene = union(allgenes, allgenes) # combine overlapped gene
  wins = slidingWindows(unique.gene, width = binsize, step = binsize)
  allbins = unlist(wins)
  allbins = allbins[width(allbins) == binsize]

  allExons = unlist(exonsBy(txdb))
  if(!IncludeIntron){
    ## get bins on the exons
    ix = (countOverlaps(allbins, allExons) > 0 )
    bins.exons = allbins[ix]
    invisible(return(list(bins=bins.exons,
                          keep.bin = ix,
                          allExons=allExons)))
  }else{
    ## get bins on the exons and introns
    allIntron = unlist(intronsByTranscript(txdb))
    ix = ((countOverlaps(allbins, allExons)+
             countOverlaps(allbins, allIntron)) > 0)

    bins.exIn = allbins[ix]
    invisible(return(list(bins=bins.exIn,
                          keep.bin = ix,
                          allExons=allExons,
                          allIntrons = allIntron)))
  }

}
