# obtain peak read count
getPeakCounts <- function(peaks, allCounts, allBins){
  ### grasp reads counts for a list of peaks
  ### allCounts: counts for all sites across the whole genome
  ### all bins cut from the whole genome
  peak.GR = GRanges(Rle(peaks$chr),
                    IRanges(peaks$start, peaks$end))
  iii = findOverlaps(peak.GR, allBins)
  query = queryHits(iii)
  subject = subjectHits(iii)
  peakCounts = matrix(0, nrow = length(peak.GR),
                      ncol = ncol(allCounts))
  for(i in seq_len(length(peak.GR))) {
    ix = query == i
    if(sum(ix) == 0) next
    peakCounts[i,] = colSums(allCounts[subject[ix],,drop = FALSE])
  }
  colnames(peakCounts) = colnames(allCounts)
  return(peakCounts)
}
