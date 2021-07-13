## log fold change of each peak
lfc.peak <- function(peak, sf, binCounts, bins.GR){
  ### peak: genomic sites
  ### sf: size factor for samples that peaks belong to
  ### binCounts: reads counts in all bins across the whloe
  ###             genome in samples that peaks belong to
  ### bin.GR: Grange object to store the genomic location of each bin
  thiscount = getPeakCounts(peaks = peak,
                            allCounts = binCounts,
                            allBins = bins.GR)
  blocks = seq(1, ncol(binCounts), 2)
  c0 = rep(0, length(blocks))
  for (j in seq_len(length(blocks))) {
    id = blocks[j]
    dat = binCounts[, id:(id+1)]
    c0[j] = mean(as.matrix(dat), na.rm = TRUE)
  }
  idx.x = seq(1, ncol(binCounts), 2)
  idx.y = seq(2, ncol(binCounts), 2)
  lfc = log(rowMeans((thiscount[, idx.y]/sf[idx.y] + c0
  )/(thiscount[, idx.x]/sf[idx.x] + c0),
  na.rm = TRUE))
  smooth.lfc <- mySmooth(lfc, windlen = 5)
  ###
  smooth.lfc
}
