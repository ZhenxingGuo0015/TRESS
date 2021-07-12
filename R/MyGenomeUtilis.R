getPeakCounts <- function(peaks, allCounts, allBins){
  ### grasp reads counts for a list of peaks
  ### allCounts: counts for all sites across the whole genome
  ### all bins cut from the whole genome
  #require(GenomicFeatures)
  #require(GenomicRanges)
  #require(matrixStats)
  peak.GR = GRanges(Rle(peaks$chr),
                    IRanges(peaks$start, peaks$end))
  iii = findOverlaps(peak.GR, allBins)
  query = queryHits(iii)
  subject = subjectHits(iii)
  peakCounts = matrix(0, nrow = length(peak.GR),
                      ncol = ncol(allCounts))
 # for(i in 1:length(peak.GR)) {
  for(i in seq_len(length(peak.GR))) {
    ix = query == i
    if(sum(ix) == 0) next
    peakCounts[i,] = colSums(allCounts[subject[ix],,drop = FALSE])
  }
  colnames(peakCounts) = colnames(allCounts)
  return(peakCounts)
}

lfc.peak <- function(peak, sf, binCounts, bins.GR){
  ### peak: genomic sites
  ### sf: size factor for samples that peaks belong to
  ### binCounts: reads counts in all bins across the whloe
  ###             genome in samples that peaks belong to
  ### bin.GR: Grange object to store the genomic location of each bin
 # require(base)
 # require(matrixStats)
  thiscount = getPeakCounts(peaks = peak,
                            allCounts = binCounts,
                            allBins = bins.GR)
  blocks = seq(1, ncol(binCounts), 2)
  c0 = rep(0, length(blocks))
 # for (j in 1:length(blocks)) {
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


mySmooth <- function(vec, windlen = 3){
  ### smooth the value of vec
 # require(base)
  cum = cumsum(vec)
  n = length(cum)
  a = c(cum[((windlen - 1)/2 + 1):n], rep(cum[n], (windlen - 1)/2))
 # b = c(rep(0, (windlen - 1)/2 + 1), cum[1:(n - (windlen - 1)/2 -1)])
  b = c(rep(0, (windlen - 1)/2 + 1), cum[seq_len((n - (windlen - 1)/2 -1))])
  res = (a-b)/windlen
  return(res)
}





ShowOnePeak <- function (onePeak, allBins, binCounts, ext = 500,
                         ylim = c(0, 1))
{

 # require(matrixStats)
  #require(base)
  #### get the methylation level of each bin
  sf = colSums(binCounts)/median(colSums(binCounts))
  allcounts.norm = sweep(binCounts, 2, sf, FUN = "/")
  input.norm = allcounts.norm[, seq(1, ncol(allcounts.norm), 2)]
  P =  allcounts.norm[, seq(2, ncol(allcounts.norm), 2)
                      ]/(allcounts.norm[, seq(1, ncol(allcounts.norm), 2)] +
                           allcounts.norm[, seq(2, ncol(allcounts.norm), 2)])

  allchr = as.character(allBins$chr)
  allpos = allBins$start
  chr = as.character(onePeak$chr)
  ix.chr = which(allchr == chr)
  thispos = allpos[ix.chr]

  #### added on July 5, 2021
  thisStrand = allBins$strand[ix.chr]
  ####

  thisInput = input.norm[ix.chr, ]
  thisP = P[ix.chr, ]



  xlim = c(onePeak$start - ext, onePeak$end + ext)
 # ix1 = which(thispos <= xlim[2] & thispos >= xlim[1] ) commentted on July 5, 2021
  ix1 = which(thispos <= xlim[2] & thispos >= xlim[1] & thisStrand == onePeak$strand)
  nSample = ncol(P)
  y.cex = 1
  #sNames = paste0("Replicate ", 1:ncol(P))#sampleNames(BSobj)
  sNames = paste0("Replicate ", seq_len(ncol(P)))
  par(mfrow = c(nSample, 1), mar = c(2.5, 2.5, 1.6, 2.5),
      mgp = c(1.5, 0.5, 0))
 # for (i in 1:ncol(P)) {
  for (i in seq_len(ncol(P))) {
    plot(thispos[ix1], thisP[ix1, i], type = "h", col = "blue",
         axes = FALSE, lwd = 1.5, xlab = "", ylab = "", ylim = ylim,
         xlim = xlim, main = sNames[i])
    box(col = "black")
    axis(1, )
    axis(2, col = "blue", col.axis = "blue")
    mtext(chr, side = 1, line = 1.33, cex = y.cex)
    mtext("methyl%", side = 2, line = 1.33, col = "blue",
          cex = y.cex)
    thisN.norm = thisInput[ix1, i]/max(thisInput[ix1, ]) * ylim[2]

    lines(thispos[ix1], thisN.norm, type = "l", col = "gray",
          lwd = 1.5)
    axis(side = 4, at = seq(0, ylim[2], length.out = 5),
         labels = round(seq(0, max(thisInput[ix1, ]), length.out = 5)))
    mtext("Input read depth", side = 4, line = 1.33, cex = y.cex)
    rect(onePeak$start, ylim[1], onePeak$end, ylim[2], col = "#FF00001A",
         border = NA)
  }
}








