### Visualize single peak
ShowOnePeak <- function(onePeak, allBins, binCounts,
                        isDMR = FALSE,
                        Sname = NULL,
                        ext = 500,
                        ylim = c(0, 1)
                        ){
  #### get the methylation level of each bin
  sf = colSums(binCounts)/median(colSums(binCounts))
  allcounts.norm = sweep(binCounts, 2, sf, FUN = "/")
  id.input = seq(1, ncol(allcounts.norm), 2)
  id.ip = seq(2, ncol(allcounts.norm), 2)

  input.norm = allcounts.norm[, id.input]
  P = allcounts.norm[, id.ip]/(allcounts.norm[, id.input] +
                                 allcounts.norm[, id.ip])

  allchr = as.character(allBins$chr)
  allpos = allBins$start
  chr = as.character(onePeak$chr)
  ix.chr = which(allchr == chr)
  thispos = allpos[ix.chr]
  thisStrand = allBins$strand[ix.chr]
  thisInput = input.norm[ix.chr, ]
  thisP = P[ix.chr, ]

  xlim = c(onePeak$start - ext, onePeak$end + ext)
  # ix1 = which(thispos <= xlim[2] & thispos >= xlim[1] )
  # commentted on July 5, 2021
  ix1 = which(thispos <= xlim[2] & thispos >= xlim[1] &
                thisStrand == onePeak$strand)
  nSample = ncol(P)
  y.cex = 1
  if(!isDMR){
    sNames = paste0("Replicate ", seq_len(ncol(P)))
  }else{
    sNames = Sname
  }
  par(mfrow = c(nSample, 1), mar = c(2.5, 2.5, 1.6, 2.5),
      mgp = c(1.5, 0.5, 0))
  for (i in seq_len(ncol(P))) {
    plot(thispos[ix1], thisP[ix1, i],
         type = "h", col = "blue",
         axes = FALSE, lwd = 1.5,
         xlab = "", ylab = "", ylim = ylim,
         xlim = xlim, main = sNames[i])
    box(col = "black")
    axis(1, )
    axis(2, col = "blue", col.axis = "blue")
    mtext(chr, side = 1, line = 1.33, cex = y.cex)
    mtext("methyl%", side = 2, line = 1.33, col = "blue",
          cex = y.cex)
    thisN.norm = thisInput[ix1, i]/max(thisInput[ix1, ])*ylim[2]
    lines(thispos[ix1], thisN.norm,
          type = "l", col = "gray",lwd = 1.5)
    axis(side = 4, at = seq(0, ylim[2], length.out = 5),
         labels = round(seq(0, max(thisInput[ix1, ]),
                            length.out = 5)))
    mtext("Input read depth",
          side = 4, line = 1.33, cex = y.cex)
    rect(onePeak$start, ylim[1], onePeak$end, ylim[2],
         col = "#FF00001A", border = NA)
  }
}
