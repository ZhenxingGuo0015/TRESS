CallPeaks.multiRep <- function(Candidates, mu.cutoff,
                               WhichThreshold = "fdr_lfc",
                               pval.cutoff = 1e-5,
                               fdr.cutoff = 0.05,
                               lfc.cutoff = 0.7
                               ){

  ### Detect and rank peaks from candidates
  ### with more sophisticated statistical model
  thispeak = Candidates$Regions
  thiscount = Candidates$Counts
  if(is.null(Candidates$sf)){
    Candidates$sf = colSums(thiscount)/median(colSums(thiscount))
  }
  if(!is.null(Candidates$lg.fc)){
    lg.fc = Candidates$lg.fc
  }else{
    lg.fc = getLogFC(Counts = thiscount, sf = Candidates$sf)
  }

  res = CallPeaks.paramEsti(mat = as.matrix(thiscount),
                            sf = Candidates$sf,
                            cutoff = mu.cutoff)
  thispeak = cbind(thispeak, lg.fc, thiscount, res)

  #### filter out insignificant peaks, added on Jan 28, 2021
  if(WhichThreshold == "pval"){
    idx <- which(thispeak$pvals < pval.cutoff)
  }else if(WhichThreshold == "fdr"){
    idx <- which(thispeak$p.adj < fdr.cutoff)
  }else if(WhichThreshold == "lfc"){
    idx <- which(thispeak$lg.fc >= lfc.cutoff)
  } else if(WhichThreshold == "pval_lfc"){
    idx <- which(thispeak$pvals < pval.cutoff &
                   thispeak$lg.fc >= lfc.cutoff)
  }else if(WhichThreshold == "fdr_lfc"){
    idx <- which(thispeak$p.adj < fdr.cutoff &
                   thispeak$lg.fc >= lfc.cutoff)
  }
  peaks = thispeak[idx, ]
  ####
  #### rank the peak list based on their statistics
  peaks = peaks[order(peaks$rScore, decreasing = TRUE), ]
  #### ended editing
  cat("The number of Peaks in this sample is: ",
      nrow(peaks), sep = "\n")
  return(peaks)
}
