M6Apeak.MultiRep.step2 <- function(Candidates, sf, mu.cutoff,
                                   WhichThreshold = "fdr_lfc",
                                   pval.cutoff = 1e-5,
                                   fdr.cutoff = 0.05,
                                   lfc.cutoff = 0.7
                                   ){

  ### Order candidate peaks from step 1
  ### with more sophisticated statistical model
  Candidates$score = NULL   ## remove score from step 1
  thispeak = Candidates
  idx = which(grepl("rep", colnames(thispeak)) |
                grepl("bam", colnames(thispeak)))
  thiscount = thispeak[, idx]
  res = M6Apeak.paramEsti(mat = as.matrix(thiscount),
                          sf = sf,
                          cutoff = mu.cutoff)
  thispeak = cbind(Candidates, res)

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
