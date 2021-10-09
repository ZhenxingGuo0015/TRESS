### peak calling for real data when there are only one replicate
CallPeaks.oneRep <- function(Counts,
                             bins,
                             sf = NULL,
                             WhichThreshold = "fdr_lfc",
                             pval.cutoff = 1e-05,
                             fdr.cutoff = 0.05,
                             lfc.cutoff = 0.7,
                             windlen = 5,
                             lowCount = 10){
  ### step 1: grasp bumps based on lfc or binomial test
  if(length(sf) == 0){
    sf = colSums(Counts)/median(colSums(Counts))
    }
  Pvals = rep(NA, nrow(Counts))
  idx = rowSums(Counts) > 0
  Pvals[idx] = 1 - pbinom(Counts[idx, 2],
                          rowSums(Counts[idx, ]),
                          prob = sf[2]/sum(sf))
  Pvals[!idx] = 1
  ### lfc
  c0 = mean(as.matrix(Counts), na.rm = TRUE)  ### pseudocount
  lfc = log((Counts[, 2]/sf[2] + c0)/(Counts[, 1]/sf[1] + c0))
  smooth.lfc <- mySmooth(lfc, windlen = windlen)
  x.vals = data.frame(pvals = Pvals,
                      fdr = p.adjust(Pvals, method = "fdr"),
                      lfc = lfc)
  tmp = findBumps(chr = bins$chr,
                  pos = bins$start,
                  strand = bins$strand,
                  x = x.vals,
                  use = WhichThreshold,
                  pval.cutoff = pval.cutoff,
                  fdr.cutoff = fdr.cutoff,
                  lfc.cutoff = lfc.cutoff,
                  count = Counts)
  Bumps = tmp[tmp$counts > lowCount, ]
  if(nrow(Bumps) >= 2){
   ### step 2: binomial test based on the counts of each bumps
    peaks = Bumps[, c("chr", "start",  "end", "strand", "summit")]
    bins.GR = GRanges(Rle(bins$chr), IRanges(bins$start, bins$end))
    count = getPeakCounts(peaks = peaks,
                          allCounts = Counts,
                          allBins = bins.GR)
    pvals = 1 - pbinom(count[, 2],
                       rowSums(count),
                       prob = sf[2]/sum(sf))
    fdr = p.adjust(pvals, method = "fdr")
    peaks = cbind(peaks, count)
    peaks$pvals = pvals
    peaks$p.adj = fdr

    ### calculate log fold change for peak regions
    thiscount = peaks[, which(grepl("bam", colnames(peaks)) |
                                grepl("rep", colnames(peaks)) ) ]
    c0 = mean(as.matrix(thiscount), na.rm = TRUE)
    lfc = log((thiscount[, 2]/sf[2] + c0)/(thiscount[, 1]/sf[1] + c0))
    if(length(lfc) > windlen){
      smooth.lfc <- mySmooth(lfc, windlen = windlen)
      }else{
        smooth.lfc = lfc
        }
    peaks$lg.fc = smooth.lfc
    peaks = peaks[order(peaks$lg.fc, decreasing = TRUE), ]
    return(peaks)
    }else{
      cat("Less than 2 peaks!", sep = "\n")
      peaks = Bumps
      return(peaks)
    }
  }
