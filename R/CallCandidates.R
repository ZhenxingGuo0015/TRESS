CallCandidates <- function(Counts, bins,
                           WhichThreshold ="fdr_lfc",
                           pval.cutoff = 1e-5,
                           fdr.cutoff = 0.05,
                           lfc.cutoff = 0.7,
                           windlen = 5,
                           lowcount = 30){
  ### 1. find bumps for each replicate based on binomial test
  sf = colSums(Counts)/median(colSums(Counts))
  sx = sf[seq(1, length(sf), 2)]; sy = sf[seq(2, length(sf), 2)]
  bn.prob = sy/(sx + sy)

  blocks = seq(1, ncol(Counts), 2)
  c0 = rep(0, length(blocks))
  Pvals = matrix(0, nrow = nrow(Counts), ncol = length(blocks))
  Bumps = vector("list", length = length(blocks))
  for (j in seq_along(blocks)) {
    #cat(j, sep = "\n")
    id = blocks[j]
    dat = Counts[, id:(id+1)]
    thissf = sf[id:(id+1)]
    ### pvals based on binomial test
    idx = rowSums(dat) > 0
    Pvals[idx, j] = 1 - pbinom(dat[idx, 2],
                               rowSums(dat[idx, ]),
                               prob = bn.prob[j])
    Pvals[!idx, j] = 1

    ### lfc
    c0[j] = mean(as.matrix(dat), na.rm = TRUE)  ### pseudocount
    lfc = log((dat[, 2]/thissf[2] + c0[j]
               )/(dat[, 1]/thissf[1] + c0[j]))
    smooth.lfc <- mySmooth(lfc, windlen = windlen)
    x.vals = data.frame(pvals = Pvals[, j],
                        fdr = p.adjust(Pvals[, j], method = "fdr"),
                        lfc = lfc)
    ### find bumps based on pvals, fdr or lfc
    tmp = findBumps(chr = bins$chr,
                    pos = bins$start,
                    strand = bins$strand,
                    x = x.vals,
                    use = WhichThreshold,
                    pval.cutoff = pval.cutoff,
                    fdr.cutoff = fdr.cutoff,
                    lfc.cutoff = lfc.cutoff,
                    count = dat)
    Bumps[[j]] = tmp[tmp$counts > 10, ]
    }

  ### 2. merge bumps from different replicates
  cat("Merge bumps from different replicates...", sep = "\n")
  thispeak = NULL
  #  for (i in 1:(length(Bumps))) {
  for (i in seq_len(length(Bumps))) {
   # cat(paste0("Bumps ", i), sep = "\n")
    thisBump = Bumps[[i]]
    thisBump$strand[thisBump$strand=="."] ="*"
    thisBump.GR = GRanges(Rle(thisBump$chr),
                          IRanges(thisBump$start, thisBump$end),
                          Rle(thisBump$strand))
    if(is.null(thispeak)){
      thispeak = thisBump.GR
      }else{
        thispeak = union(thispeak, thisBump.GR)
      }
    }
  thispeak = as.data.frame(thispeak)
  colnames(thispeak) = c("chr", "start", "end", "width", "strand")

  if(nrow(thispeak) >= 1){
    ### add summit in bumps to peak
    peak.summit = rep(NA, nrow(thispeak))
    for (i in seq_len(length(Bumps))) {
      thisBump = Bumps[[i]]
      thisBump$strand[thisBump$strand=="."] = "*"
      bump.GR = GRanges(Rle(thisBump$chr),
                        IRanges(thisBump$start, thisBump$end),
                        Rle(thisBump$strand))
      peak.GR = GRanges(Rle(thispeak$chr),
                        IRanges(thispeak$start, thispeak$end),
                        Rle(thispeak$strand))
      hits = findOverlaps(peak.GR, bump.GR)
      if(length(hits) > 0){
        idx.query = unique(queryHits(hits))
        idx.subject = subjectHits(hits)
        for (j in seq_along(idx.query)) {
          thisquery = idx.query[j]
          thissubject = idx.subject[which(queryHits(hits) == thisquery)]
          if(is.na(peak.summit[idx.query[j]]) ){
            peak.summit[idx.query[j]] =
              paste0(Bumps[[i]]$summit[thissubject],
                     collapse = "_")
          }
        }
      }
    }
    thispeak$summit = peak.summit
    thispeak$width = NULL

    #### 3. get reads count in candidate regions
    bins.GR = GRanges(Rle(bins$chr), IRanges(bins$start,
                                             bins$end), Rle(bins$strand))
    thiscount = getPeakCounts(peaks = thispeak,
                              allCounts = Counts,
                              allBins = bins.GR)
    ### add lfc, only used for peak calling
    if(nrow(thispeak) > 1){
      myMeans = rowMeans
    }else{
      myMeans = mean
    }
    idx.x = seq(1, ncol(thiscount), 2)
    idx.y = seq(2, ncol(thiscount), 2)
    lfc = log(myMeans((thiscount[, idx.y]/sf[idx.y] + c0
    )/(thiscount[, idx.x]/sf[idx.x] + c0),
    na.rm = TRUE))
    if(length(lfc) > windlen){
      smooth.lfc <- mySmooth(lfc, windlen = windlen)
    }else{
      smooth.lfc = lfc
    }
    #####

    ##########
    Candidates = list(Regions = thispeak, Counts = thiscount,
                      lg.fc = smooth.lfc,
                      sf = sf)
    if(nrow(thiscount) >= 2){
      idx = which(rowSums(thiscount[, seq(1, ncol(thiscount), 2)])
                  > lowcount &
                    rowSums(thiscount[, seq(2, ncol(thiscount), 2)])
                  > lowcount)
      tmp1 = thispeak[idx, ];tmp2 = thiscount[idx, ]; tmp3 = smooth.lfc[idx]
      rownames(tmp1) = rownames(tmp2) = as.factor(seq_len(nrow(tmp1)))
      names(tmp3) = as.factor(seq_len(nrow(tmp1)))
      Candidates = list(Regions = tmp1, Counts = tmp2,
                        lg.fc = tmp3, sf = sf)
    }
   
  }else{
    Candidates = list()
    cat("No candidates from any samples!", sep = "\n")
  }
  return(Candidates)
}
