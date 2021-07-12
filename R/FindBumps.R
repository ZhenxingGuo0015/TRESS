## should consider chr, but this should be okay most of the time
#require(base)
#require(stats)
#require(matrixStats)
findRegion <- function(chr, pos, sep=1000) {
  pos.diff <- abs(c(as.integer(0), diff(pos)))
  idx.jump <- which(pos.diff>sep)
  regions <- rbind(c(1, idx.jump),
                   c(idx.jump-1, length(pos)))
  regions
}

findBumps <- function(chr, pos, strand, x, count,
                      use = "pval",
                      pval.cutoff,
                      fdr.cutoff,
                      lfc.cutoff,
                      sep = 2000, minlen=100,
                      minCount=3, dis.merge=100,
                      scorefun = mean, sort=TRUE) {
  ### This function is used to combine significant bins to form bumps
  ### given the location and p-value of each bin

  if(sep < dis.merge)
    sep = dis.merge + 1

  if(use == "pval"){
    flag <- as.numeric(x$pvals < pval.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "fdr"){
    flag <- as.numeric(x$fdr < fdr.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "lfc"){
    flag <- as.numeric(x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  } else if(use == "pval_lfc"){
    flag <- as.numeric(x$pvals < pval.cutoff &
                         x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "fdr_lfc"){
    flag <- as.numeric(x$fdr < fdr.cutoff &
                         x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  }

  ## divide the whole genome into consecutive regions, which were sequenced
  regions <- findRegion(chr, pos, sep)
  ## loop on regions: within each region,
  ## combine significant bins to form a bump.
  ## there may be multiple bumps (candidate peaks)
  ## within one region
  initn <- 100000
  result <- data.frame(chr=rep("chr1",initn), start=rep(0, initn),
                       end=rep(0, initn), length=rep(0, initn),
                       strand = rep(NA, initn),
                       summit = rep(0, initn), counts = rep(0, initn),
                       score=rep(0, initn))
  levels(result[,1]) <- unique(chr)
  result.idx <- 0
  Peakstrand <- NULL
 # for(i in 1:ncol(regions)) {
  for(i in seq_len(ncol(regions))) {
    idx <- regions[1,i]:regions[2, i]
    pos.region <- pos[idx]
    strand.region <- strand[idx]
    count.region <- count[idx, ]
    if(length(idx)<minCount) next
    nn <- length(idx)
    flag.region <- flag[idx]
    ## get start/end position
    # original code
    # startidx <- which(flag.region[-nn]==0 & flag.region[-1]==1)+1
    ### modified on Jan 29, to consider strand
    startidx <- which((flag.region[-nn]==0 & flag.region[-1]==1)|
                        (strand.region[-1]!=strand.region[-nn] &
                           flag.region[-1]==1)
                      )
    startidx <- startidx + 1
    ####
    if(flag.region[1]==1)
      startidx <- c(1, startidx)
    if(length(startidx)==0)
      next

    ### original end idx
    #endidx <- which(flag.region[-nn]==1 & flag.region[-1]==0)
    #### modified on Jan 29, 2021 to consider strand
    endidx <- which((flag.region[-nn]==1 & flag.region[-1]==0) |
                      (strand.region[-nn]!= strand.region[-1] &
                         flag.region[-nn]==1))
    ####

    if(flag.region[nn]==1)
      endidx <- c(endidx, nn)




    ## remove if there are less than minCount probes
    idx.keep <- (endidx-startidx+1)>=minCount
    startidx <- startidx[idx.keep]
    endidx <- endidx[idx.keep]
    if(length(endidx)==0) next

    # ## merge if they are really close
    nbump <- length(startidx)
    if(nbump>1) {
      bumppos <- cbind(pos[idx][startidx], pos[idx][endidx])
      dis <- bumppos[-1,1]>(bumppos[-nbump,2]+dis.merge)
      idx.start <- which(c(1,dis)==1)
      idx.end <- which(c(dis,1)==1)
      ## merged
      startidx <- startidx[idx.start]
      endidx <- endidx[idx.end]
    }
    nbump <- length(startidx)
    ll <- pos.region[endidx] - pos.region[startidx] + 1 ### length of each bump
    tmpn <- length(ll)
    # ## make bump scores
    x.thisregion <- x[idx, ]
    scores.thisregion <- rep(0, nbump)
    counts.thisregion <- rep(0, nbump)
    summit.thisregion <- rep(0, nbump)
    strand.thisregion <- rep(0, nbump)
   # for(ibump in 1:nbump){
    for(ibump in seq_len(nbump)){
      thisrange <- startidx[ibump]:endidx[ibump]
      scores.thisregion[ibump] <- scorefun(x.thisregion$pval[thisrange])
      counts.thisregion[ibump] <- sum(count.region[thisrange,])
      thispos <- pos.region[thisrange]
      thistrand <- strand.region[thisrange]
      strand.thisregion[ibump] <- defineStrand(thistrand)
      summit.thisregion[ibump] <- thispos[
        which.max(count.region[thisrange, 2])
        ]+24
    }

    #result[result.idx+(1:tmpn),] <-
    result[result.idx+(seq_len(tmpn)),] <-
      data.frame(chr = as.character(chr[idx][startidx]),
                 start = pos[idx][startidx],
                 end = pos[idx][endidx]+ 49,
                 length=ll+49,
                 strand = as.character(strand.thisregion),
                 summit = summit.thisregion,
                 count = counts.thisregion,
                 score = scores.thisregion)
    Peakstrand = c(Peakstrand, strand.thisregion)

    result.idx <- result.idx + tmpn
  }

  #result <- result[1:result.idx,]
  result <- result[seq_len(result.idx),]
  result$strand <- Peakstrand
  ## remove really short ones
  result <- result[result[,4]>minlen,]
  ## sort according to score

  ii <- sort(result$score, decreasing = FALSE, index=TRUE)
  result <- result[ii$ix,]
  result
}






defineStrand <- function(strand){
  if(all(strand== "+")){
    peakStrand = "+"
  }else if (all(strand == "-")) {
    peakStrand = "-"
  } else if (sum(strand == "*")>0 |
             (sum(strand == "+") >0 &&
              sum(strand == "-") >0)){
    peakStrand = "*"
  }else if(sum(strand == ".") >0){
    peakStrand = "."
  }
  return(peakStrand)
}


