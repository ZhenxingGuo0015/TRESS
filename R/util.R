## function to read BAM file and create GRanges for the reads
### read in the BAM file and create GRanges for all reads
read.BAM <- function(fn){
  param = ScanBamParam(what=c("rname","strand","pos","qwidth"))
  bam = scanBam(fn, param=param)[[1]]
  ix = !is.na(bam$rname) & !is.na(bam$pos)
  qwidth = bam$qwidth[ix]
  IRange.reads <- GRanges(seqnames=Rle(bam$rname[ix]),
                          ranges=IRanges(bam$pos[ix],
                                         width=bam$qwidth[ix]),
                          strand=Rle(bam$strand[ix]))
  IRange.reads
}


## get read counts for given genomic windows
getWinCounts <- function(files, wins,
                         filetype = c("bed", "bam", "GRanges")) {
  if((!is.data.frame(wins)) & is(wins)[1]!="GRanges")
    stop("Input genomic intervals must be a GRanges or data frame!")
  if(is.data.frame(wins)) wins = import(wins)
  counts = matrix(0, nrow = length(wins), ncol = length(files))
  #for(i in 1:length(files)) {
  for(i in seq_len(length(files))) {
    if(filetype!= "GRanges"){
      if(filetype == "bam")
        reads = read.BAM(files[i])
      else if(filetype == "bed")
        reads = import(files[i])
    }else if(filetype == "GRanges"){
      ### data in datasetTRES package
      reads = get(files[i])
    }
    ### added my Zhenxing to avoid repeated counting
    # width(reads) = 2
    #### ended
    counts[,i] = countOverlaps(wins,reads)
  }
  colnames(counts) = files
  counts
}



exonBins.byTXDB <- function(txdb,
                            binsize = 50,
                            IncludeIntron = FALSE){
  #### each transcript will be divided into bins of length "binsize"
  ### when the annotation file are provided by
  ###    user instead of downloading from internet
  tmp = unique(transcripts(txdb))
  # idx <- split(1:length(as.data.frame(tmp)$seqnames),
  #              as.data.frame(tmp)$seqnames)

  idx <- split(seq_len(length(as.data.frame(tmp)$seqnames)),
               as.data.frame(tmp)$seqnames)
  iii = (grepl("random", names(idx)) |
           grepl("M", names(idx)) |
           grepl("hap", names(idx)) |
           grepl("_", names(idx)))
  idx[iii] = NULL

  allchr = names(idx)
  ALL = vector("list", length = length(allchr))
 # for (ichr in 1:length(allchr)) {
  for (ichr in seq_len(length(allchr))) {
    cat(allchr[ichr], sep = "\t")
    thischr = allchr[ichr]
    thiset  = idx[[ichr]]
    thiset.strand = as.character(strand(tmp[thiset]))
    ### split negative and positive bins
    ### positive
    ii.pos = which(thiset.strand == "+")
    pos.min = min(start(tmp[thiset[ii.pos]]))
    pos.max = max(end(tmp[thiset[ii.pos]]))

    posbins.start = seq(pos.min, pos.max, by = binsize)
    posbins.end = c(posbins.start[-1] - 1,pos.max)
    posbin.strand = unique(tmp@strand[thiset[ii.pos]])
    pos.bin = GRanges(Rle(thischr),
                      IRanges(posbins.start, posbins.end),
                      Rle("+"))
    ### negative
    ii.neg = which(thiset.strand == "-")
    neg.min = min(start(tmp[thiset[ii.neg]]))
    neg.max = max(end(tmp[thiset[ii.neg]]))

    negbins.start = seq(neg.min, neg.max, by = binsize)
    negbins.end = c(negbins.start[-1] - 1, neg.max)
    negbin.strand = unique(tmp@strand[thiset[ii.neg]])
    neg.bin = GRanges(Rle(thischr),
                      IRanges(negbins.start, negbins.end),
                      Rle("-"))

    #### combine positive and negtive bins
    thisbins = c(pos.bin, neg.bin)
    ALL[[ichr]] = thisbins
  }

  allbins = unlist(GRangesList(ALL))


  ####
 # require(GenomicFeatures)
  allExons = unlist(exonsBy(txdb))
  allIntron = unlist(intronsByTranscript(txdb))

  if(!IncludeIntron){
    ## get bins on the exons
    #ix = allbins %over% allExons
    ix = (countOverlaps(allbins, allExons) > 0 &
            width(allbins) == 50)
    bins.exons = allbins[ix]
    ## return
    invisible(return(list(bins=bins.exons,
                          keep.bin = ix,
                          allExons=allExons)))
  }else{
    #ix = (allbins %over% allExons) | (allbins %over% allIntron)

    ix = ((countOverlaps(allbins, allExons)+
             countOverlaps(allbins, allIntron)) > 0 &
           width(allbins) == 50)
      #(allbins %over% allExons) | (allbins %over% allIntron)

    bins.exIn = allbins[ix]
    invisible(return(list(bins=bins.exIn,
                          keep.bin = ix,
                          allExons=allExons,
                          allIntrons = allIntron)))
  }

}




exonBins.byTXDB.2 <- function(txdb,
                              binsize = 50,
                              IncludeIntron = FALSE){
  #### each transcript will be divided into bins of length "binsize"
  ### when the annotation file are provided by
  ###    user instead of downloading from internet
  tmp = suppressMessages(unique(genes(txdb)))
  idx = (grepl("random", tmp@seqnames)|
           grepl("M", tmp@seqnames)|
           grepl("hap", tmp@seqnames)|
           grepl("_", tmp@seqnames))
  allgenes = tmp[!idx]
  unique.gene = union(allgenes, allgenes) # combine overlapped gene
 # tile.bin = tile(unique.gene, n = floor(width(unique.gene)/(binsize)))
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



getStrand <- function(bins, anno_TXDB){
  ### get bin strand given bin sites and genome
  # anno_TXDB: an annotation file in form of txdb
  txdb = anno_TXDB
  allGenes = genes(txdb)
  allGenes.strand = as.character(strand(allGenes))

  tmp = findOverlaps(bins, allGenes)
  length(tmp) ## there are bins overlapping more than one gene.
  binIdx = queryHits(tmp)
  geneIdx = subjectHits(tmp)

  tt = table(binIdx)
  binNames = as.integer(names(tt))

  ### start to get strands
  binStrand = rep(".", length(bins))

  ## unique ones
  ix.unique = binNames[which(tt==1)]
  tmpidx = which(binIdx %in% ix.unique)
  geneidx.unique = geneIdx[tmpidx]
  binStrand[ix.unique] = allGenes.strand[geneidx.unique]

  ## deal with bins overlapping multiple genes
  ix.dup = binNames[which(tt>1)]
  tmpidx = which(binIdx %in% ix.dup)
  geneidx.dup = geneIdx[tmpidx]
  strands.dup = allGenes.strand[geneidx.dup]
  binIdx.dup = binIdx[tmpidx]
  res <- tapply(strands.dup, binIdx.dup,
                function(x) {
                  if(all(x == "+")) return("+")
                  else if(all(x == "-")) return("-")
                  else return("*")
                })
  binStrand[ix.dup] = res
  return(binStrand)
}

