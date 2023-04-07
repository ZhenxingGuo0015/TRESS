### read BAM file
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
