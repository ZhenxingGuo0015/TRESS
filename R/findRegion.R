## divide the whole genome into consecutive regions,
## which were sequenced
findRegion <- function(chr, pos, sep=1000) {
  pos.diff <- abs(c(as.integer(0), diff(pos)))
  idx.jump <- which(pos.diff>sep)
  regions <- rbind(c(1, idx.jump),
                   c(idx.jump-1, length(pos)))
  regions
}
