BgMethylevel <- function(CandidCounts, BinCounts){
  ### calculate background methylation level
  # CandidCounts: read counts from candidate regions
  # BinCounts: read counts from whole transcriptome-wide bins
  # the sample order in "CandidCounts"  and "BinCounts" must be
  # exactly the same as: Input1, IP1, Input2, IP2, ...
  sf = colSums(BinCounts)/median(colSums(BinCounts))
  bgCount = colSums(BinCounts) - colSums(CandidCounts)
  bg.Input = bgCount[seq(1, length(bgCount), 2)]
  bg.IP = bgCount[seq(2, length(bgCount), 2)]
  bg.mu = mean((bg.IP/sf[seq(2, length(bgCount), 2)]
  )/(bg.IP/sf[seq(2, length(bgCount), 2)]
     + bg.Input/sf[seq(1, length(bgCount), 2)]),
  na.rm = TRUE)
  return(bg.mu)
}
