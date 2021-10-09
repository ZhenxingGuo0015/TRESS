getLogFC <- function(Counts, sf, c0 = NULL, windlen = 5){
  ####  calculate log fold change between IP and input
  ### Counts: read counts matrix in order of input1, ip1, input2, ip2,...
  ### sf: size factor of each sample
  ### c0: speudo counts for each paired ip and input replicate

  if(is.null(c0)){
    c0 = rep(5, ncol(Counts)/2)
  }
  if(nrow(Counts) > 1){
    myMeans = rowMeans
  }else{
    myMeans = mean
  }
  idx.x = seq(1, ncol(Counts), 2)
  idx.y = seq(2, ncol(Counts), 2)
  lfc = log(myMeans((Counts[, idx.y]/sf[idx.y] + c0
  )/(Counts[, idx.x]/sf[idx.x] + c0),
  na.rm = TRUE))
  if(length(lfc) > windlen){
    smooth.lfc <- mySmooth(lfc, windlen = windlen)
  }else{
    smooth.lfc = lfc
  }

  return(smooth.lfc)
}
