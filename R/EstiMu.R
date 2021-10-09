### estimate methylation level
EstiMu <- function(counts, sf){
  ### function to estimate mu
  ### counts: expression count, input1, ip1,
  ###                 input2, ip2, input3, ip3, ...
  x = counts[, seq(1, ncol(counts), 2)] # input
  y = counts[, seq(2, ncol(counts), 2)] # ip
  n = ncol(y)
  if(is.null(sf)){
    sf = colSums(counts, na.rm = TRUE)/median(
      colSums(counts, na.rm = TRUE))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]
  lambda_x.hat = sweep(x, 2, sf.x, FUN = "/")
  lambda_y.hat = sweep(y, 2, sf.y, FUN = "/")

  #### version 1
  # p.hat = lambda_y.hat/(lambda_x.hat + lambda_y.hat)
  # mu.hat = rowMeans(p.hat)
  #
  #### version 2
  mu = rowMeans(lambda_y.hat)/rowMeans(
    lambda_x.hat + lambda_y.hat)

  return(mu)
}
