meRatio <- function(counts, sf){
  ### calculate IP/Input ratio for each sample
  ### counts: read counts for input1, ip1, input2, ip2, input3, ip3, ...
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
  ratio = (lambda_y.hat)/(lambda_x.hat + lambda_y.hat)
  ratio[is.na(ratio)] = runif(sum(is.na(ratio)), 0.1, 0.2) #0.001
  ratio[ratio == 0] = runif(sum(ratio == 0), 0.1, 0.2) #0.01
  ratio[ratio == 1] = runif(sum(ratio == 1), 0.9, 0.95) #0.99
  return(ratio)
}
