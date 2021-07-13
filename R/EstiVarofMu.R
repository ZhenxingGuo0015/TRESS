EstiVarofMu <- function(counts, sf, mu, phi, theta){
  #require(stats)
  # require(matrixStats)
  ### variance estimate of mu.hat
  n = ncol(counts)/2
  if(is.null(sf)){
    sf = colSums(counts)/median(colSums(counts))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]

  tmp = sum(sf.y^{-1}) + n*theta
  mu.var = (mu/(n^2*theta))*(phi/(1 -  phi))*tmp
  return(mu.var)
}
