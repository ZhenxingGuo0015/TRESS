Posterior.phi <- function(counts, sf, D, R,
                          phi.mom = NULL,
                          theta.mom = NULL){

  ### this function is used to get posterior estimate of phi
  x = counts[, seq(1, ncol(counts), 2)] # input
  y = counts[, seq(2, ncol(counts), 2)] # ip
  if(is.null(sf)){
    sf = colSums(counts, na.rm = TRUE)/median(
      colSums(counts, na.rm = TRUE))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]

  ### get the posterior of phi
  phi.bar = mean(phi.mom, na.rm = TRUE)
  s2phi = var(phi.mom, na.rm = TRUE)
  mlphi = log(phi.bar + 0.00001) - (1/2)*log(s2phi/(phi.bar)^2 + 1)
  sdlphi = sqrt(log(s2phi/(phi.bar)^2 + 1))

  ### get posterior of phi
  res = PostPhi.NBlogN(y = y, x = x, sy = sf.y, sx = sf.x,
                       D = D, R = R, theta = theta.mom,
                       mlphi = mlphi, sdlphi = sdlphi)

  return(res)
}
