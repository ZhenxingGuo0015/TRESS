Posterior.phi <- function(counts, sf, D, R,
                          phi.mom, theta.mom){

  ### this function is used to get posterior estimate of phi
  x = counts[, seq(1, ncol(counts), 2)] # input
  y = counts[, seq(2, ncol(counts), 2)] # ip
  if(is.null(sf)){
    sf = colSums(counts, na.rm = TRUE)/median(
      colSums(counts, na.rm = TRUE))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]
  mu = t(exp(D%*%t(R))/(1+exp(D%*%t(R))))

  #### get the posterior of phi
  ######### added on Nov 10, 2021
  id = which(phi.mom > 1e-4)
  phi.bar = mean(phi.mom[id], na.rm = TRUE)
  s2phi = var(phi.mom[id], na.rm = TRUE)
  mlphi = log(phi.bar + 0.00001) - (1/2)*log(s2phi/(phi.bar)^2 + 1)
  sdlphi.old = IQR(log(phi.mom[id]))/1.349
  #### estimate nu in the prior of phi
  s0.bar = esti.nullsd.lphi(mlphi = mlphi, mu = mu,
                            sf.x = sf.x, sf.y = sf.y,
                            theta = theta.mom)
  sdlphi = max(sdlphi.old - s0.bar, 0.1)
  ##########

  ### get posterior of phi
  res = PostPhi.NBlogN(y = y, x = x, sy = sf.y, sx = sf.x,
                       D = D, R = R, theta = theta.mom,
                       mlphi = mlphi, sdlphi = sdlphi)

  return(res)
}
