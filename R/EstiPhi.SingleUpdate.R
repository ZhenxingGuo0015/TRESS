### only update phi
EstiPhi.SingleUpdate <- function(y, x,
                                 sy, sx,
                                 mu, theta,
                                 mlphi, sdlphi){
  # require(stats)
  #require(matrixStats)
  ### only update phi while theta is fixed as its moment estimate
  log.lik <- function(phi, yy, xx, sy, sx, mmu,
                      ttheta, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*ttheta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*ttheta), log = TRUE)
    ) + dlnorm(phi, meanlog = mlphi, sdlog = sdlphi,
               log = TRUE)
    loglik
  }

  theta[is.na(theta)] = mean(theta, na.rm = TRUE)
  res = matrix(0, nrow = nrow(y), ncol = 3)
  ix = which(is.na(mu))
  # for (i in 1:length(mu)) {
  for (i in seq_along(mu)) {
    cat(i, sep = "\n")
    if(!is.na(mu[i])){
      tmp = optimize(log.lik, interval = c(0, 1),
                     yy = y[i, ], xx = x[i, ],
                     sy = sy, sx = sx,
                     mmu = mu[i], ttheta = theta[i],
                     mlphi = mlphi, sdlphi = sdlphi,
                     maximum = TRUE, tol = 1e-05)
      res[i, ] = c(tmp$maximum, theta[i], tmp$objective)
    }
    res[ix, ] = NA
  }
  colnames(res) = c("phi", "theta", "obj")
  res = as.data.frame(res)
  #res$convergence = 0
  res$obj = NULL
  return(res)
}
