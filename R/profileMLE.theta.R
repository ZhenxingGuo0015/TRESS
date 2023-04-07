profileMLE.theta <- function(x, y, sx, sy, mu0, phi0, maximum){
  #### profile MLE of theta given all other parameters
  log.lik.theta <- function(theta, phi,yy, xx, mmu, sx, sy){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*theta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*theta), log = TRUE)
    )
    loglik
  }

  res = optimize(log.lik.theta,
                 c(0, 100000), tol = 0.0001,
                 phi = phi0,
                 yy = y, xx = x,
                 mmu = mu0,
                 sx = sx, sy = sy,
                 maximum = maximum)

  return(res)
}
