profileMLE.phi <- function(x, y, sx, sy, mu0, theta0, maximum){
  #### profile MLE of Phi given all other parameters
  log.lik.phi <- function(phi, theta, yy, xx, mmu, sx, sy){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*theta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*theta), log = TRUE)
    )
    loglik
  }

  res = optimize(log.lik.phi, c(0, 1), tol = 0.0001,
                 theta = theta0,
                 yy = y, xx = x,
                 mmu = mu0,
                 sx = sx, sy = sy,
                 maximum = maximum)

  return(res)

}
