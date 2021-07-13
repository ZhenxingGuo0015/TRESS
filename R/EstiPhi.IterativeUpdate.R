EstiPhi.IterativeUpdate <- function(y, x, sy, sx,
                                    mu, theta0,
                                    mlphi, sdlphi,
                                    maxit = 100, eps = 1e-03){
  ### iteratively update phi and theta
  # require(stats)
  # require(matrixStats)
  log.lik.phi <- function(phi, yy, xx, sy, sx, mmu,
                          ttheta, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*ttheta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*ttheta), log = TRUE)
    ) + dlnorm(phi, meanlog = mlphi, sdlog = sdlphi, log = TRUE)
    loglik
  }

  log.lik.theta <- function(theta, yy, xx, sy, sx,
                            mmu, pphi, mlphi, sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(pphi^{-1} - 1),
                         prob = 1/(1+sy*theta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(pphi^{-1} - 1),
                           prob = 1/(1+sx*theta), log = TRUE)
    ) + dlnorm(pphi, meanlog = mlphi, sdlog = sdlphi,
               log = TRUE)
    loglik
  }

  res = matrix(0, nrow = nrow(y), ncol = 4)
  ix = which(is.na(mu))
  # for (i in 1:length(mu)) {
  for (i in seq_along(mu)) {
    cat(i, sep = "\n")
    if(!is.na(mu[i])){
      delta = 1
      theta.old = theta0
      iter = 0
      while(delta > eps && iter <= maxit){
        iter = iter + 1
        ### update phi while fixing theta
        tmp = optimize(log.lik.phi, interval = c(0, 1),
                       yy = y[i, ], xx = x[i, ],
                       sy = sy, sx = sx,
                       mmu = mu[i], ttheta = theta.old,
                       mlphi = mlphi, sdlphi = sdlphi,
                       maximum = TRUE, tol = 1e-05)
        phi.old = tmp$maximum
        obj.old = tmp$objective
        tmp = optimize(log.lik.theta, interval = c(0, exp(700)),
                       yy = y[i, ], xx = x[i, ],
                       sy = sy, sx = sx,
                       mmu = mu[i], pphi = phi.old,
                       mlphi = mlphi, sdlphi = sdlphi,
                       maximum = TRUE, tol = 1e-05)
        theta.old = tmp$maximum
        obj.new = tmp$objective
        delta = obj.new - obj.old
      }
      res[i, ] = c(phi.old, theta.old, obj.new,
                   as.numeric(iter >= maxit))
    }
    res[ix, ] = NA
  }
  colnames(res) = c("phi", "theta", "obj", "convergence")
  res = as.data.frame(res)
  return(res)
}
