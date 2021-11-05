PostPhi.NBlogN <- function(y, x, sy, sx, D, R, theta,
                           mlphi, sdlphi){
  ### only update phi while theta is fixed as its moment estimate
  log.lik <- function(phi, yy, xx, sy, sx, mmu, ttheta, mlphi,
                      sdlphi){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*ttheta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*ttheta), log = TRUE)
    ) + dlnorm(phi, meanlog = mlphi, sdlog = sdlphi,
               log = TRUE)
    loglik
  }
  mu = t(exp(D%*%t(R))/(1+exp(D%*%t(R))))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  theta[is.na(theta)] = mean(theta, na.rm = TRUE)
  res = matrix(NA, nrow = nrow(y), ncol = 3)
  #ix = which(rowSums(is.na(mu)) == ncol(mu))
  #for (i in 1:nrow(mu)) {
  for (i in seq_len(nrow(mu))) {
    #cat(i, sep = "\n")
    if(sum(is.na(mu[i, ])) != ncol(mu)){
      tmp = optimize(log.lik, interval = c(0, 1),
                     yy = y[i, ], xx = x[i, ],
                     sy = sy, sx = sx,
                     mmu = mu[i, ], ttheta = theta[i],
                     mlphi = mlphi, sdlphi = sdlphi,
                     maximum = TRUE, tol = 1e-05)
      res[i, ] = c(tmp$maximum, theta[i], tmp$objective)
    }
  }
  colnames(res) = c("phi", "theta", "obj")
  res = as.data.frame(res)
  res$convergence = 0
  return(res)
}
