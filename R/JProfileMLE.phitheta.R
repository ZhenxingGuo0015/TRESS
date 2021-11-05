JProfileMLE.phitheta <- function(x, y, mu0, sx, sy){
  #### jointly obtain the MLE of phi and theta
  log.lik <- function(para, yy, xx, mmu, sx, sy){
    s = para[1]
    t = para[2]
    loglik = sum(dnbinom(yy, size = mmu*(s^{-1} - 1),
                         prob = 1/(1+sy*t), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(s^{-1} - 1),
                           prob = 1/(1+sx*t), log = TRUE)
    )
    loglik
  }

  res = optim(c(exp(-5), exp(1)),
              log.lik,
              yy = y, xx = x,
              mmu = mu0, sx = sx, sy = sy,
              method = "L-BFGS-B",
              lower = rep(0.00001, 2),
              upper = c(0.999, exp(700)),
              control = list(fnscale = -1),
              hessian = FALSE)

  return(res)
}
