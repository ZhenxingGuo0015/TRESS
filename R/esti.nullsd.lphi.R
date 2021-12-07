esti.nullsd.lphi <- function(mlphi, nreps, mu, sf.x, sf.y, theta,
                             nsim.use = 10
                             ){

  log.lik.phi <- function(phi, theta, yy, xx, mmu, sx, sy){
    loglik = sum(dnbinom(yy, size = mmu*(phi^{-1} - 1),
                         prob = 1/(1+sy*theta), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(phi^{-1} - 1),
                           prob = 1/(1+sx*theta), log = TRUE)
    )
    loglik
  }
  call.rgamma <- function(x, n){
    rgamma(n, shape = x[seq_len(n)], scale = x[length(x)])
  }
  phi = rep(exp(mlphi), length(theta))
  Phi.esti = matrix(NA, nrow = nrow(mu), ncol = nsim.use)
  for (isi in seq_len(nsim.use)) {
    cat(paste0("The ", isi, "th sampling of data..."), sep = "\n")
    #### simulate X and Y
   # set.seed(isi)
    plx = cbind((1- mu)*(1/phi - 1), theta)
    lambda_x = t(apply(plx, 1, call.rgamma, n = ncol(mu)))
    ply = cbind(mu*(1/phi - 1), theta)
    lambda_y = t(apply(ply, 1, call.rgamma, n = ncol(mu)))
    para.x = sweep(lambda_x, 2, sf.x, FUN = "*")
    para.y = sweep(lambda_y, 2, sf.y, FUN = "*")
    count.Input = matrix(rpois(n = nrow(mu)*ncol(mu), as.vector(para.x)),
                         nrow = nrow(mu), byrow = FALSE)
    count.IP = matrix(rpois(n = nrow(mu)*ncol(mu), as.vector(para.y)),
                      nrow = nrow(mu), byrow = FALSE)
    X.na = sum(is.na(count.Input))
    if(length(X.na) > 0)
      count.Input[is.na(count.Input)] = rpois(X.na, 1e+06)
    Y.na = sum(is.na(count.IP))
    if(length(Y.na) > 0)
      count.IP[is.na(count.IP)] = rpois(Y.na, 1e+06)

    ####### estimate phi
    for (i in seq_len(nrow(mu))) {
      if(!is.na(theta[i])){
        phi.null = optimize(log.lik.phi, c(0, 1),
                            tol = 0.0001,
                            theta = theta[i],
                            yy = count.IP[i, ], xx = count.Input[i, ],
                            mmu = mu[i, ],
                            sx = sf.x, sy = sf.y,
                            maximum = TRUE)$maximum
        Phi.esti[i, isi] = phi.null
      }
    }
  }
  # #######  added on Nov 10, 2021
  S.BAR = mean(colIQRs(log(Phi.esti)/1.349), na.rm = TRUE)
  return(S.BAR)
}
