### Estimate Phi and theta jointly
EstiPhi.JointUpdate <- function(y, x,
                                sy, sx,
                                mu,
                                mlphi,
                                sdlphi,
                                trans,
                                optM,
                                myfscale,
                                eps = 1e-05){
  ### function to get estimate of phi and theta with optim function
  ### log likelihood
  log.lik <- function(para, yy, xx, mmu, sx, sy,
                      mlphi, sdlphi){
    s = para[1]
    t = para[2]
    loglik = sum(dnbinom(yy, size = mmu*(s^{-1} - 1),
                         prob = 1/(1+sy*t), log = TRUE) +
                   dnbinom(xx, size = (1 - mmu)*(s^{-1} - 1),
                           prob = 1/(1+sx*t), log = TRUE)
    ) + dlnorm(s, meanlog = mlphi, sdlog = sdlphi,
               log = TRUE)
    loglik
  }

  log.lik_1 <- function(para, yy, xx, mmu, sx, sy,
                        mlphi, sdlphi, trans){
    s = para[1]
    t = para[2]
    if(trans == "sin"){
      loglik = sum(dnbinom(yy, size = mmu*(((sin(s)+1)/2)^{-1} - 1),
                           prob = 1/(1+sy*exp(t)), log = TRUE) +
                     dnbinom(xx, size = (1 - mmu)*(((sin(s)+1)/2)^{-1} - 1),
                             prob = 1/(1+sx*exp(t)), log = TRUE)
                   ) + dlnorm((sin(s)+1)/2, meanlog = mlphi, sdlog = sdlphi,
                              log = TRUE)
    }else if(trans == "exp"){
      loglik = sum(dnbinom(yy, size = mmu*((exp(-exp(s)))^{-1} - 1),
                           prob = 1/(1+sy*exp(t)), log = TRUE) +
                     dnbinom(xx, size = (1 - mmu)*((exp(-exp(s)))^{-1} - 1),
                             prob = 1/(1+sx*exp(t)), log = TRUE)
                   ) + dlnorm(exp(-exp(s)), meanlog = mlphi, sdlog = sdlphi,
                              log = TRUE)
    }
    loglik
  }

  ### estimate with optim function
  res = matrix(0, nrow = nrow(y), ncol = 4)
  ix = which(is.na(mu))
  # for(i in 1:length(mu)){
  for(i in seq_along(mu)){
    if(!is.na(mu[i])){
      if(optM != "L-BFGS-B"){
        tmp = optim(c(-0.01, 0.001), log.lik_1,
                    yy = y[i, ], xx = x[i, ],
                    mmu = mu[i], sx = sx, sy = sy,
                    mlphi = mlphi, sdlphi = sdlphi,
                    trans = trans,
                    method = optM,
                    control = list(fnscale = -1),
                    hessian = FALSE)
      }else if(optM == "L-BFGS-B"){
        tmp = optim(c(exp(-6.9), exp(2.3)), log.lik,
                    yy = y[i, ], xx = x[i, ],
                    mmu = mu[i], sx = sx, sy = sy,
                    mlphi = mlphi, sdlphi = sdlphi,
                    method = optM,
                    lower = rep(0.00001, 2),  # rep(0.00001, 2)
                    upper = c(0.999, exp(700)), #c(0.999, exp(700))
                    control = list(fnscale = myfscale),
                    hessian = FALSE)
      }

      res[i, ] = c(tmp$par, tmp$value, tmp$convergence)
    }
  }
  res[ix, ] = NA

  ### grasp results
  if(optM != "L-BFGS-B"){
    if(trans == "sin"){
      return(data.frame(phi = (sin(res[, 1]) + 1)/2,
                        theta = exp(res[, 2]),
                        obj = res[, 3],
                        convergence = res[, 4]))
    }else if(trans == "exp"){
      return(data.frame(phi = exp(-exp(res[, 1])),
                        theta = exp(res[, 2]),
                        obj = res[, 3],
                        convergence = res[, 4]))
    }
  }else if(optM == "L-BFGS-B"){
    return(data.frame(phi = res[, 1],
                      theta = res[, 2],
                      obj = res[, 3],
                      convergence = res[, 4]))
  }


}
