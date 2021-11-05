iMLE <- function(i, X, Y, sx, sy, Ratio, D,
                #iniby = "lm",
                 max.iter = 10,  eps = 1e-5){
  #### conduct MLE estimate based on NB model for a single site.
  #Y.norm = round(sweep(Y, 2, sy, FUN = "/"))
  #X.norm = round(sweep(X, 2, sx, FUN = "/"))
  #Total = X.norm + Y.norm
  # R.ini = InitialCoef(x.norm = X.norm[i, ],
  #                     y.norm = Y.norm[i, ],
  #                     ratio = Ratio[i, ],
  #                     D, iniby = iniby)
  #
  R.ini = solve(t(D)%*%D)%*%t(D)%*%mylogit(Ratio[i, ])
  mu0 = exp(D%*%R.ini)/(1+exp(D%*%R.ini))
  if(!all(mu0 <= 0.51)){
    tmp = JProfileMLE.phitheta(x = X[i, ], y = Y[i, ],
                               mu0 = mu0, sx = sx, sy = sy)
    phi.ini = tmp$par[1];  theta.ini = tmp$par[2]
    ################################### re-iterate phi and theta
    theta.tmp = tmp$par[2]; obj.tmp = tmp$value
    for (iter in seq_len(100)) {
      tmp2 = profileMLE.phi(x = X[i, ], y = Y[i, ],
                            sx = sx, sy = sy,
                            mu0 = mu0,
                            theta0 = theta.tmp,
                            maximum = TRUE)
      cat(iter, tmp2$objective, sep = "\n")
      phi.tmp = tmp2$maximum
      if(abs(tmp2$objective - obj.tmp) > eps){
        tmp3 = profileMLE.theta(x = X[i, ], y = Y[i, ],
                                sx = sx, sy = sy,
                                mu0 = mu0,
                                phi0 = phi.tmp,
                                maximum = TRUE)
        theta.tmp = tmp3$maximum
        obj.tmp = tmp3$objective
      }else{
        break
      }
    }
    phi.NB = phi.tmp; theta.NB = theta.tmp
    ############################### update R
    res = profileMLE.coef(x = X[i, ], y = Y[i, ],
                          sx = sx, sy = sy, D = D,
                          R.ini = R.ini,
                          phi0 = phi.NB, theta0 = theta.NB,
                          max.iter = max.iter)
    R.NB = res$R
    loglik.NB = res$loglik
    ####

  }else{
    phi.ini = theta.ini = phi.NB = theta.NB = loglik.NB =  NA
    R.ini = R.NB = rep(NA, ncol(D))
  }
  return(c(R.ini, phi.ini, theta.ini, R.NB, phi.NB, theta.NB, loglik.NB))
}
