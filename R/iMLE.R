iMLE <- function(i, X, Y, sx, sy, Ratio, D,
                 max.iter = 10,  eps = 1e-5){
  #### conduct MLE estimate based on NB model for a single site.
  X.norm = sweep(X, 2, sx, FUN = "/")
  Y.norm = sweep(Y, 2, sy, FUN = "/")
  tmp.x = sweep(X^2 - X, 2, sx^2, FUN = "/")
  tmp.y = sweep(Y^2 - Y, 2, sy^2, FUN = "/")
  tmp.xy = X.norm*Y.norm
  a = rowMeans(X.norm + Y.norm)
  b = rowMeans( tmp.x + tmp.y + 2*tmp.xy )
  phi.mom = 1- a^2/b;
  phi.mom[phi.mom <= 0] = 0.0001
  theta.mom = a/(phi.mom^{-1} - 1)
  ###
  R.ini = solve(t(D)%*%D)%*%t(D)%*%mylogit(Ratio[i, ])
  mu.ini = exp(D%*%R.ini)/(1+exp(D%*%R.ini))
  mu.ini[mu.ini==1] = 0.99 # added on Aug 28, 2022
  mu.ini[mu.ini==0] = 0.01
  if(!all(mu.ini <= 0.51)){
    R.old = R.ini; phi.old = phi.mom[i]; theta.old = theta.mom[i]
    ############################## first round of iteration
    for (iter in seq_len(2)) {
      #### update R
      res = profileMLE.coef(x = X[i, ], y = Y[i, ],
                            sx = sx, sy = sy, D = D,
                            R.ini = R.old,
                            phi0 = phi.old,
                            theta0 = theta.old,
                            max.iter = max.iter)
      R.old = res$R
      loglik.old = res$loglik
      mu.old = exp(D%*%R.old)/(1+exp(D%*%R.old))
      mu.old[mu.old==1] = 0.99 # added on Aug 28, 2022
      mu.old[mu.old==0] = 0.01

      ##### update Phi and theta
      tmp = JProfileMLE.phitheta(x = X[i, ], y = Y[i, ],
                                 mu0 = mu.old, sx = sx, sy = sy)
      phi.old = tmp$par[1];  theta.old = tmp$par[2]
    }
    ## re-iterate phi and theta
    phi.tmp = phi.old; theta.tmp = theta.old; obj.tmp = tmp$value
    for (iter in seq_len(200)) {
      tmp2 = profileMLE.phi(x = X[i, ], y = Y[i, ],
                            sx = sx, sy = sy,
                            mu0 = mu.old,
                            theta0 = theta.tmp,
                            maximum = TRUE)
     # cat(iter, tmp2$objective, sep = "\n")
      if(abs(tmp2$objective - obj.tmp) > eps){
        phi.tmp = tmp2$maximum
        tmp3 = profileMLE.theta(x = X[i, ], y = Y[i, ],
                                sx = sx, sy = sy,
                                mu0 = mu.old,
                                phi0 = phi.tmp,
                                maximum = TRUE)
        theta.tmp = tmp3$maximum
        obj.tmp = tmp3$objective
      }else{
        break
      }
    }

    phi.ini = phi.mom[i]; theta.ini = theta.mom[i];
    R.NB = R.old; phi.NB = phi.tmp; theta.NB = theta.tmp; loglik.NB = obj.tmp
  }else{
    phi.ini = theta.ini = phi.NB = theta.NB = loglik.NB =  NA
    R.ini = R.NB = rep(NA, ncol(D))
  }
  return(c(R.ini, phi.ini, theta.ini, R.NB, phi.NB, theta.NB, loglik.NB))
}
