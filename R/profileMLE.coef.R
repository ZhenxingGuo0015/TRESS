profileMLE.coef <- function(x, y, sx, sy, D, R.ini,
                            phi0, theta0,
                            max.iter = 10, eps = 1e-3){
  #### Newton-Raphson to obtain the profile MLE of R
  R.old = R.ini
  loglik.old = log.lik_NB(x = x, y = y, sx = sx, sy = sy, D = D,
                          R = R.old, s = mylogit(phi0),
                          t = log(theta0))

  for (iter in seq_len(max.iter)) {
    cat("Iteration: ", iter, "log.lik: ", loglik.old, "\n")
    Grad.old = dL.R(x = x, y = y, sx = sx, sy = sy, D = D,
                    R = R.old, phi = phi0,
                    theta = theta0)
    Hessian.old = d2L.R2(x = x, y = y, sx = sx, sy = sy, D = D,
                         R = R.old, phi = phi0,
                         theta = theta0)
    rank.m = rankMatrix(Hessian.old)[1]
    if(rank.m < ncol(Hessian.old)){
      tmp.svd = svd(Hessian.old)
      Hessian.old = (tmp.svd$u%*%diag(tmp.svd$d)%*%t(tmp.svd$v) -
                       0.01*tmp.svd$u%*%diag(
                         c(rep(0, rank.m),
                           rep(1, ncol(Hessian.old) - rank.m))
                       )%*%t(tmp.svd$v)
      )
      R.new = R.old - solve(Hessian.old)%*%Grad.old
    }else{
      R.new = R.old - solve(Hessian.old)%*%Grad.old
    }
    ########
    loglik.new = log.lik_NB(x = x, y = y, sx = sx, sy = sy,
                            D = D, R = R.new,
                            s = mylogit(phi0),
                            t = log(theta0))
    if(is.na(loglik.new)| sqrt(sum(Grad.old^2)) < eps |
       loglik.new - loglik.old < 0){
      break
    }else{
      if(loglik.new - loglik.old > 0){
        R.old = R.new
        loglik.old = loglik.new
      }
    }
  }
  R.NB = R.old
  #return(R.NB)

  ##### added on August 08, 2021
  loglik.NB = loglik.old
  #####

  return(list(R = R.NB, loglik = loglik.old))
}

