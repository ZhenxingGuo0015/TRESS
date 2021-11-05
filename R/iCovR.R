iCovR <- function(i, X, Y, sx, sy, D, Coef, phi, theta){
  ### calculate Cov(R) for a particular region
  if(!is.na(phi[i])){
    R.old = Coef[i, ]
    Hessian.old = E.d2L.R2(x = X[i, ], y = Y[i, ],
                           sx = sx, sy = sy, D = D,
                           R = R.old, phi = phi[i],
                           theta = theta[i])
    Hessian.old = adjustHessian(Hessian.old)
    Cov.R = solve(-Hessian.old)
    
    ### make sure Cov(R) is at least semi-positive
    if((!all(diag(Cov.R) > 0)) | (!all(eigen(Cov.R)$values >=0))){
      M0 = max(max(X[i, ]), max(Y[i, ]))+10000
      Hessian.old = E.d2L.R2(x = X[i, ], y = Y[i, ],
                             sx = sx, sy = sy, D = D,
                             R = R.old, phi = phi[i],
                             theta = theta[i], M = max(1e+5, M0))
      Hessian.old = adjustHessian(Hessian.old)
      Cov.R = solve(-Hessian.old)
      ### one more time
      if((!all(diag(Cov.R) > 0)) | (!all(eigen(Cov.R)$values >=0))){
        Hessian.old = E.d2L.R2(x = X[i, ], y = Y[i, ],
                               sx = sx, sy = sy, D = D,
                               R = R.old, phi = phi[i],
                               theta = theta[i], M = max(1e+6, M0))
        Hessian.old = adjustHessian(Hessian.old)
        Cov.R = solve(-Hessian.old)
      }
    }
    R = R.old
  }else{
    R = rep(NA, ncol(D))
    Cov.R = matrix(NA, ncol = ncol(D), nrow = ncol(D))
  }
  colnames(Cov.R) = rownames(Cov.R) = colnames(D)
  return(list(R = R, Cov = Cov.R))
}
