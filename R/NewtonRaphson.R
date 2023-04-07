############### Functions related to Newton-Raphson based on NB
log.lik_NB <- function(x, y, sx, sy, D, R, s, t){
  ### log likelihood for a single region
  ###### To make sure phi ~ (0,1), theta >0, we reparameterize
  ###### phi = exp(s)/(1+exp(s)), theta = exp(t)
  # mu = exp(D%*%R)/(1 + exp(D%*%R))
  # mu[mu == 1] = 0.9
  # mu[mu == 0] = 0.01
  # phi = max(min(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
  # theta = max(min(exp(t), 1000), 0.0001)
  ##### modified on Aug. 08, 2021
  # sum(dnbinom(x, size = (1-mu)*(phi^{-1} -1),
  #             prob = 1/(1+theta*sx), log = TRUE) +
  #       dnbinom(y, size = mu*(phi^{-1} -1),
  #               prob = 1/(1+theta*sy), log = TRUE))

  if(length(s) == 1){ # a single region
    mu = exp(D%*%R)/(1 + exp(D%*%R))
    mu[mu == 1] = 0.9
    mu[mu == 0] = 0.01
    phi = max(min(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
    theta = max(min(exp(t), 1000), 0.0001)
    loglik = sum(dnbinom(x, size = (1-mu)*(phi^{-1} -1),
                         prob = 1/(1+theta*sx), log = TRUE) +
                   dnbinom(y, size = mu*(phi^{-1} -1),
                           prob = 1/(1+theta*sy), log = TRUE))

  }else if(length(s) > 1){ ### all regions
    mu = t(exp(D%*%t(R))/(1+exp(D%*%t(R))))
    mu[mu == 1] = 0.9
    mu[mu == 0] = 0.01
    phi = pmax(pmin(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
    theta = pmax(pmin(exp(t), 1000), 0.0001)
    loglik = rep(NA, length(s))
    for (i in seq_len(nrow(x))) {
      loglik[i] = sum(dnbinom(x[i, ], size = (1-mu[i, ])*(phi[i]^{-1} -1),
                              prob = 1/(1+theta[i]*sx), log = TRUE) +
                        dnbinom(y[i, ], size = mu[i,]*(phi[i]^{-1} -1),
                                prob = 1/(1+theta[i]*sy), log = TRUE))

    }
    return(loglik)
  }
}


Grad.NB <- function(x, y, sx, sy, D, R, s, t){
  ## Score of log.lik
  ###### To make sure phi ~ (0,1), theta >0, we reparameterize
  ###### phi = exp(s)/(1+exp(s)), theta = exp(t)

  phi = max(min(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
  theta = max(min(exp(t), 1000), 0.0001)

  U.dL.R = dL.R(x, y, sx, sy, D, R, phi, theta)
  U.dL.phi = dL.phi(x, y, sx, sy, D, R, phi, theta)
  U.dL.theta = dL.theta(x, y, sx, sy, D, R, phi, theta)


  ### derivative of log.lik with respect to s, t
  U.dL.s = U.dL.phi*(phi*(1-phi))
  U.dL.t = U.dL.theta*theta
  ####

  U = c(U.dL.R, U.dL.s, U.dL.t)
  names(U) = c(colnames(D), "phi.s", "theta.t")
  return(U)
}

Hessian.NB <- function(x, y, sx, sy, D, R, s, t){
  ### Hessian matrix of log.lik
  phi = max(min(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
  theta = max(min(exp(t), 1000), 0.0001)

  H.d2L.R2 <- d2L.R2(x, y, sx, sy, D, R, phi, theta)
  H.d2L.Rphi <- d2L.Rphi(x, y, sx, sy, D, R, phi, theta)
  H.d2L.Rtheta <- d2L.Rtheta(x, y, sx, sy, D, R, phi, theta)
  H.d2L.phi2 <- d2L.phi2(x, y, sx, sy, D, R, phi, theta)
  H.d2L.phitheta <- d2L.phitheta(x, y, sx, sy, D, R, phi, theta)
  H.d2L.theta2 <- d2L.theta2(x, y, sx, sy, D, R, phi, theta)

  ####
  H.d2L.Rs <- H.d2L.Rphi*(phi*(1-phi))
  ####

  ####
  H.d2L.Rt <- H.d2L.Rtheta*theta
  ####

  ####
  U.dL.phi = dL.phi(x, y, sx, sy, D, R, phi, theta)
  H.d2L.s2 <- ( H.d2L.phi2*(phi*(1-phi)) + U.dL.phi*(1-2*phi)
                )*(phi*(1-phi))
  ####

  ####
  H.d2L.st <- H.d2L.phitheta*(phi*(1-phi))*theta
  ####

  ####
  U.dL.theta = dL.theta(x, y, sx, sy, D, R, phi, theta)
  H.d2L.t2 <- ( H.d2L.theta2*theta + U.dL.theta )*theta
  ####


  p = length(R)
  H = matrix(NA, ncol = p + 2, nrow = p + 2)
  # H[1:p, 1:p] = H.d2L.R2
  # H[1:p, (p+1):(p+2)] = cbind(H.d2L.Rs, H.d2L.Rt)
  # H[(p+1):(p+2), 1:p] = rbind(H.d2L.Rs, H.d2L.Rt)
  H[seq_len(p), seq_len(p)] = H.d2L.R2
  H[seq_len(p), p+seq_len(2)] = cbind(H.d2L.Rs, H.d2L.Rt)
  H[p+seq_len(2), seq_len(p)] = rbind(H.d2L.Rs, H.d2L.Rt)


  H[p+1, p+1] = H.d2L.s2
  H[p+2, p+2] = H.d2L.t2
  H[p+2, p+1] = H[p+1, p+2] = H.d2L.st
  colnames(H) = rownames(H) = c(colnames(D), "phi.s", "theta.t")
  return(H)
}


E.Hessian.NB <- function(x, y, sx, sy, D, R, s, t, M = 1e+05){
  ### Calculate the expectation of Hessian matrix of log.lik
  phi = max(min(exp(s)/(1+exp(s)), 0.99), exp(-20)/(1+exp(-20)))
  theta = max(min(exp(t), 1000), 0.0001)

  EH.d2L.R2 <- E.d2L.R2(x, y, sx, sy, D, R, phi, theta, M = M)
  EH.d2L.Rphi <- E.d2L.Rphi(x, y, sx, sy, D, R, phi, theta, M = M)
  EH.d2L.Rtheta <- E.d2L.Rtheta(x, y, sx, sy, D, R, phi, theta)
  EH.d2L.phi2 <- E.d2L.phi2(x, y, sx, sy, D, R, phi, theta, M = M)
  EH.d2L.phitheta <- E.d2L.phitheta(x, y, sx, sy, D, R, phi, theta)
  EH.d2L.theta2 <- E.d2L.theta2(x, y, sx, sy, D, R, phi, theta, M = M)

  ####
  EH.d2L.Rs <- EH.d2L.Rphi*(phi*(1-phi))
  ####

  ####
  EH.d2L.Rt <- EH.d2L.Rtheta*theta
  ####

  ####
  #U.dL.phi = dL.phi(x, y, sx, sy, D, R, phi, theta)
   # H.d2L.s2 <- ( H.d2L.phi2*(phi*(1-phi)) +
   #                 U.dL.phi*(1-2*phi) )*(phi*(1-phi))
  EH.d2L.s2 <- EH.d2L.phi2*(phi*(1-phi))*(phi*(1-phi))
  ####

  ####
  EH.d2L.st <- EH.d2L.phitheta*(phi*(1-phi))*theta
  ####

  ####
  #U.dL.theta = dL.theta(x, y, sx, sy, D, R, phi, theta)
  #H.d2L.t2 <- ( H.d2L.theta2*theta + U.dL.theta )*theta
  EH.d2L.t2 <- EH.d2L.theta2*theta*theta
  ####


  p = length(R)
  EH = matrix(NA, ncol = p + 2, nrow = p + 2)
  # EH[1:p, 1:p] = EH.d2L.R2
  # EH[1:p, (p+1):(p+2)] = cbind(EH.d2L.Rs, EH.d2L.Rt)
  # EH[(p+1):(p+2), 1:p] = rbind(EH.d2L.Rs, EH.d2L.Rt)
  EH[seq_len(p), seq_len(p)] = EH.d2L.R2
  EH[seq_len(p), p+seq_len(2)] = cbind(EH.d2L.Rs, EH.d2L.Rt)
  EH[p+seq_len(2), seq_len(p)] = rbind(EH.d2L.Rs, EH.d2L.Rt)

  EH[p+1, p+1] = EH.d2L.s2
  EH[p+2, p+2] = EH.d2L.t2
  EH[p+2, p+1] = EH[p+1, p+2] = EH.d2L.st
  colnames(EH) = rownames(EH) = c(colnames(D), "phi.s", "theta.t")
  return(EH)
}


dL.phi <- function(x, y, sx, sy, D, R, phi, theta){
  ### partial derivative of loglik with respect to phi
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- phi^{-2}*sum( (1-mu)*log(1+theta*sx)
                        ) + phi^{-2}*sum( mu*log(1+theta*sy)
                                          )
  tmp2 <- -phi^{-2}*sum( (1-mu)*digamma((1-mu)*(phi^{-1} -1) + x)
                         )
  tmp3 <- phi^{-2}*sum( (1-mu)*digamma((1-mu)*(phi^{-1} -1))
                        )
  tmp4 <- -phi^{-2}*sum( mu*digamma(mu*(phi^{-1} - 1)+y)
                         )
  tmp5 <- phi^{-2}*sum( mu*digamma(mu*(phi^{-1} - 1))
                        )
  tmp1 + tmp2 + tmp3 + tmp4 + tmp5
}


dL.theta <- function(x, y, sx, sy, D, R, phi, theta){
  ### partial derivative of loglik with respect to theta
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- sum( x/(theta*(1+theta*sx)) + y/(theta*(1+theta*sy))
               )
  tmp2 <- -(phi^{-1} - 1)*sum( (sx*(1-mu))/(1+theta*sx) +
                                 (sy*mu)/(1+theta*sy)
                               )
  tmp1 + tmp2
  }



d2L.phi2 <- function(x, y, sx, sy, D, R, phi, theta){
  #### second order derivative
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- -2*phi^{-3}*sum( (1-mu)*log(1+theta*sx) +
                             mu*log(1+theta*sy)
                           )
  tmp2 <- sum( trigamma((1-mu)*(phi^{-1} - 1) + x)*(phi^{-2}*(1-mu))^2
               ) + 2*phi^{-3}*sum(
                 (1-mu)*digamma((1-mu)*(phi^{-1} - 1) + x)
                 )
  tmp3 <- -sum( trigamma((1-mu)*(phi^{-1} - 1))* ( phi^{-2}*(1-mu) )^2
                ) - 2*phi^{-3}*sum(
                  (1-mu)*digamma((1-mu)*(phi^{-1} - 1))
                  )
  tmp4 <- sum( trigamma(mu*(phi^{-1} - 1) + y)* ( phi^{-2}*mu )^2
               ) + 2*phi^{-3}*sum(
                 mu*digamma(mu*(phi^{-1} - 1) + y)
                 )
  tmp5 <- -sum( trigamma(mu*(phi^{-1} - 1))*(phi^{-2}*mu)^2
                ) - 2*phi^{-3}*sum(
                  mu*digamma(mu*(phi^{-1} - 1))
                  )
  tmp1 + tmp2 + tmp3 + tmp4 + tmp5
  }

E.d2L.phi2 <- function(x, y, sx, sy, D, R, phi, theta,
                       M = 1e+05){
  ### expectation of the second derivative
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- -2*phi^{-3}*sum((1-mu)*log(1+theta*sx) +
                            mu*log(1+theta*sy)
                          )

  fx <- function(thisx, thismu, thissx, phi, theta){
    (trigamma((1-thismu)*(phi^{-1} - 1) + thisx)*(phi^{-2}*(1-thismu) )^2
     + 2*phi^{-3}* (1-thismu)*digamma((1-thismu)*(phi^{-1} - 1) + thisx)
     )*dnbinom(
       thisx, size = (1-thismu)*(phi^{-1} - 1),
       prob =  1/(1+thissx*theta)
     )
    }
  tmp2 = 0
 # for (j in 1:length(x)) {
  for (j in seq_len(length(x))) {
    tmp2 = tmp2 + sum(fx(thisx = 0:M, thismu = mu[j],
                         thissx = sx[j], phi = phi, theta = theta))
  }


  #
  tmp3 <- -sum( trigamma((1-mu)*(phi^{-1} - 1))* ( phi^{-2}*(1-mu) )^2
                + 2*phi^{-3}* (1-mu)*digamma((1-mu)*(phi^{-1} - 1))
  )

  fy <- function(thisy, thismu, thissy, phi, theta){
    (trigamma(thismu*(phi^{-1} - 1) + thisy)* ( phi^{-2}*thismu )^2
     + 2*phi^{-3}* thismu*digamma(thismu*(phi^{-1} - 1) + thisy))*
      dnbinom(thisy, size = thismu*(phi^{-1} - 1),
              prob = 1/(1+thissy*theta))
    }
  tmp4 = 0
  #for (j in 1:length(x)) {
  for (j in seq_len(length(x))) {
    tmp4 = tmp4 + sum(fy(thisy = 0:M, thismu = mu[j],
                         thissy = sy[j], phi = phi,
                         theta = theta))
    }

  ###
  tmp5 <- -sum( trigamma(mu*(phi^{-1} - 1))*(phi^{-2}*mu)^2
                + 2*phi^{-3}* mu*digamma(mu*(phi^{-1} - 1))
                )
  tmp1 + tmp2 + tmp3 + tmp4 + tmp5
  }

d2L.theta2 <- function(x, y, sx, sy, D, R, phi, theta){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- -sum( (x*(1+2*theta*sx))/(theta*(1+theta*sx))^2  +
                  (y*(1+2*theta*sy))/(theta*(1+theta*sy))^2
                )
  tmp2 <- (phi^{-1} -1)*sum( (1-mu)*(sx/(1+theta*sx))^2 +
                               mu*(sy/(1+theta*sy))^2
                             )
  tmp1 + tmp2
}


E.d2L.theta2 <- function(x, y, sx, sy, D, R, phi, theta,
                         M = 1e+05){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  ## tmp1
  # tmp1 <- -sum( (x*(1+2*theta*sx))/(theta*(1+theta*sx))^2  +
  #                 (y*(1+2*theta*sy))/(theta*(1+theta*sy))^2)

  fx <- function(this.x, this.mu, this.sx, phi, theta){
    ((this.x*(1+2*theta*this.sx))/(theta*(1+theta*this.sx))^2)*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fy <- function(this.y, this.mu, this.sy, phi, theta){
    ((this.y*(1+2*theta*this.sy))/(theta*(1+theta*this.sy))^2)*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
    }
  tmp1 = 0
  #for (j in 1:length(x)) {
  for (j in seq_len(length(x))) {
    tmp1 = tmp1 + sum(fx(this.x = 0:M, this.mu = mu[j],
                         this.sx = sx[j], phi = phi, theta = theta)
                      ) + sum(fy(this.y = 0:M, this.mu = mu[j],
                                 this.sy = sy[j], phi = phi, theta = theta)
                      )
    }
  tmp1 = - tmp1

  ### tmp2
  tmp2 <- (phi^{-1} -1)*sum( (1-mu)*(sx/(1+theta*sx))^2 +
                               mu*(sy/(1+theta*sy))^2
                             )
  tmp1 + tmp2

}


d2L.phitheta <- function(x, y, sx, sy, D, R, phi, theta){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  phi^{-2}*sum( (sx*(1-mu))/(1+theta*sx) + (sy*mu)/(1+theta*sy)
                )
  }

E.d2L.phitheta <- function(x, y, sx, sy, D, R, phi, theta){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 0] = 0.01
  mu[mu == 1] = 0.9
  phi^{-2}*sum( (sx*(1-mu))/(1+theta*sx) + (sy*mu)/(1+theta*sy)
                )
  }

dmuj.R <- function(d, R){
  ### partial derivative of mu with respect to coefficient
  #### modified on Mar 17, 2021
  if(sum(d*R) > 700){
    ( 1 + exp(700) )^{-2}*exp(700)*d
  }else{
    ( 1 + exp(sum(d*R)) )^{-2}*exp(sum(d*R))*d
  }
  ####
}

dmu.R <- function(D, R){
  ### derivative of mu with respect to coefficient
  # tmp = (1+exp(D%*%R))^{-2}*exp(D%*%R)
  # diag(as.vector(tmp))%*%D

  ### modified on Mar17,2021
  tmp.s = D%*%R
  tmp.s[tmp.s > 700] = 700
  tmp = (1+exp(tmp.s))^{-2}*exp(tmp.s)
  diag(as.vector(tmp))%*%D
  ###
}

dL.R <- function(x, y, sx, sy, D, R, phi, theta){
  #### partial derivative of log.lik with respect to coefficient
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  (phi^{-1} -1)*colSums(diag(as.vector(
    log(1+theta*sx) - log(1+theta*sy) -
      digamma((1-mu)*(phi^{-1}-1) + x)+
      digamma((1-mu)*(phi^{-1} - 1)) +
      digamma(mu*(phi^{-1} - 1) + y) -
      digamma(mu*(phi^{-1} - 1))
  ) )%*%dmu.R(D, R)
  )
}



Fisher.R <- function(x, y, sx, sy, D, R, phi, theta, M = 1e+05){
  ### fisher information for coefficient given phi and theta

  fx.1 <- function(this.x, this.mu, this.sx, phi, theta){
    (digamma((1-this.mu)*(phi^{-1}-1) + this.x))*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fx.2 <- function(this.x, this.mu, this.sx, phi, theta){
    (digamma((1-this.mu)*(phi^{-1}-1) + this.x))^2*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }

  fy.1 <- function(this.y, this.mu, this.sy, phi, theta){
    (digamma(this.mu*(phi^{-1}-1) + this.y))*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
    }
  fy.2 <- function(this.y, this.mu, this.sy, phi, theta){
    (digamma(this.mu*(phi^{-1}-1) + this.y))^2*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
  }

  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  c = log(1+theta*sx) - log(1+theta*sy) +
    digamma((1-mu)*(phi^{-1} - 1)) - digamma(mu*(phi^{-1} - 1))

  E_ttT = matrix(NA, ncol = length(x), nrow = length(x))
  #for (j in 1:length(x)) {
  for (j in seq_len(length(x))) {
    #for (l in j:length(y)) {
    for (l in setdiff(seq_len(length(y)), seq_len(j-1))) {
      if(l == j){
        E_ttT[j,l] = c[j]^2 + 2*c[j]*(
          sum(fy.1(this.y = 0:M, this.mu = mu[j],
                   this.sy = sy[j], phi = phi, theta = theta)) -
            sum(fx.1(this.x = 0:M, this.mu = mu[j],
                     this.sx = sx[j], phi = phi, theta = theta))
          ) +
          sum(fy.2(this.y = 0:M, this.mu = mu[j],
                   this.sy = sy[j], phi = phi, theta = theta) ) -
          2*sum(fy.1(this.y = 0:M, this.mu = mu[j],
                     this.sy = sy[j], phi = phi, theta = theta))*
          sum(fx.1(this.x = 0:M, this.mu = mu[j],
                   this.sx = sx[j], phi = phi, theta = theta)) +
          sum(fx.2(this.x = 0:M, this.mu = mu[j],
                   this.sx = sx[j], phi = phi, theta = theta))
      }else{
        E_ttT[j,l] = (c[j] + sum(
          fy.1(this.y = 0:M, this.mu = mu[j],
               this.sy = sy[j], phi = phi, theta = theta)
          ) -
            sum(fx.1(this.x = 0:M, this.mu = mu[j],
                     this.sx = sx[j], phi = phi, theta = theta)
                )
          )*(c[l] + sum(fy.1(this.y = 0:M, this.mu = mu[l],
                             this.sy = sy[l], phi = phi, theta = theta)
                        ) -
               sum(fx.1(this.x = 0:M, this.mu = mu[l],
                        this.sx = sx[l], phi = phi, theta = theta)
                   )
             )
        E_ttT[l,j] = E_ttT[j,l]
        }
      }
    }
  W = diag(as.vector((1+exp(D%*%R))^{-2}*exp(D%*%R)))
  fisher.R = (phi^{-1} - 1)^2*t(W%*%D)%*%E_ttT%*%(W%*%D)
  return(fisher.R)
  }




d2L.Rtheta <- function(x, y, sx, sy, D, R, phi, theta){
  #mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  (phi^{-1} - 1)*colSums( diag(as.vector(sx/(1+theta*sx) -
                                           sy/(1+theta*sy))
                               )%*%dmu.R(D, R)
                          )

}

E.d2L.Rtheta <- function(x, y, sx, sy, D, R, phi, theta){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  (phi^{-1} - 1)*colSums( diag(as.vector(sx/(1+theta*sx) -
                                           sy/(1+theta*sy))
                               )%*%dmu.R(D, R)
                          )

}


d2L.Rphi <- function(x, y, sx, sy, D, R, phi, theta){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  tmp1 <- -phi^{-2}*colSums(diag(
    as.vector(
      log(1+theta*sx) -
        log(1+theta*sy) -
        digamma((1-mu)*(phi^{-1}-1) + x)+
        digamma((1-mu)*(phi^{-1} - 1)) +
        digamma(mu*(phi^{-1} - 1) + y) -
        digamma(mu*(phi^{-1} - 1))
      )
    )%*%dmu.R(D, R)
    )
  tmp2 <- (phi^{-1} -1)*colSums(diag(
    as.vector(
      trigamma((1-mu)*(phi^{-1} -1 ) + x)*( phi^{-2}*(1-mu) ) -
        trigamma((1-mu)*(phi^{-1} -1 ))*( phi^{-2}*(1-mu)) -
        trigamma(mu*(phi^{-1} -1 ) + y)*( phi^{-2}*mu) +
        trigamma(mu*(phi^{-1} -1 ))*( phi^{-2}*mu)
      )
    )%*% dmu.R(D, R)
    )

  tmp1+tmp2
}




E.d2L.Rphi <- function(x, y, sx, sy, D, R, phi, theta, M = 1e+05){
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  ## digamma
  fx.1 <- function(this.x, this.mu, this.sx, phi, theta){
    (digamma((1-this.mu)*(phi^{-1}-1) + this.x))*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fy.1 <- function(this.y, this.mu, this.sy, phi, theta){
    (digamma(this.mu*(phi^{-1} - 1) + this.y))*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
    }

  ### trigamma
  fx.2 <- function(this.x, this.mu, this.sx, phi, theta){
    (trigamma((1-this.mu)*(phi^{-1} -1 ) + this.x))*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fy.2 <- function(this.y, this.mu, this.sy, phi, theta){
    (trigamma(this.mu*(phi^{-1} - 1) + this.y))*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
    }

  Ex.1 = Ey.1 =
    Ex.2 = Ey.2 = rep(NA, length(x))
 # for (j in 1:length(x)) {
  for (j in seq_len(length(x))) {
    Ex.1[j] = sum(fx.1(this.x = 0:M, this.mu = mu[j],
                       this.sx = sx[j], phi = phi, theta = theta))
    Ey.1[j] = sum(fy.1(this.y = 0:M, this.mu = mu[j],
                       this.sy = sy[j], phi = phi, theta = theta))
    Ex.2[j] = sum(fx.2(this.x = 0:M, this.mu = mu[j],
                       this.sx = sx[j], phi = phi, theta = theta))
    Ey.2[j] = sum(fy.2(this.y = 0:M, this.mu = mu[j],
                       this.sy = sy[j], phi = phi, theta = theta))
  }
  ####

  tmp1 <- -phi^{-2}*colSums(diag(
    as.vector(log(1+theta*sx) - log(1+theta*sy) -
                Ex.1 +
                digamma((1-mu)*(phi^{-1} - 1)) +
                Ey.1 -
                digamma(mu*(phi^{-1} - 1))
              )
    )%*%dmu.R(D, R)
    )

  tmp2 <- (phi^{-1} -1)*colSums(diag(as.vector(
    Ex.2 *( phi^{-2}*(1-mu)) -
      trigamma((1-mu)*(phi^{-1} -1 ))*( phi^{-2}*(1-mu)) -
      Ey.2*( phi^{-2}*mu) +
      trigamma(mu*(phi^{-1} -1 ))*( phi^{-2}*mu)
    )
    )%*% dmu.R(D, R)
    )
  tmp1+tmp2
}


d2muj.R2 <- function(d, R){
  ### partial derivative of mu with respect to coefficient
  d = matrix(d, nrow = length(R), ncol = 1)
  #### modified on Mar 17, 2021
  if(sum(d*R) < 700 ){
    A = exp(sum(d*R) - 2*log(1+ exp(sum(d*R))))
    B = exp(log(2) + 2*sum(d*R) - 3*log(1+ exp(sum(d*R))))
    (A-B)*(d%*%t(d))
  }else{
    A = exp(700 - 2*log(1+ exp(700)))
    B = exp(log(2) + 2*700 - 3*log(1+ exp(700)))
    (A-B)*(d%*%t(d))
  }

  ####
}

d2L.R2 <- function(x, y, sx, sy, D, R, phi, theta){
  ### second order derivative of log.lik with respect to R
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  tmp1 <- tmp2 <- vector("list", length = nrow(D))
  c1 <- c2 <- rep(NA, nrow(D))
  Tsum <- matrix(0, ncol = ncol(D), nrow = ncol(D))
  #for (j in 1:nrow(D)) {
  for (j in seq_len(nrow(D))) {
    c1[j] = trigamma((1-mu[j])*(phi^{-1} - 1) + x[j]) -
      trigamma((1-mu[j])*(phi^{-1} - 1)) +
      trigamma(mu[j]*(phi^{-1} - 1) + y[j]) -
      trigamma(mu[j]*(phi^{-1} - 1))
    v = dmuj.R(d = D[j, ], R = R)
    tmp1[[j]] = matrix(v, ncol = 1, nrow = length(R))%*%t(v)

    c2[j] = log(1+theta*sx[j]) - log(1+theta*sy[j]) -
      digamma((1-mu[j])*(phi^{-1} - 1) + x[j]) +
      digamma((1-mu[j])*(phi^{-1} - 1)) +
      digamma(mu[j]*(phi^{-1} - 1) + y[j]) -
      digamma(mu[j]*(phi^{-1} - 1))
    tmp2[[j]] = d2muj.R2(d = D[j, ], R = R)


    Tsum <- Tsum +
      (phi^{-1} -1)^2*c1[j]*tmp1[[j]] +
      (phi^{-1} -1)*c2[j]*tmp2[[j]]
  }

  Tsum
}



E.d2L.R2 <- function(x, y, sx, sy, D, R, phi, theta, M = 1e+4){
  ### second order derivative of log.lik with respect to R


  ## digamma
  fx.1 <- function(this.x, this.mu, this.sx, phi, theta){
    (digamma((1-this.mu)*(phi^{-1}-1) + this.x))*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fy.1 <- function(this.y, this.mu, this.sy, phi, theta){
    (digamma(this.mu*(phi^{-1} - 1) + this.y))*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
    }


  ### trigamma
  fx.2 <- function(this.x, this.mu, this.sx, phi, theta){
    (trigamma((1-this.mu)*(phi^{-1} -1 ) + this.x))*
      dnbinom(this.x, size = (1-this.mu)*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sx))
    }
  fy.2 <- function(this.y, this.mu, this.sy, phi, theta){
    (trigamma(this.mu*(phi^{-1} - 1) + this.y))*
      dnbinom(this.y, size = this.mu*(phi^{-1} - 1),
              prob = 1/(1+theta*this.sy))
  }

  ########
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp1 <- tmp2 <- vector("list", length = nrow(D))
  E.c1 <- E.c2 <- rep(NA, nrow(D))
  E.Tsum <- matrix(0, ncol = ncol(D), nrow = ncol(D))
  #for (j in 1:nrow(D)) {
  for (j in seq_len(nrow(D))) {
    E.c1[j] =
      sum(fx.2(this.x = 0:M, this.mu = mu[j],
               this.sx = sx[j], phi = phi, theta = theta)) -
      trigamma((1-mu[j])*(phi^{-1} - 1)) +
      sum(fy.2(this.y = 0:M, this.mu = mu[j],
               this.sy = sy[j], phi = phi, theta = theta)) -
      trigamma(mu[j]*(phi^{-1} - 1))

    v = dmuj.R(d = D[j, ], R = R)
    tmp1[[j]] = matrix(v, ncol = 1, nrow = length(R))%*%t(v)

    E.c2[j] = log(1+theta*sx[j]) - log(1+theta*sy[j]) -
      sum(fx.1(this.x = 0:M, this.mu = mu[j],
               this.sx = sx[j], phi = phi, theta = theta)) +
      digamma((1-mu[j])*(phi^{-1} - 1)) +
      sum(fy.1(this.y = 0:M, this.mu = mu[j],
               this.sy = sy[j], phi = phi, theta = theta)) -
      digamma(mu[j]*(phi^{-1} - 1))

    tmp2[[j]] = d2muj.R2(d = D[j, ], R = R)

    E.Tsum <- E.Tsum +
      (phi^{-1} -1)^2*E.c1[j]*tmp1[[j]] +
      (phi^{-1} -1)*E.c2[j]*tmp2[[j]]
  }

  E.Tsum
}



d.d2L_R2.theta <- function(x, y, sx, sy, D, R, phi, theta){
  ### derivative of d2L.R2 with respect to theta, the result is a matrix
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01
  tmp = (sx/(1+theta*sx) - sy/(1+theta*sy))
  s = 0
  #for (j in 1:ncol(D)) {
  for (j in seq_len(ncol(D))) {
   s = s + tmp[j]*d2muj.R2(d = D[j, ], R = R)
  }
  (phi^{-1} - 1)*s
}

d.d2L_R2.phi <- function(x, y, sx, sy, D, R, phi, theta){
  ### derivative of d2L.R2 with respect to phi, the result is a matrix
  mu = exp(D%*%R)/(1 + exp(D%*%R))
  mu[mu == 1] = 0.9
  mu[mu == 0] = 0.01

  Tsum <- matrix(0, ncol = ncol(D), nrow = ncol(D))
  c1 = c2 = c3 = c4 = rep(NA, ncol(D))
  tmp1 = tmp2 = vector("list", ncol(D))
  #for (j in 1:nrow(D)) {
  for (j in seq_len(nrow(D))) {
    c1[j] = 2*(phi^{-1} - 1)*(-phi^{-2})*(
      trigamma((1-mu[j])*(phi^{-1} - 1) + x[j]) -
        trigamma((1-mu[j])*(phi^{-1} - 1)) +
        trigamma(mu[j]*(phi^{-1} - 1) + y[j]) -
        trigamma(mu[j]*(phi^{-1} - 1)) )

    c2[j] = phi^{-2}*( phi^{-1} - 1)^2*(
      -(1-mu[j])*psigamma((1-mu[j])*(phi^{-1} - 1) + x[j], deriv = 2) +
        (1-mu[j])*psigamma((1-mu[j])*(phi^{-1} - 1), deriv = 2) -
        mu[j]*psigamma(mu[j]*(phi^{-1} - 1) + y[j], deriv = 2) +
        mu[j]*psigamma(mu[j]*(phi^{-1} - 1), deriv = 2)
    )

    v = dmuj.R(d = D[j, ], R = R)
    tmp1[[j]] = matrix(v, ncol = 1, nrow = length(R))%*%t(v)


    c3[j] = (-phi^{-2})*(
      log(1+theta*sx[j]) - log(1+theta*sy[j]) -
        digamma((1-mu[j])*(phi^{-1} - 1) + x[j]) +
        digamma((1-mu[j])*(phi^{-1} - 1)) +
        digamma(mu[j]*(phi^{-1} - 1) + y[j]) -
        digamma(mu[j]*(phi^{-1} - 1))
      )
    c4[j] = phi^{-2}*(phi^{-1} -1)*(
      (1-mu[j])*trigamma((1-mu[j])*(phi^{-1} - 1) + x[j]) -
        (1-mu[j])*trigamma((1-mu[j])*(phi^{-1} - 1) ) -
        mu[j]*trigamma(mu[j]*(phi^{-1} - 1) + y[j])+
        mu[j]*trigamma(mu[j]*(phi^{-1} - 1))
      )
    tmp2[[j]] = d2muj.R2(d = D[j, ], R = R)


    Tsum <- Tsum + (c1[j] + c2[j])*tmp1[[j]] +
      (c3[j] + c4[j])*tmp2[[j]]
  }
  Tsum
}











