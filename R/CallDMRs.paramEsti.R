### DMR modeling and estimate
CallDMRs.paramEsti <- function(counts, sf,
                               model, variable,
                               shrkPhi = TRUE,
                               addsuedo = FALSE){
  #### parameter estimation based on NB model
  if(addsuedo) counts = counts + 5  ### to avoid ratio==1
  Ratio = meRatio(counts = counts, sf = sf) ##methylatinon ratio
  #### 1. Preliminary estimate MLE based on NB model
  print("Start to estimate preliminary MLE ...")
  res.MLE = MLE.parallel(mat = as.matrix(counts),
                         sf = sf,
                         D = model.matrix(model, variable)
                         )
  idx = !is.na(res.MLE$phi.nb)
  ### 2. Update phi, theta with their posterior given other parameters
  print("Update phi with its posterior ...")
  D = model.matrix(model, variable)
  mu = t(exp(D%*%t(res.MLE$Coef.nb))/(1+exp(D%*%t(res.MLE$Coef.nb))))
  if(shrkPhi){
    PostPhi = Posterior.phi(counts = counts[idx, ],
                            sf = sf, D = D,
                            R = res.MLE$Coef.nb[idx, ],
                            phi.mom = res.MLE$phi.nb[idx],
                            theta.mom = res.MLE$theta.nb[idx])
    phi = rep(NA, nrow(counts));
    phi[idx] = as.vector(PostPhi$phi)
    theta = rep(NA, nrow(counts));
    theta[idx] = as.vector(PostPhi$theta)
  }else{
    phi = res.MLE$phi.nb
    theta = res.MLE$theta.nb
  }
  ### 3. Calculate Cov(R)
  print("Calculate variance-covariance of coefficients...")
  res.Cov = CovR.parallel(mat = counts, sf = sf,
                          model = model, variable = variable,
                          Coef = res.MLE$Coef.nb,
                          phi = phi, theta = theta)

  ### 4. Calculate log.likelihood, newly added on Aug.08, 2021
  loglik = log.lik_NB(x = as.matrix(counts[,seq(1, ncol(counts), 2)]),
                      y = as.matrix(counts[,seq(2, ncol(counts), 2) ]),
                      sx = sf[seq(1, length(sf), 2)],
                      sy = sf[seq(2, length(sf), 2)],
                      D = D,
                      R = res.MLE$Coef.nb,
                      s = mylogit(phi),
                      t = log(theta))

  #### 5. output all estimate
  print("Output all estimate...")
  DMR = list(Coef = res.Cov$Coef, Cov = res.Cov$Cov,
             Ratio = Ratio, loglik = loglik
             #shrkPhi = round(phi, 5),
             #shrkTheta = round(theta,2),
            # id.rm = which(!idx)
            )
  return(DMR)
}
