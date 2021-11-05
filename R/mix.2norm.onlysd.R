mix.2norm.onlysd <- function(Y, pi = 0.5,
                             nmaxiter=100, TOL=1e-10*length(Y)){
  ### 2-mixed guassin: mu1 = mu2 = 0, but estimate pi, sd1 and sd2
  n = length(Y)
  ### set initial value for sd1 and sd2
  sd1 = 1; sd2 = 2
  p = c(pi, 1-pi)
  f = cbind(dnorm(Y, mean = 0, sd = sd1),
            dnorm(Y, mean = 0, sd = sd2) )
  tmp = sweep(f, 2, p, FUN = "*")
  loglik_obs.old = sum(log(rowSums(tmp)))

  # E-step
  omega = sweep(tmp, 1, rowSums(tmp), FUN = "/")

  # M-step
  pi = sum(omega[, 1])/n
  sd1 = sqrt(sum(omega[, 1]*Y^2)/sum(omega[, 1]))
  sd2 = sqrt(sum(omega[, 2]*Y^2)/sum(omega[, 2]))

  ### likelihood
  p = c(pi, 1-pi)
  f = cbind(dnorm(Y, mean = 0, sd = sd1),
            dnorm(Y, mean = 0, sd = sd2) )
  tmp = sweep(f, 2, p, FUN = "*")
  loglik_obs.new = sum(log(rowSums(tmp)))
  cat( "loglik=", loglik_obs.new, "sd1=",
       sd1, "sd2=",sd2,"\n")

  while (abs(loglik_obs.new - loglik_obs.old) > TOL) {
    loglik_obs.old = loglik_obs.new
    # E-step
    omega = sweep(tmp, 1, rowSums(tmp), FUN = "/")
    # M-step
    pi = sum(omega[, 1])/n
    sd1 = sqrt(sum(omega[, 1]*Y^2)/sum(omega[, 1]))
    sd2 = sqrt(sum(omega[, 2]*Y^2)/sum(omega[, 2]))

    ### likelihood
    p = c(pi, 1-pi)
    f = cbind(dnorm(Y, mean = 0, sd = sd1),
              dnorm(Y, mean = 0, sd = sd2) )
    tmp = sweep(f, 2, p, FUN = "*")
    loglik_obs.new = sum(log(rowSums(tmp)))
    cat( "loglik=", loglik_obs.new,
         "sd1=", sd1, "sd2=",sd2,"\n")
  }

  ## return a list
  list(sd1=sd1, sd2=sd2, pi=pi)
}
