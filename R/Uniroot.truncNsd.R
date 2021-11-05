Uniroot.truncNsd <- function(Y, a = -2, b = 2){
  #### X ~ N(0, sigma_x^2)
  #### Y = X*I(a<= X <= b) # Y is observed statistics
  ###### The relationship between E(Y) and E(X),
  ######  var(Y) and var(X) can be found through wiki

  ### In fact, based on the relationship,
  ### if E(X) = 0, a = -b, then E(Y) = 0, var(Y) = E(Y^2)
  ### Given E(Y), and var(Y) estimated using mean(Y^2),
  ### we can use root function to estimate sigma_x^2
  myfun <-function(x, c0, a, b){
    alpha = a/x
    beta = b/x
    x^2*( 1+(alpha*dnorm(alpha)-beta*dnorm(beta)
             )/(pnorm(beta)-pnorm(alpha)) ) - c0
  }

  Y1 = Y[which(Y<= b & Y >= a)]
  c0 = mean(Y1^2, na.rm = TRUE)

  # myfun(x = 0.001, c0 = c0, a=-2, b = 2);
  # myfun(x = 10, c0 = c0, a=-2, b = 2)
  res = uniroot(myfun, interval = c(0.001, 100),
                extendInt="yes", c0 = c0,
                a = a,  b = b)
  return(res$root)
}
