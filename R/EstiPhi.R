### Estimate dispersion
EstiPhi <- function(counts, sf,
                    update = c("OnlyPhi", "Iterative", "Joint"),
                    optM = "L-BFGS-B",
                    myfscale = -1e+6,
                    trans = "sin"){
  # require(stats)
  # require(matrixStats)
  x = counts[, seq(1, ncol(counts), 2)]
  y = counts[, seq(2, ncol(counts), 2)]
  n = ncol(y)
  if(is.null(sf)){
    sf = colSums(counts)/median(colSums(counts))
  }
  sf.x = sf[seq(1, ncol(counts), 2)]
  sf.y = sf[seq(2, ncol(counts), 2)]
  lambda_x.hat = sweep(x, 2, sf.x, FUN = "/")
  lambda_y.hat = sweep(y, 2, sf.y, FUN = "/")
  p.hat = lambda_y.hat/(lambda_x.hat + lambda_y.hat)
  mu.hat = rowMeans(p.hat)

  # require(matrixStats)
  phi.mom = ((n-1)/n)*rowVars(p.hat)/(mu.hat*(1 - mu.hat))
  theta.mom = rowMeans(lambda_x.hat + lambda_y.hat)/(phi.mom^{-1} -1)
  phi.bar = mean(phi.mom, na.rm = TRUE)
  s2phi = var(phi.mom, na.rm = TRUE)
  mlphi = log(phi.bar + 0.00001) - (1/2)*log(s2phi/(phi.bar)^2 + 1)
  sdlphi = sqrt(log(s2phi/(phi.bar)^2 + 1))

  ### bayes estimate of phi and theta
  update <- match.arg(update)
  switch(update,
         OnlyPhi = EstiPhi.SingleUpdate(y = y,
                                        x = x,
                                        sy = sf.y,
                                        sx = sf.x,
                                        mu = mu.hat,
                                        theta = theta.mom,
                                        mlphi = mlphi,
                                        sdlphi = sdlphi),
         Iterative = EstiPhi.IterativeUpdate(y = y,
                                             x = x,
                                             sy = sf.y,
                                             sx = sf.x,
                                             mu = mu.hat,
                                             theta0 = mean(
                                               theta.mom[theta.mom < Inf],
                                               na.rm = TRUE),
                                             mlphi = mlphi,
                                             sdlphi = sdlphi),
         Joint = EstiPhi.JointUpdate(y, x,
                                     sy = sf.y,
                                     sx = sf.x,
                                     mu = mu.hat,
                                     mlphi = mlphi,
                                     sdlphi = sdlphi,
                                     optM = optM,
                                     myfscale = myfscale,
                                     trans = trans)
         )
}
