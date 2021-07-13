# Parameter estimation and inference in step2 of peak calling
M6Apeak.paramEsti <- function(mat,
                              sf = NULL,
                              cutoff = NULL,
                              update = "Joint",
                              trans = NULL,
                              optM = "L-BFGS-B",
                              myfscale = -1e+6){
  ### get statistic for each peak
  mu.hat <- EstiMu(counts = mat, sf = sf)
  para <- EstiPhi(counts = mat, sf = sf,
                  update = update, trans = trans,
                  optM = optM, myfscale = myfscale)
  mu.var <- EstiVarofMu(counts = mat, sf = sf,
                        mu = mu.hat, phi = para$phi,
                        theta = para$theta)
  Infer <- getInfer(mu = mu.hat, mu0 = cutoff,
                    var.mu = mu.var)

  ### ranking score, added on April 14, 2020
  rankScore <- (mu.hat - mean(mu.hat,
                              na.rm = TRUE))/sqrt(mu.var)
  ###
  return(data.frame(mu = mu.hat,
                    mu.var = mu.var,
                    stats = Infer$stat,
                    shrkPhi = para$phi,
                    shrkTheta = para$theta,
                    pvals = Infer$pval,
                    p.adj = Infer$fdr,
                    rScore = rankScore)
         )
  }
