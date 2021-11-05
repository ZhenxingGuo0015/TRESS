MLE.parallel <- function(mat, sf, D){
  ### Perform MLE estimation for NB model
  ### mat: read count matrix: input1, ip1, input2, ip2, ....
  ### sf: size factor of all samples
  ### D: design matrix
  ### variable: a dataframe containing the predictor
  ### variable where the inference would be conducted for,
  ###           and other covariates (if exists)
  X = mat[,seq(1, ncol(mat), 2)]
  sx = sf[seq(1, ncol(mat), 2)]
  Y = mat[, seq(2, ncol(mat), 2)]
  sy = sf[seq(2, ncol(mat), 2)]
  Ratio = meRatio(counts = mat, sf = sf)
  ###### estimate mu, phi by the MLE either under Beta or NB model
  # res.MLE = mclapply(1:nrow(Ratio), iMLE,
  #                    X, Y, sx, sy,
  #                    Ratio, D, iniby,
  #                    mc.set.seed = TRUE,
  #                    mc.cores = 2 #max(1, detectCores() - 1)
  #                    )
  res.MLE = mclapply(seq_len(nrow(Ratio)), iMLE,
                     X, Y, sx, sy,
                     Ratio, D,
                     mc.set.seed = TRUE,
                     mc.cores = 2
                     )
  #####
  res.MLE = matrix(unlist(res.MLE), nrow = length(res.MLE), byrow = TRUE)
 # R.ini = res.MLE[, 1:ncol(D)];
  R.ini = as.matrix(res.MLE[, seq_len(ncol(D))]);
  phi.ini = res.MLE[, (ncol(D)+1)];
  theta.ini = res.MLE[, (ncol(D)+2)]
  R.nb = as.matrix(res.MLE[, ((ncol(D)+2) + 1):((ncol(D)+2) + ncol(D))]);
  phi.nb = res.MLE[, (ncol(D)+2) + ncol(D) + 1];
  theta.nb = res.MLE[, (ncol(D)+2) + ncol(D) + 2]

  #### added on August 08, 2021
  loglik.nb = res.MLE[, 2*(ncol(D)+2) + 1]
  ####

  colnames(R.nb) = colnames(R.ini) = colnames(D)
  res = list(Coef.nb = R.nb, phi.nb = phi.nb,
             theta.nb = theta.nb, loglik.nb = loglik.nb,
             Coef.ini = R.ini, phi.ini = phi.ini, theta.ini = theta.ini)
  return(res)
}
