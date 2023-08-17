CovR.parallel <- function(mat, sf, model, variable,
                          Coef, phi, theta){
  #### Calculate Cov(R) in parallel
  X = mat[,seq(1, ncol(mat), 2)]
  sx = sf[seq(1, ncol(mat), 2)]
  Y = mat[, seq(2, ncol(mat), 2)]
  sy = sf[seq(2, ncol(mat), 2)]
  D = model.matrix(model, variable)

  if(.Platform$OS.type == "windows" | Sys.info()['sysname'] == "Windows"){
    ## Windows, use single core
    res = vector("list", length = nrow(X))
    for (i in 1:nrow(X)) {
      res[[i]] = iCovR(i, X, Y, sx, sy, D,
                       Coef, phi, theta)
    }
  }else{
    ncores = 2#max(detectCores() - 3, 1)
    res = mclapply(seq_len(nrow(mat)), iCovR,
                   X, Y, sx, sy, D,
                   Coef, phi, theta,
                   mc.cores = ncores)
  }
  ########
  R = matrix(NA, nrow = nrow(mat), ncol = ncol(D))
  Cov = vector("list", length = nrow(mat))
  for (i in seq_len(length(res))) {
    R[i, ] = res[[i]]$R
    Cov[[i]] = res[[i]]$Cov
  }
  colnames(R) = colnames(D)
  #########
  return(list(Coef = R, Cov = Cov))
}
