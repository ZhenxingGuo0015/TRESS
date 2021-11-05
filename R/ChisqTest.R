ChisqTest <- function(Coef, Cov, Contrast){
  ### Performe chisquare test x^2_{nrow(Contrast)}
  ###          for composite hypothesis: h(beta) = 0
  #### Coef: Regression coefficients for all sites
  #### Cov: Covariance matrix set for all sites
  #### Contrast: contrast matrix of interest
  TW = rep(NA, nrow(Coef))
  #for (i in 1:nrow(Coef)) {
  for (i in seq_len(nrow(Coef))) {
    tmp = Cov[[i]]
    h_beta = Contrast%*%as.matrix(Coef[i,], ncol = 1)
    TW[i] = t(h_beta)%*%solve(Contrast%*%Cov[[i]]%*% t(Contrast))%*%h_beta
    }
  pval = 1- pchisq(TW, df = nrow(Contrast))
  fdr = p.adjust(pval, method = "fdr")
  return(res = data.frame(stat = round(TW,3),
                          pvalue = as.numeric(pval),
                          padj = as.numeric(fdr)))
}




