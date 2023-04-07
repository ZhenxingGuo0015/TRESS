WaldTest <- function(Coef, Cov, Contrast, nullModel ="standN"){
  # simple Wald test for: Contrast^T%*%beta = 0
  #### nullModel: which null model to use to calculate p-value
  Contrast = as.matrix(Contrast, ncol = 1)
  se = rep(NA, nrow(Coef))
  #for (i in 1:nrow(Coef)) {
  for (i in seq_len(nrow(Coef))) {
    if(!all(is.na(Cov[[i]]))){
      tmp = Cov[[i]]
      se[i] = sqrt( t(Contrast)%*%Cov[[i]]%*% Contrast)
    }
  }
  lor = t(Contrast)%*%t(Coef)
  stat = (t(Contrast)%*%t(Coef))/se
  tmp.infer = DMRInfer(stat, nullModel = nullModel)
  return(res = data.frame(logOR = as.numeric(lor),
                          lorSE = as.numeric(se),
                          stat = round(as.numeric(stat), 3),
                          pvalue = as.numeric(tmp.infer$pvalue),
                          padj = as.numeric(tmp.infer$padj)))
}

