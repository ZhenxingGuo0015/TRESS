### conduct inference
getInfer <- function(mu, mu0, var.mu){
  ### formulate statistics
  #require(stats)
  # require(matrixStats)
  #require(base)
  stat = (mu - mu0)/sqrt(var.mu)
  pval = 1 - pnorm(stat)
  fdr = p.adjust(pval, method = "fdr")
  return(list(stats = stat, pval = pval, fdr = fdr))
}
