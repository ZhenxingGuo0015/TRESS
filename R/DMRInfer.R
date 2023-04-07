# conduct inference for DMR calling
DMRInfer <- function(stat, nullModel = "standN"){
  ## calculate pval either by mixture normal, truncated normal
  if(nullModel == "2mix"){
    ## 2-mixed gaussian with mean1=mean2=0, but sd1 !=sd2
    res.EM = mix.2norm.onlysd(Y = stat[!is.na(stat)], pi = 0.1)
    sd0 = res.EM$sd1
    pval = 2*(1- pnorm(abs(stat), mean = 0, sd = sd0) )
    fdr = p.adjust(pval, method = "fdr")
    }else if(nullModel == "trunN"){
      ## truncated normal
      bounds = seq(1.5, 2, 0.1)
      sd0.range = rep(NA, length(bounds))
      for (ib in seq_len(length(bounds))) {
        sd0.range[ib] = Uniroot.truncNsd(Y = stat,
                                         a = -bounds[ib],
                                         b = bounds[ib])
        }
      if(max(sd0.range) - min(sd0.range) >= 0.5){
        sd0 = min(sd0.range)
        }else{
          sd0 = sd0.range[length(sd0.range)]
          }
      pval = 2*(1- pnorm(abs(stat), mean = 0, sd = sd0) )
      fdr = p.adjust(pval , method = "fdr")
      }else if(nullModel == "standN"){
        ## standard normal
        pval = 2*(1-pnorm(abs(stat)))
        fdr = p.adjust(pval , method = "fdr")
      }

  return(res = data.frame(pvalue = as.numeric(pval),
                          padj = as.numeric(fdr)))
}

