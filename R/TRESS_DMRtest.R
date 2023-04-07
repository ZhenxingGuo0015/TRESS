TRESS_DMRtest <- function(DMR, contrast, nullModel = "standN"){
  ### A function to extract results based on test of interest
  if(is.vector(contrast)){
    ### Wald Test
    res = WaldTest(Coef = DMR$Coef, Cov = DMR$Cov,
                   Contrast = contrast, nullModel = nullModel)
  }else if(is.matrix(contrast)){
    ### Chi-square Test
    res = ChisqTest(Coef = DMR$Coef, Cov = DMR$Cov,
                    Contrast = contrast)
  }
  baseMean = rowMeans(DMR$Ratio, na.rm = TRUE)
  res = data.frame(baseMean = baseMean, res)
  return(res)
}





