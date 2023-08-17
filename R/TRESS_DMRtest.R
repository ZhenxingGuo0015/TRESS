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

  ########## editted on August 16, 2023
  #res = data.frame(baseMean = baseMean, res)

  if("gene_Symbol" %in% colnames(DMR$Candidates$Regions)){
    res = data.frame(chr = DMR$Candidates$Regions$chr,
                     start = DMR$Candidates$Regions$start,
                     end = DMR$Candidates$Regions$end,
                     strand = DMR$Candidates$Regions$strand,
                     gene_Symbol = DMR$Candidates$Regions$gene_Symbol,
                     baseMean = baseMean,
                     res)
  }else{
    res = data.frame(chr = DMR$Candidates$Regions$chr,
                     start = DMR$Candidates$Regions$start,
                     end = DMR$Candidates$Regions$end,
                     strand = DMR$Candidates$Regions$strand,
                     baseMean = baseMean,
                     res)
  }
  ##########

  return(res)
}






