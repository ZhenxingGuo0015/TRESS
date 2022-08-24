### filter out regions based on coefficient of variation
filterRegions <- function(Candidates, quant = 0.25){
  Ratio = meRatio(counts = as.matrix(Candidates$Counts),
                  sf = Candidates$sf)
  CV = rowSds(Ratio, na.rm = TRUE)/rowMeans(Ratio, na.rm = TRUE)
  Q.point = as.numeric(quantile(CV, prob = quant))
  idx = which(CV > Q.point)
  counts = Candidates$Counts[idx, ]
  regions = Candidates$Regions[idx, ]
  rownames(regions) = rownames(counts) = 1:nrow(regions)
  Candidates = list(Regions = regions, Counts = counts,
                    sf = Candidates$sf)
  return(Candidates)
}
