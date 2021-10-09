### smooth the value of vec
mySmooth <- function(vec, windlen = 3){
  cum = cumsum(vec)
  n = length(cum)
  a = c(cum[((windlen - 1)/2 + 1):n],
        rep(cum[n], (windlen - 1)/2))
  b = c(rep(0, (windlen - 1)/2 + 1),
        cum[seq_len((n - (windlen - 1)/2 -1))])
  res = (a-b)/windlen
  return(res)
}
