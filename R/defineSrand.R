### obtain strand in bump-finding
defineStrand <- function(strand){
  if(all(strand== "+")){
    peakStrand = "+"
  }else if (all(strand == "-")) {
    peakStrand = "-"
  } else if (sum(strand == "*")>0 |
             (sum(strand == "+") >0 &&
              sum(strand == "-") >0)){
    peakStrand = "*"
  }else if(sum(strand == ".") >0){
    peakStrand = "."
  }
  return(peakStrand)
}