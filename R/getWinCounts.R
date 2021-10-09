## get read counts for given genomic windows
getWinCounts <- function(files, wins,
                         filetype = c("bed", "bam", "GRanges")
                         ) {
  if((!is.data.frame(wins)) & is(wins)[1]!="GRanges")
    stop("Input genomic intervals must be a GRanges or data frame!")
  if(is.data.frame(wins)) wins = import(wins)
  counts = matrix(0, nrow = length(wins), ncol = length(files))
  for(i in seq_len(length(files))) {
    if(filetype!= "GRanges"){
      if(filetype == "bam")
        reads = read.BAM(files[i])
      else if(filetype == "bed")
        reads = import(files[i])
      }else if(filetype == "GRanges"){
        ### data in datasetTRES package
        reads = get(files[i])
        }
    #### ended
    counts[,i] = countOverlaps(wins,reads)
    }
  colnames(counts) = files
  counts
  }
