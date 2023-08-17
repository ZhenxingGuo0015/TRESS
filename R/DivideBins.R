## Obtains bins and corresponding counts
DivideBins <- function(IP.file, Input.file,
                       Path_To_AnnoSqlite,
                       InputDir,
                       OutputDir = NA,
                       experimentName,
                       binsize = 50,
                       filetype = "bam",
                       IncludeIntron = FALSE){
  if(length(IP.file) ==0){
    stop("IP samples are missing!",call. = TRUE, domain = NULL)
    }
  if(length(Input.file) == 0){
    stop("Input samples are missing!",call. = TRUE, domain = NULL)
  }
  if(length(Input.file) != length(IP.file)){
    stop("IP and Input samples are not paired!",
         call. = TRUE, domain = NULL)
  }
  if(is.na(Path_To_AnnoSqlite)){
    stop("Must provide the annotation file for TRES to work!",
         call. = TRUE, domain = NULL)
  }
  TXDB = loadDb(Path_To_AnnoSqlite)
  bins = exonBins.byTXDB(txdb = TXDB, binsize = binsize,
                         IncludeIntron = IncludeIntron)
  if(all(bins$bins@strand == "*"))
    bins$bins@strand = Rle(getStrand(anno_TXDB = TXDB, bins = bins$bins))

  allBins = as.data.frame(bins$bins)
  colnames(allBins)[1] = "chr"
  ## bin-level read counts
 # t1 = Sys.time()
  datafiles = rep(NA, 2*length(IP.file))
  datafiles[seq(1, length(datafiles), 2)] = Input.file
  datafiles[seq(2, length(datafiles), 2)] = IP.file
  #cat("Get bin counts for samples of all groups", sep = "\n")
  if(length(InputDir) == 0){
    allCounts = getWinCounts(files = file.path(paste0(InputDir,
                                                      datafiles)),
                             wins = bins$bins, filetype = filetype)
  }else{
    allCounts = getWinCounts(
      files = file.path(paste0(InputDir, "/", datafiles)),
      wins = bins$bins, filetype = filetype)
  }
 # t2 = Sys.time()
  # cat("Time used for obtaining bin counts in all samples is: ",
  #     t2 - t1, sep = "\n")
  colnames(allCounts) = datafiles
  if (!is.na(OutputDir)) {
    save(bins, allCounts,
         file = paste0(OutputDir, "/",
                       experimentName, "_BinsAndCounts.rda")
         )
    }
  return(res = list(bins = allBins, binCount = allCounts))
}
