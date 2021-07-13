## main function of peak calling in TRESS
TRESS_peak <- function(IP.file, Input.file,
                       Path_To_AnnoSqlite,
                       binsize = 50,
                       sf0 = NULL,
                       WhichThreshold = "fdr_lfc",
                       pval.cutoff0 = 1e-5,
                       fdr.cutoff0 = 0.05,
                       lfc.cutoff0 = 0.7,
                       lowcount = 30,
                       InputDir,
                       OutputDir = NA,
                       experiment_name,
                       filetype = "bam",
                       IncludeIntron = FALSE){
  ## 1.divide the genome into bins and get bin counts and strand
  t.0 = Sys.time()

  t1 = Sys.time()
  if(length(IP.file) ==0){
    stop("IP samples are missing!",
         call. = TRUE, domain = NULL)
    }

  if(length(Input.file) == 0){
    stop("Input samples are missing!",
         call. = TRUE, domain = NULL)
    }

  if(length(Input.file) != length(IP.file)){
    stop("IP and Input samples are not paired!",
         call. = TRUE, domain = NULL)
    }

  cat(paste0("Divid the genome into bins..." ), sep = "\n")

  if(is.na(Path_To_AnnoSqlite)){
    stop("Must provide the annotation file for TRES to work!",
         call. = TRUE, domain = NULL)
    }

  TXDB = loadDb(Path_To_AnnoSqlite)
  bins = exonBins.byTXDB(txdb = TXDB, binsize = binsize,
                         IncludeIntron = IncludeIntron)

  if(all(bins$bins@strand == "*")){
    binStrand = getStrand(anno_TXDB = TXDB, bins = bins$bins)
    }else{
      binStrand = bins$bins@strand
      }

  t2 = Sys.time()
  t2 - t1
  cat("Time used for dividing genome is: ", t2 - t1, sep = "\n")

  datafiles = rep(NA, 2*length(IP.file))
  datafiles[seq(1, length(datafiles), 2)] = Input.file
  datafiles[seq(2, length(datafiles), 2)] = IP.file

  cat("Get bin counts...", sep = "\n")
  t1 =  Sys.time()

  if(length(InputDir) == 0){
    allCounts = getWinCounts(
      files = file.path(paste0(InputDir,datafiles)),
      wins = bins$bins, filetype = filetype
      )
    }else{
      allCounts = getWinCounts(
        files = file.path(paste0(InputDir, "/", datafiles)),
        wins = bins$bins, filetype = filetype)
      }

  t2 =  Sys.time()
  t2 - t1
  cat("Time used for obtaining bin counts is: ",
      t2 - t1, sep = "\n")

  ### 2. peak calling
  if(is.null(sf0))
    sf0 = colSums(allCounts)/median(colSums(allCounts))
  allBins = as.data.frame(bins$bins)
  colnames(allBins)[1] = "chr"
  allBins$strand = binStrand

  cat("Start to call peaks...", sep = "\n")
  if(length(IP.file) > 1){
    ##### two-step procedure

    ### step 1
    cat("###### Step 1:...", sep = "\n")
    t1 = Sys.time()
    Peak.candidates = M6Apeak.MultiRep.step1(
      Counts = allCounts,
      bins = allBins,
      sf = sf0,
      WhichThreshold = WhichThreshold,
      pval.cutoff = pval.cutoff0,
      fdr.cutoff = fdr.cutoff0,
      lfc.cutoff = lfc.cutoff0,
      lowcount = lowcount
      )

    t2 = Sys.time()
    t2 - t1
    cat("Time used in Step 1 is: ", t2 - t1, sep = "\n")

    if(nrow(Peak.candidates) >= 2){
      ###  step 2
      cat("###### Step 2:...", sep = "\n")
      t1 = Sys.time()

      ### estimate background methylation level bg.mu
      idx = which(grepl("rep", colnames(Peak.candidates)) |
                    grepl("bam", colnames(Peak.candidates)))
      PeakCount = Peak.candidates[, idx]
      bgCount = colSums(allCounts) - colSums(PeakCount)
      bg.Input = bgCount[seq(1, length(bgCount), 2)]
      bg.IP = bgCount[seq(2, length(bgCount), 2)]
      bg.mu = mean((bg.IP/sf0[seq(2, length(bgCount), 2)]
                    )/(bg.IP/sf0[seq(2, length(bgCount), 2)]
                       + bg.Input/sf0[seq(1, length(bgCount), 2)]),
      na.rm = TRUE)

      #### stop editing
      Peaks = M6Apeak.MultiRep.step2(
        Candidates = Peak.candidates,
        sf = sf0,
        mu.cutoff = bg.mu,
        WhichThreshold = WhichThreshold,
        pval.cutoff = pval.cutoff0,
        fdr.cutoff = fdr.cutoff0,
        lfc.cutoff = lfc.cutoff0
        )
      t2 = Sys.time()
      t2 - t1
      cat("Time used in Step 2 is: ", t2 - t1, sep = "\n")
    }else{
      Peaks = Peak.candidates
      cat("Less than 2 candidates!", sep = "\n")
    }

    if(!is.null(Peaks$width)){
      Peaks$width = NULL
    }

    if (!is.na(OutputDir)) {
      ### save results
      save(bins, binStrand, allCounts,
           file = paste0(OutputDir, "/",
                         experiment_name, ".rda"))
      write.table(Peaks,
                  file = paste0(OutputDir, "/",
                                experiment_name, "_peaks.xls"),
                  row.names = FALSE, sep = "\t")
    }

  }else if(length(IP.file) == 1){
    Peaks = M6Apeak.oneRep(Counts = allCounts,
                           bins = allBins,
                           sf = sf0,
                           WhichThreshold = WhichThreshold,
                           pval.cutoff = pval.cutoff0,
                           fdr.cutoff = fdr.cutoff0,
                           lfc.cutoff = lfc.cutoff0,
                           lowCount = 10)
    if(!is.null(Peaks$width)){
      Peaks$width = NULL
    }

    if(!is.na(OutputDir)) {
      ### save results
      save(bins, binStrand, allCounts,
           file = paste0(OutputDir, "/",
                         experiment_name, ".rda"))

      write.table(Peaks,
                  file = paste0(OutputDir, "/",
                                experiment_name, "_peaks.xls"),
                  row.names = FALSE, sep = "\t")
    }

  }
  cat("###### Done!", sep = "\n")
  ###
  t.1 = Sys.time()
  t.1 - t.0
  cat("Total time used for peak calling is: ",
      t.1 - t.0, sep = "\n")
}
