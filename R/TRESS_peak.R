## main function of peak calling in TRESS
TRESS_peak <- function(IP.file, Input.file,
                       Path_To_AnnoSqlite,
                       Path_To_OrgdbSqlite = NA,
                       binsize = 50,
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


  if(is.na(Path_To_AnnoSqlite)){
    stop("Must provide the annotation file for TRESS to work!",
         call. = TRUE, domain = NULL)
  }

  ## 1.divide the genome into bins and get bin counts and strand
  cat("##### Divid the genome into bins and obtain bin counts...",
      sep = "\n")
  Bins.Counts = DivideBins(
    IP.file = IP.file,
    Input.file = Input.file,
    Path_To_AnnoSqlite = Path_To_AnnoSqlite,
    InputDir = InputDir,
    OutputDir = OutputDir,
    experimentName = experiment_name,
    binsize = binsize,
    filetype = filetype,
    IncludeIntron = IncludeIntron)
  t2 = Sys.time()
  cat("Time used to obtain bin-level data is: ",
      t2 - t1, sep = "\n")

  ### 2. peak calling
  if(length(IP.file) > 1){
    cat("##### Step 1: Call candidate peaks... ", sep = "\n")
    # 2.1 obtain candidates and estimate background mu
    t1 = Sys.time()
    Peak.candidates = CallCandidates(
      Counts = Bins.Counts$binCount,
      bins = Bins.Counts$bins,
      WhichThreshold = WhichThreshold,
      pval.cutoff = pval.cutoff0,
      fdr.cutoff = fdr.cutoff0,
      lfc.cutoff = lfc.cutoff0,
      lowcount = lowcount)
    t2 = Sys.time()
    cat("Time used to obtain candidates is: ", t2 - t1, sep = "\n")
    cat("The number of candidates is: ", nrow(Peak.candidates$Regions),
        sep = "\n")
    if(nrow(Peak.candidates$Counts) >= 2){
      ### 2.2 peak calling from candidates
      cat("###### Step 2: Peak calling from candidates...",
          sep = "\n")
      t1 = Sys.time()
      bg.mu = BgMethylevel(CandidCounts = Peak.candidates$Counts,
                           BinCounts = Bins.Counts$binCount)

      Peaks = CallPeaks.multiRep(
        Candidates = Peak.candidates,
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
      Peaks = c(Peak.candidates$Regions,
                Peak.candidates$lg.fc,
                Peak.candidates$Counts)
      cat("Less than 2 candidates!", sep = "\n")
    }
  }else if(length(IP.file) == 1){
    sf = colSums(Bins.Counts$binCount)/median(
      colSums(Bins.Counts$binCount))
    Peaks = CallPeaks.oneRep(Counts = Bins.Counts$binCount,
                             bins = Bins.Counts$bins,
                             sf = sf,
                             WhichThreshold = WhichThreshold,
                             pval.cutoff = pval.cutoff0,
                             fdr.cutoff = fdr.cutoff0,
                             lfc.cutoff = lfc.cutoff0,
                             lowCount = lowcount/2)
    if(!is.null(Peaks$width)){
      Peaks$width = NULL
    }

  }

  #### added on August 18, 2023
  # add gene name to each region
  if(nrow(Peaks) > 1 ){
    if(!is.na(Path_To_OrgdbSqlite)){
      GeneName = getGeneID(RegionList = Peaks,
                           Path_To_AnnoSqlite = Path_To_AnnoSqlite,
                           Path_To_OrgdbSqlite = Path_To_OrgdbSqlite)$gene_Symbol
      Peaks$gene_Symbol = GeneName
    }
  }
  ####

  if(!is.na(OutputDir)) {
    ### save results
    write.table(Peaks,
                file = paste0(OutputDir, "/",
                              experiment_name, "_peaks.xls"),
                row.names = FALSE, sep = "\t")
    }
  ###
  t.1 = Sys.time()
  t.1 - t.0
  cat("##### Done! Total time used to call peaks is: ",
      t.1 - t.0, sep = "\n")
}
