TRESS_DMRfit <- function(IP.file, Input.file,
                         Path_To_AnnoSqlite,
                         variable = NULL,
                         model = NULL,
                         InputDir,
                         OutputDir = NA,
                         experimentName = NA,
                         binsize = 50,
                         filetype = "bam",
                         IncludeIntron = FALSE,
                         filterRegion = TRUE,
                         shrkPhi = TRUE,
                         addsuedo = FALSE
                         ){
  #### A wrapper function for differential peak calling

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
    stop("Must provide the annotation file for TRES to work!",
         call. = TRUE, domain = NULL)
  }

  if(length(variable)==0 | length(model) == 0){
    stop("Must provide variable and model for model fitting!",
         call. = TRUE, domain = NULL)
  }
  if( !is.na(OutputDir) & is.na(experimentName)){
    stop("Must provide a name for your comparison to save results!",
         call. = TRUE, domain = NULL)
  }
  t.1 = Sys.time()
  cat("##### Divid the genome into bins and obtain bin counts...",
      sep = "\n")
  t1 = Sys.time()
  allBins = DivideBins(IP.file = IP.file,
                       Input.file = Input.file,
                       Path_To_AnnoSqlite = Path_To_AnnoSqlite,
                       InputDir = InputDir,
                       OutputDir = OutputDir,
                       experimentName = experimentName,
                       binsize = binsize,
                       filetype = filetype,
                       IncludeIntron = IncludeIntron)
  t2 = Sys.time()
  cat("Time used to obtain bin-level data is: ", t2 - t1, sep = "\n")
  cat("##### Step 1: Call candidate DMRs...", sep = "\n")
  t1 = Sys.time()
  Candidates = CallCandidates(Counts = allBins$binCount,
                              bins = allBins$bins)
  Candidates$lg.fc = NULL   ## not necessary for DMR calling
  if(filterRegion){
    Candidates = filterRegions(Candidates)
  }
  t2 = Sys.time()
  cat("The number of candidates is: ",nrow(Candidates$Regions),
      sep = "\n")
  cat("Time used in Step 1 is: ", t2 - t1, sep = "\n")
  if(!is.na(OutputDir)){
    save(Candidates,
         file = paste0(OutputDir, "/",
                       experimentName,"_Candidates.rda"))
  }

  cat("##### Step 2: Model fitting on candidates...", sep = "\n")
  ### model and parameter estimation for candidate DMRs
  t1 = Sys.time()
  DMRfit = CallDMRs.paramEsti(counts = Candidates$Counts,
                              sf = Candidates$sf,
                              variable = variable,
                              model = model,
                              shrkPhi = shrkPhi,
                              addsuedo = addsuedo)
  t2 = Sys.time()
  cat("Time used in Step 2 is: ", t2 - t1, sep = "\n")

  t.2 = Sys.time()
  cat("##### Done! Total time used is : ",
      t.2 - t.1, sep = "\n")
  return(DMRfit)
}
