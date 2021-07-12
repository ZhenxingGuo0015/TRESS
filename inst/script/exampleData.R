#### here is a simple pipeline of how to generate data saved in "inst/extdata/"

##
### take GRanges data from package "datasetTRES"
library(datasetTRES)
IP.file = c("cb_6wk_ip_rep1_chr19", "cb_6wk_ip_rep2_chr19")
Input.file = c("cb_6wk_input_rep1_chr19", "cb_6wk_input_rep2_chr19")
annoDir = file.path(system.file(package = "datasetTRES"), "extdata/mm9_chr19_knownGene.sqlite")
OutDir = "/Users/zhenxingguo/Documents/research/m6a/packagetest"
TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = NULL,
           OutputDir = OutDir,
           experiment_name = "example",
           filetype = "GRanges")
peaks = read.table(paste0(OutDir, "/","example_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)


## or, directly take bam files in datasetTRES
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
OutDir = "/Users/zhenxingguo/Documents/research/m6a/packagetest"
TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = BamDir,
           OutputDir = OutDir,
           experiment_name = "examplebyBam",
           filetype = "bam")
peaks = read.table(paste0(OutDir, "/examplebyBam_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)
