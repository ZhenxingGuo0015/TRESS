Analyzing MeRIP-seq data with TRESS
=================

`TRESS` is an R package desinged for the RNA methylation sequencing data analysis. 

The post-transcriptional epigenetic modiﬁcation on mRNA is an emerging ﬁeld to study the 
gene regulatory mechanism and their association with diseases. 
Recently developed high-throughput sequencing technology named Methylated RNA Immunoprecipitation Sequencing (MeRIP-seq) 
enables one to proﬁle mRNA epigenetic modiﬁcation transcriptome-wide. Two major tasks in the analysis of MeRIP-seq 
data is to identify transcriptome-wide m6A regions (namely "peak calling") and differential m6A regions (differential peak calling). 

Our package TRESS provides functions for peak calling of MeRIP-seq data, 
based on empirical Bayesian hierarchical models. 
The method accounts for various sources of variations in the data through rigorous modeling, 
and achieves shrinkage estimation by 
borrowing information from transcriptome-wide data to stabilize the parameter estimation.

Here, we briefly describe how to install TRESS package through GitHub. For detailed usage of TRESS, 
please refer to the vignette file.



## Installation
From GitHub: 

```r
install.packages("devtools") # if you have not installed "devtools" package
library(devtools)
install_github("https://github.com/ZhenxingGuo0015/TRESS", build_vignettes = TRUE)
```

To view the package vignette in HTML format, run the following lines in R

```r
library(TRESS)
browseVignettes("TRESS")
```

## Quick start on peak calling
Here we provide quick examples of how TRESS performs peak 
calling and differential peak calling.
Prior to analysis, TRESS requires paired 
input control and IP BAM files for each replicate of all samples: 
"input1.bam \& ip1.bam", "input2.bam \& ip2.bam", .... 
The BAM files contain mapped reads sequenced from 
respective samples and are the output of sequence alignment tools 
like ``Bowtie2``. In addition to BAM files, 
TRESS also needs the genome annotation of reads saved 
in format of ``*.sqlite``.

For illustration purpose, we include four example BAM files 
and one corresponding genome annotation file in 
our publicly available data package ``datasetTRES``on github, 
which can be installed with
```{r, eval= FALSE}
install_github("https://github.com/ZhenxingGuo0015/datasetTRES")
```
The BAM files contain sequencing reads (only on chromosome 19) 
from two input \& IP mouse brain cerebellum samples.
Given both BAM and annotation files, 
peak calling in TRESS is conducted 
by:

```{r, eval= FALSE}
## Directly take BAM files in "datasetTRES" available on github
library(TRESS)
library(datasetTRES)
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
annoDir = file.path(system.file(package = "datasetTRES"),
                    "extdata/mm9_chr19_knownGene.sqlite")
OutDir = "/directory/to/output"  
TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = BamDir,
           OutputDir = OutDir, # specify a directory for output
           experiment_name = "examplebyBam", # name your output 
           filetype = "bam")
```
```{r, eval= TRUE}
### example peaks
peaks = read.table(file.path(system.file(package = "TRESS"),
                           "extdata/examplebyBam_peaks.xls"),
                 sep = "\t", header = TRUE)
head(peaks[, -c(5, 14, 15)], 3)
```

To replace the example BAM files with your BAM files, the codes are:
```{r, eval=FALSE}
## or, take BAM files from your path
Input.file = c("input_rep1.bam", "input_rep2.bam")
IP.file = c("ip_rep1.bam", "ip_rep2.bam")
BamDir = "/directory/to/BAMfile"
annoDir = "/path/to/xxx.sqlite"
OutDir = "/directory/to/output"
TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = BamDir,
           OutputDir = OutDir,
           experiment_name = "example",
           filetype = "bam")
peaks = read.table(paste0(OutDir, "/", 
                          "example_peaks.xls"), 
                   sep = "\t", header = TRUE)
head(peaks, 3)
```

For detailed usage of the package, please refer to the vignette file through

```r
browseVignettes("TRESS")
```
