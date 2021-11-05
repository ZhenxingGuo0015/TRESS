Analyzing MeRIP-seq data with TRESS
=================

`TRESS` is an R package desinged for the RNA methylation sequencing data analysis. 

The post-transcriptional epigenetic modiﬁcation on mRNA is an emerging ﬁeld to study the 
gene regulatory mechanism and their association with diseases. 
Recently developed high-throughput sequencing technology named Methylated RNA Immunoprecipitation Sequencing (MeRIP-seq) 
enables one to proﬁle mRNA epigenetic modiﬁcation transcriptome-wide. Two major tasks in the analysis of MeRIP-seq 
data is to identify transcriptome-wide m6A regions (namely "peak calling") and differential m6A regions (differential peak calling). 

Our package TRESS provides functions for peak calling and differential peak calling of MeRIP-seq data, 
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

## Quick start on differential peak calling
If one has paired input and IP ("input1.bam \& ip1.bam", 
"input2.bam \& ip2.bam", ..., "inputN.bam \& ipN.bam")
BAM files for samples from 
different conditions, then one can apply TRESS to call
differential m6A methylation regions (DMRs). Note that,
the input order of BAM files from 
different conditions should be appropriately 
listed in case that samples from different conditions 
are mistakenly treated as one group.

As TRESS is designed for differential analysis under 
general experimental design, then in addition to BAM and 
genome annotation files, sample 
attributes determined by all factors in study should also be 
provided to construct a design matrix for model fitting.
For this, TRESS requires a dataframe (taken by ``variable``) 
containing, for each factor, the attribute value of 
all samples (the 
order of sample should be exactly the same as BAM files
taken by TRESS).   
A particular model (taken by ``model``) 
determining which factor will be 
included into design matrix should also be provided.



All aforementioned input requirements 
are for model fitting in TRESS.
For hypothesis testing, TRESS requires a contrast of 
coefficients.
The contrast should be in line with the name and order of all
coefficients in the design matrix. 
It can be a vector for 
simple linear relationship detection 
or a matrix for composite relationship detection.


With all required information prepared, do,
```{r, eval=FALSE, message= FALSE, warning= FALSE}
InputDir = "/directory/to/BAMfile"
Input.file = c("input1.bam", "input2.bam",..., "inputN.bam")
IP.file = c("ip1.bam", "ip2.bam", ..., "ipN.bam")
OutputDir = "/directory/to/output"
Path_sqlit = "/path/to/xxx.sqlite"
variable = "YourVariable" # a dataframe containing both
# testing factor and potential covariates, 
# e.g., for two group comparison with balanced samples
# variable = data.frame(Trt = rep(c("Ctrl", "Trt"), each = N/2))
model = "YourModel"     # e.g. model = ~1 + Trt
DMR.fit = TRESS_DMRfit(IP.file = IP.file,
                       Input.file = Input.file,
                       Path_To_AnnoSqlite = Path_sqlit,
                       variable = variable,
                       model = model,
                       InputDir = InputDir,
                       OutputDir = OutputDir,
                       experimentName = "example"
                       )
CoefName(DMR.fit)# show the name of and order of coefficients 
                 # in the design matrix
Contrast = "YourContrast" # e.g., Contrast = c(0, 1)
DMR.test = TRESS_DMRtest(DMR = DMR.fit, contrast = Contrast)
```
As shown above, TRESS separates the model fitting 
(implemented by function ``TRESS_DMRfit()``), which is the most 
computationally heavy part, from the hypothesis testing 
(implemented by function ``TRESS_DMRtest()``). 
Given an experimental design with multiple factors, 
the parameter estimation (model fitting) only 
needs to be performed once, 
and then the hypothesis testing for 
DMR calling can be performed for different factors efficiently. 

For detailed usage of the package, please refer to the vignette file through

```r
browseVignettes("TRESS")
```
