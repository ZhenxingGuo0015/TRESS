\name{TRESS_peak}
\alias{TRESS_peak}
\title{
Detecting m6A methylation regions from
Methylated RNA Immunoprecipitation Sequencing.
}
\description{
This is a wrapper function to call m6A peaks transcriptome wide. When there
are multiple biological replicates, it
\itemize{
  \item Divides the whole genome to obtain bin-level read counts:
  \code{\link{DivideBins}}
  \item Calls candidate m6A methylation regions:
  \code{\link{CallCandidates}}
  \item Model fitting on candidate peaks based on
  Negative Binomial distribution: \code{\link{CallPeaks.multiRep}}
}
If there is only one replicate, it calls \code{\link{CallPeaks.oneRep}}
to detect m6A methylation regions.
}
\usage{
TRESS_peak(IP.file, Input.file, Path_To_AnnoSqlite,
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
           IncludeIntron = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP.file}{
  A vector of characters containing the name of BAM files
  for all IP samples.
}
\item{Input.file}{
  A vector of characters containing the name of BAM files
  for all input control samples.
}

\item{Path_To_AnnoSqlite}{
A character to specify the path to a "*.sqlite" file used
for genome annotation.
}

\item{binsize}{
  A numerical value to specify the size of window to
  bin the genome and get bin-level read counts.
  Default value is 50.}

\item{WhichThreshold}{
  A character to specify which criterion to select significant bins
  in the first step, and also significant m6A regions in the second step.
  It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": The inference is only based on P-values;
  "fdr": The inference is only based on FDR;
  "lfc": The inference is only based on log fold changes
  between normalized IP and normalized input read counts;
  "pval_lfc": The inference is based on both p-values and log fold changes;
  "fdr_lfc": The inference is based on both FDR and log fold changes.
  Default is "fdr_lfc".
}
  \item{pval.cutoff0}{
  A numerical value to specify a cutoff for p-value.
  Default is 1e-5.
  }
  \item{fdr.cutoff0}{
  A numerical value to specify a cutoff for fdr. Default is 0.05.
  }
 \item{lfc.cutoff0}{
  A numerical value to specify a cutoff for log fold change
  between normalized IP and input read counts.
  Default is 0.7 for fold change of 2.
  }

\item{lowcount}{
An integer to filter out regions with total input counts < lowcount.
Default is 30.
}

\item{InputDir}{
 A character to specify the input directory of all BAM files.
}

\item{OutputDir}{
 A character to specify an output directory save all results.
 Default is NA, which will not save any results.
}

\item{experiment_name}{
 A character to specify the name of results.
}

\item{filetype}{
 A character to specify the format of input data.
 Possible choices are: "bam", "bed" and "GRanges".
 Default is "bam".
}

\item{IncludeIntron}{
A logical value indicating whether to include (TRUE)
intronic regions or not (False). Default is FALSE.
}
}
\details{
TRESS implements a two-step procedure
to conduct peak calling for MeRIP-seq data with
multiple biological replicates.
In the first step, it quickly divide the whole
genome into equal sized bins and loosely indentifies
candidate peak regions using an ad hoc procedure.
In the second step, it detects high confident peaks among
candidate regions
and ranks them with more rigorous statistical modeling
based on an empirical Bayesian hierarchical model.

When there is only one biological replciate,
candidate regions from the above
two-step procedure
will be output as the final list of peaks.
P-values come from
binomial test, which are further adjusted using
Benjamini-Hochberg procedure.
}

\value{
If directory OutputDir is specified, this function will
output two sets of results.
One is saved as ".rda", which contains all bin-level data
(genome coordinates and read counts matrix).
The other one is an ".xls" file, which contains
information of all peaks. The columns of the peak excel files are:
\item{chr}{Chromosome number of each peak.}
\item{start}{The start of genomic position of each peak.}
\item{end}{The end of genomic position of each peak.}
\item{strand}{The strand of each peak.}
\item{summit}{The summit of each peak.}
\item{lg.fc}{Log fold change between normalized IP
and normalized input read counts for each peak.}
\item{pvals}{P-values calculated based on the Wald-test.}
\item{p.adj}{Adjusted p-values using Benjamini-Hochberg procedure.}

If there are multiple replicates,
the excel will also include following columns:
\item{mu}{Estimated methylation level of each peak.}
\item{mu.var}{Estimated variance for
methylation level of each peak}
\item{stats}{Wald test statistics of each peak}
\item{shrkPhi}{The shrinkage estimation of dispersion
for mehtylation levels of each peak.}
\item{shrkTheta}{The shrinkage estimation
for scale parameter theta in the gamma distribution.}
\item{rSocre}{A score defined by TRESS to rank each peak.
The higher the score, the higher the rank would be.}
Note, there are additional columns regardless of the number of replicates. Those columns contain read counts from respective samples and have names "*.bam".
}

\author{
Zhenxing Guo <zhenxing.guo@emory.edu>
}

\references{
Guo, Z., Shafik, A. M., Jin, P., Wu, Z., and Wu, H. (2021) Detecting m6A methylation regions from Methylated RNA Immunoprecipitation Sequencing. Bioinformatics, 1-7. \url{https://doi-org.proxy.library.emory.edu/10.1093/bioinformatics/btab181}
}



\examples{
## Use BAM files in datasetTRES
# install_github("https://github.com/ZhenxingGuo0015/datasetTRES")
\dontrun{
library(datasetTRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
annoDir = file.path(
  system.file(package = "datasetTRES"),
  "extdata/mm9_chr19_knownGene.sqlite"
  )
OutDir = "/directory/to/output"
 TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = BamDir,
           OutputDir = OutDir,
           experiment_name = "examplebyBam",
           filetype = "bam")
peaks = read.table(paste0(OutDir, "/", "c"),
                   sep = "\t", header = TRUE)
}
}
