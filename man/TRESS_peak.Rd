\name{TRESS_peak}
\alias{TRESS_peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detecting m$^6$A methylation regions from Methylated RNA Immunoprecipitation Sequencing.
}
\description{
This function implements a two-step procedure to conduct peak calling for MeRIP-seq data.
In the first step, it loosely indentify candidate peak regions using an ad hoc procedure. In the second step, it detects high confident peaks and ranks them with more rigorous statistical modeling based on an empirical Bayesian hierarchical model.
}
\usage{
TRESS_peak(IP.file, Input.file, Path_To_AnnoSqlite, binsize = 50,
sf0 = NULL, WhichThreshold = "fdr_lfc", pval.cutoff0 = 1e-5,
fdr.cutoff0 = 0.05, lfc.cutoff0 = 0.7,lowcount = 30,
InputDir, OutputDir = NA, experiment_name, filetype = "bam", IncludeIntron = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP.file}{
  A vector of characters containing the name of bam files for IP samples.
}
\item{Input.file}{
  A vector of characters containing the name of bam files for input control samples.
}

\item{Path_To_AnnoSqlite}{
A character to specify the path to a "Sqlite" file used for genome annotation.
}

\item{binsize}{
  A numerical value to specify the size of window to bin the genome and get bin-level read counts. Default value is 50.}

  \item{sf0}{
  Numerical vectors to specify the size factors of each sample. Default is NULL.
}

\item{WhichThreshold}{
  A character to specify the name of threshod to select significant bins in the first step, and also significant m6A regions in the second step. It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": The inference is only based on P-values;
  "fdr": The inference is only based on FDR;
  "lfc": The inference is only based on log fold changes between normalized IP and normalized input read counts;
  "pval_lfc": The inference is based on both P-values and log fold changes;
  "fdr_lfc": The inference is based on both FDR and log fold changes.
  Default value of ``WhichThreshold" is "fdr_lfc".
}
  \item{pval.cutoff0}{
  A numerical value to specify the cutoff of p-values. Default is 1e-5.
}
  \item{fdr.cutoff0}{
  A numerical value to specify the cutoff of fdr. Default is 0.05.
}
 \item{lfc.cutoff0}{
  A numerical value to specify the cutoff of log foldchanges between normalized IP and input read counts. Default is 0.7 for fold change of 2.
}

\item{lowcount}{
An iteger to filter out regions with total input counts < ``lowcount". Default is 30.
}

\item{InputDir}{
 A character to specify the input directory of all bam files.
}

\item{OutputDir}{
 A character to specify the output directory of all results. Default is NA, which will not save any results.
}

\item{experiment_name}{
 A character to specify the name of results folder.
}

\item{filetype}{
 A character to specify the format of input data. Possible choices are: ``bam", ``bed" and ``GRanges". Note, ``Granges" only works for the example data saved in the data package ``datasetTRES". Default is "bam".
}

\item{IncludeIntron}{
A logical value indicating whether to include (TRUE) intronic regions or not (False). Default is FALSE.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output of this function contains two results. One is saved as ".rda", which contains all the bins and theirs counts. The other is one ".xls" file, which contains a list of peaks called under the two-step procedure.
The columns of the peak excel files are:
\item{chr}{Chromosome number of each peak.}
\item{start}{The start of genomic position of each peak.}
\item{end}{The end of genomic position of each peak.}
\item{strand}{The strand of genomic position of each peak.}
\item{summit}{The summit position of each peak.}
\item{lg.fc}{The log foldchange between normalized IP and normalized input read counts.}
\item{mu}{The methylation level of each peak if there are more than one replicate.}
\item{mu.var}{The estimated variance of for methylation level of each peak, if there are more than one replicate.}
\item{stats}{The Wald test statistics of each peak, if there are more than one replicate.}
\item{shrkPhi}{The shrinkage estimation of dispersion for mehtylation levels of each peak, if there are more than one replicate.}
\item{shrkTheta}{The shrinkage estimation for scale parameter theta in the gamma distribution, if there are more than one replicate.}
\item{pvals}{ P-value calculated based on the Wald-test.}
\item{p.adj}{Adjusted p-values using FDR.}
\item{rSocre}{A score used to rank each region. The higher the score, the higher the rank would be.}
Note, there are additional columns whose name involve the character ".bam". These columns contain the read counts in respective samples.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Zhenxing Guo \email{zhenxing.guo@emory.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## Use bam files in "datasetTRES"
# install_github("https://github.com/ZhenxingGuo0015/datasetTRES")
library(datasetTRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
annoDir = file.path(system.file(package = "datasetTRES"), "extdata/mm9_chr19_knownGene.sqlite")
 TRESS_peak(IP.file = IP.file,
           Input.file = Input.file,
           Path_To_AnnoSqlite = annoDir,
           InputDir = BamDir,
           experiment_name = "examplebyBam",
           filetype = "bam")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
