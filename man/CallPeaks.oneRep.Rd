\name{CallPeaks.oneRep}
\alias{CallPeaks.oneRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
m6A peak calling with only one replicate.
}
\description{
This function conducts peak calling for data when
there is only one biological replicate of input and IP sample.
}
\usage{
CallPeaks.oneRep(Counts, bins, sf = NULL,
                 WhichThreshold = "fdr_lfc",
                 pval.cutoff = 1e-05, fdr.cutoff = 0.05,
                 lfc.cutoff = 0.7, windlen = 5, lowCount = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
  A two-column data matrix containing bin-level read counts
  for both IP and input samples.
}
  \item{sf}{
 A numerical vector containg size factors of both IP and input samples. It can be provided by the user, or automatically estimated using "Counts". Default is NULL.
}
  \item{bins}{
  A dataframe containing the genomic locations (chr, start, end, strand) of each bin.
}
  \item{WhichThreshold}{
   A character specifying a criterion to select significant bins
   in bump finding using an ad hoc algorithm.
   There are five options: "pval" (only use p-values),
   "fdr" (only use FDR), "lfc" (only use log fold change),
   "pval_lfc" (use both p-values and log fold changes) and
   "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".
}
  \item{pval.cutoff}{
A constant indicating a cutoff for p-value.
Default is 1e-05.
}
  \item{fdr.cutoff}{
A constant indicating a cutoff for FDR. Default is 0.05.
}
  \item{lfc.cutoff}{
A constant indicating a cutoff for log fold change.
Default is 0.7 for fold change of 2.
}
  \item{windlen}{
An integer specifying the length of consecutive bins
used in simple moving average smooth of log fold change.
Default is 5.
}
  \item{lowCount}{
  An integer to filter out m6A regions with lower read counts.
  Default is 10.
}
}
\details{
When there is only one replicate, TRESS assigns a p-value
for each bin based on the binomial test. Then it calls candidates
with the same algorithm used when there are
multiple biological replicates.
Binomal tests are performed one more time to
select significant candidates as
final list of peaks.
}
\value{
It returns an excel containing the information for each peak:
\item{chr}{Chromosome number of each peak.}
\item{start}{The start of genomic position of each peak.}
\item{end}{The end of genomic position of each peak.}
\item{strand}{The strand of each peak.}
\item{summit}{The summit of each peak.}
\item{pvals}{P-value for each peak calculated based on binomial test.}
\item{p.adj}{Adjusted p-values using Benjamini-Hochberg procedure.}
\item{lg.fc}{Log fold change between normalized IP
and normalized input read counts.}
Note, there are additional columns with name "*.bam".
These columns contain the read counts from IP and input samples.
}


\examples{
## A toy example
data("Basal")
bincounts = Basal$Bins$Counts[, 1:2]
sf0 = Basal$Bins$sf[1:2]
bins = Basal$Bins$Bins
peaks = CallPeaks.oneRep(Counts = bincounts,
                         sf = sf0, bins = bins)
head(peaks, 3)
}
