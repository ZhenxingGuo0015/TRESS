\name{CallPeaks.multiRep}
\alias{CallPeaks.multiRep}
\title{
m6A peak calling with multiple replicates.
}
\description{
This function identifies and ranks significant m6A peaks, given candidate regions obtained from multiple paired of input \& IP replicates.}
\usage{
CallPeaks.multiRep(Candidates, mu.cutoff,
                   WhichThreshold = "fdr_lfc",
                   pval.cutoff = 1e-5,
                   fdr.cutoff = 0.05,
                   lfc.cutoff = 0.7)
}
\arguments{
   \item{Candidates}{
  A list containing: genomic coordinates of each
  candidate region,
  read counts and log fold change between IP and input in
  each candidate region.
  It also contains the size factor of each sample.
  }

  \item{mu.cutoff}{
  A constant specifying the background methylation levels.
  This is estimated automatically based on the first
  step of peak calling.
  }
  \item{WhichThreshold}{
  A character specifying a threshold for significant peaks.
  There are three options: "pval" (only use p-values),
  "fdr" (only use FDR), "lfc" (only use log fold change),
  "pval_lfc" (use both p-values and log fold changes) and
  "fdr_lfc" (use FDR and log fold changes).
  Default is "fdr_lfc".}
  \item{pval.cutoff}{
  A constant indicating the cutoff for p-value. Default is 1e-05.
  }
  \item{fdr.cutoff}{
  A constant indicating the cutoff for FDR. Default is 0.05.
  }
  \item{lfc.cutoff}{
  A constant indicating the cutoff for log fold change.
  Default is 0.7 for fold change of 2.
  }
}

\details{
This function first calls \code{\link{CallPeaks.paramEsti}}
to conduct parameter estimation and hypothesis testing for all
candidate m6A regions. Then it filters and ranks
candidate regions using respective criteria to obtain a list of
significant m6A peaks.
}
\value{
The output is a dataframe whose columns are:
\item{chr}{Chromosome number of each peak.}
\item{start}{The start of genomic position of each peak.}
\item{end}{The end of genomic position of each peak.}
\item{strand}{The strand of each peak.}
\item{summit}{The summit of each peak.}
\item{lg.fc}{Log fold change between normalized IP and
normalized input read counts.}
\item{mu}{Methylation level of each peak if
there are more than one replicates.}
\item{mu.var}{Estimated variance of
methylation level for each peak, when there are more
than one replicates.}
\item{stats}{Wald test statistics of each peak,
when there are more than one replicate.}
\item{shrkPhi}{Shrinkage estimation of mehtylation dispersion
for each peak, when there are
more than one replicates.}
\item{shrkTheta}{Shrinkage estimation for
scale parameter theta in the gamma distribution,
when there are more than one replicates.}
\item{pvals}{ P-value calculated based on the Wald-test.}
\item{p.adj}{Adjusted p-values using Benjamini-Hochberg procedure.}
\item{rSocre}{A score used to rank each peak.
The higher the score, the higher the rank would be.}
Note, there are additional columns with name "*.bam".
These columns contain the read counts from respective samples.
}


\examples{
### A toy example
data("Basal")
CallPeaks.multiRep(
    Candidates = Basal$Candidates,
    mu.cutoff = 0.5
    )
}
