\name{M6Apeak.MultiRep.step2}
\alias{M6Apeak.MultiRep.step2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dectect and rank m6A candiate regions.
}
\description{
This is a wrap function, which calls another function
to identify and rank significant m6A regions
using Bayesian hierarchical negative binomial models.
}
\usage{
M6Apeak.MultiRep.step2(Candidates, sf, mu.cutoff,
WhichThreshold = "fdr_lfc", pval.cutoff = 1e-5,
fdr.cutoff = 0.05,lfc.cutoff = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Candidates}{
  A dataframe which contain genomic location and
  read counts of candidate m6A regions.
  All detailed columns can be found from
  the output of function "M6Apeak.MultiRep.step1()".
}
  \item{sf}{
  A vector of size factors for each sample,
  which is the same as that used in the
  first step of identifying candidate regions.
}
  \item{mu.cutoff}{
  A constant specifying the background methylation levels.
  This is estimated automatically based on the first step of peak calling.
}
  \item{WhichThreshold}{
  A character specifying a threshold for significant peaks.
  There are three options: "pval" (only use p-values),
  "fdr" (only use FDR), "lfc" (only use log fold change),
  "pval_lfc" (use both p-values and log fold changes) and
  "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".}
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
This functions maninly call funciton "M6APeak" to rank candidate regions.
}
\value{
The output is a dataframe whose columns are:
\item{chr}{Chromosome number of each peak.}
\item{start}{The start of genomic position of each peak.}
\item{end}{The end of genomic position of each peak.}
\item{strand}{The strand of genomic position of each peak.}
\item{summit}{The summit position of each peak.}
\item{lg.fc}{The log foldchange between normalized IP and
normalized input read counts.}
\item{mu}{The methylation level of each peak if
there are more than one replicate.}
\item{mu.var}{The estimated variance of for
methylation level of each peak, if there are more than one replicate.}
\item{stats}{The Wald test statistics of each peak,
if there are more than one replicate.}
\item{shrkPhi}{The shrinkage estimation of dispersion
for mehtylation levels of each peak, if there are more than one replicate.}
\item{shrkTheta}{The shrinkage estimation for
scale parameter theta in the gamma distribution,
if there are more than one replicate.}
\item{pvals}{ P-value calculated based on the Wald-test.}
\item{p.adj}{Adjusted p-values using FDR.}
\item{rSocre}{A score used to rank each region.
The higher the score, the higher the rank would be.}
Note, there are additional columns whose name involve the character ".bam".
These columns contain the read counts in respective samples.
}

\examples{
### A toy example
data("Basal_regionlevel")
data("Basal_binlevel") ### only use sf
M6Apeak.MultiRep.step2(Candidates = Basal_regionlevel,
sf = Basal_binlevel$sf, mu.cutoff = 0.5)
}
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
