\name{M6Apeak.MultiRep.step1}
\alias{M6Apeak.MultiRep.step1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The first step of m6A peak calling to identify candidates regions.
}
\description{ This step quickly scan the whole transcriptome,
perform and binomial test and loosely identify candidate regions
using an ad hoc algorithm.
}
\usage{
M6Apeak.MultiRep.step1(Counts, sf, bins, WhichThreshold = "fdr_lfc",
pval.cutoff = 1e-05, fdr.cutoff = 0.05,
lfc.cutoff = 0.7, windlen = 5, lowcount = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{ A data matrix containing bin-level (50bp) read counts
  in both IP and input samples. The number of column depends on the number of
  replicates, where the sample order is: input1, ip1, input2, ip2, ...
}
  \item{sf}{ A vector specifying the size factor of each sample,
  which is estimated using "Counts": colSums(Counts)/median(colSums(Counts)),
  or provided by the user.
}
  \item{bins}{
  A data frame containing the genomic position of each bin of fixed length.
}
  \item{WhichThreshold}{
  A character specifying a threshold for significant bins in bump finding using
  an ad hoc algorithm. There are three options: "pval" (only use p-values),
  "fdr" (only use FDR), "lfc" (only use log fold change),
  "pval_lfc" (use both p-values and log fold changes)
  and "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".
  }
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
  \item{windlen}{
  An integer specifying the length of consecutive bins used in simple
  moving average smooth of log fold change. Default is 5.
}
  \item{lowcount}{
  An iteger to filter out candidate regions with lower read counts
  in input. Default is 30.
}
}
\details{
The function involves three steps:
1. Perform binomial test for each bin based bin-level counts.
   Calculate log fold change for each bin.
2. Use ad hoc algorithm to combine significant bins to form candidate regions.
}
\value{
 \item{chr}{Chromosome number of each m6A candidate region.}
 \item{start}{Start position of each m6A candidate region.}
 \item{end}{End position of each m6A candidate region.}
 \item{summit}{The summit of each m6A candidate region.}
 \item{score}{Adjusted P-value of each candidate region.
              The P-value is derived using on the binomial test
              based on the region-level read
              counts.
              }
 \item{lg.fc}{Log fold changes between IP and input read counts
  within each m6A region.
  }
 \item{Other}{All the other columns corresponds to read
              counts of IP/input within each region from each replicate.
              }
}
\examples{
### A toy example, whose results do not have real applications.
data("Basal_binlevel")
Peak.candidates =
  M6Apeak.MultiRep.step1(
    Counts = Basal_binlevel$bincounts,
    bins = Basal_binlevel$bins,
    sf = Basal_binlevel$sf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
