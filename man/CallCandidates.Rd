\name{CallCandidates}
\alias{CallCandidates}
\title{
Call candidate m6A regions or candidate differential m6A regions.
}
\description{
This function first calls m6A bumps from each pair of input and IP
sample using bin-level data. Then, bumps from all input and IP pairs
are unioned together to obtain a list of candidate regions.
}
\usage{
CallCandidates(Counts, bins,
               WhichThreshold ="fdr_lfc", pval.cutoff = 1e-5,
               fdr.cutoff = 0.05, lfc.cutoff = 0.7,
               windlen = 5,lowcount = 30)
}

\arguments{
\item{Counts}{ A data matrix containing
  bin-level (default 50bp) read counts in both
  IP and input samples,
  where the sample order is: input1, ip1, input2, ip2, ...
  }
  \item{bins}{
  A data frame containing the genomic coordinate of
  each bin of fixed length.
  }
  \item{WhichThreshold}{
  A character specifying a criterion to select significant bins
  in bump finding using an ad hoc algorithm.
  There are five options: "pval" (only use p-values),
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
  A constant indicating the cutoff for
  log fold change. Default is 0.7 for fold change of 2.
  }
  \item{windlen}{
  An integer specifying the length of
  consecutive bins used in simple
  moving average smooth of log fold change. Default is 5.
  }
  \item{lowcount}{
  An integer to filter out candidate regions with lower read counts
  in input. Default is 30.
  }
}

\details{
  The function involves three steps:
\itemize{
 \item Perform binomial test for each bin based bin-level counts
 \item Merge significant bins in each input \& IP pair
  to form bumps usng: \code{\link{findBumps}}
 \item Combine bumps from all input \& IP pairs to construct a list of
 candidate regions.
 }
}

\value{
A list containing
  \item{Regions}{A data frame containng genomic coordinate
  for each candidate region.}
  \item{Counts}{A data matrix containing read counts of
  all samples for each candidate region.}
}



\examples{
### A toy example, whose results do not have real applications.
data("Basal")
Candidates = CallCandidates(
    Counts = Basal$Bins$Counts,
    bins = Basal$Bins$Bins
    )
}
