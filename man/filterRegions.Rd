\name{filterRegions}
\alias{filterRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pre-filtering of candidate DMRs.
}
\description{
This function filters out candidate DMRs who have small
(e.g., less than 25\% quantile)
marginal coefficient of variation (CV) in methylation ratio.
}
\usage{
filterRegions(Candidates, quant = 0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Candidates}{
A list containing "Counts", "Regions", "sf" for read counts,
genomic coordinate of each candidate DMR
and size factor of all samples. See output of \code{\link{CallCandidates}}.
}
\item{quant}{A percentage to specify a quantile cutoff.
Regions with CV smaller than this quantile would be filtered out.
Default is 25\%.}
}

\value{
A list with the same structure with
input "Candidates" but with a smaller number of candidate DMRs.
}

\examples{
# A toy example
data(DMR_M3vsWT) # data from TRESS
sub.DMR_M3vsWT = filterRegions(DMR_M3vsWT)
}

