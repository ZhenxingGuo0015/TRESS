\name{EstiMu}
\alias{EstiMu}
\title{
Estimation of m6A methylation.
}
\description{
This function calculates, for each candidate region,
the averaged enrichment of normalized IP read counts versus the sum of
normalized IP and input control read counts.
}
\usage{
EstiMu(counts, sf)
}
\arguments{
  \item{counts}{
 A data matrix containing read counts in each region across
sample input1, ip1, input2, ip2, input3, ip3, ...
}
  \item{sf}{
  A numerical vector containing the size factor of each sample,
  which is used for sequencing depth normalization.
  The sample order here is the same as that in counts.
}
}

\value{
\item{mu}{
A numerical vector containing the methylation
level of all candidate regions.
}
}


\examples{
data("Basal")
## methylatinon level
mu = EstiMu(
    counts = Basal$Candidates$Counts,
    sf = Basal$Bins$sf
    )
head(mu, 3)
}
