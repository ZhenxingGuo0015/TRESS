\name{meRatio}
\alias{meRatio}
\title{
Observed m6A methylation ratio.
}
\description{
This function calculates, for each candidate region,
the enrichment of normalized IP read counts versus the sum of
normalized IP and input control read counts.
}
\usage{
meRatio(counts, sf)
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
\item{ratio}{
A numerical data matrix containing the methylation
ratio of each candidate region across all samples.
Here, the number of columns is half of the number of columns in
read count matrix.
}
}


\examples{
data("Basal")
## methylatinon ratio
Ratio = meRatio(
    counts = Basal$Candidates$Counts,
    sf = Basal$Bins$sf
    )
head(Ratio, 3)
}
