\docType{data}
\name{DMR_SixWeekvsTwoWeek}
\alias{DMR_SixWeekvsTwoWeek}
\title{
Transcriptome location and read counts of 200 candidate DMRs.}
\format{
A list with 4 elements
\describe{
  \item{Regions}{A dataframe of 200 obs and 5 variables,
  containing the transcriptome location for 200 candidate
  DMRs from two mouse brain regions at two time points.}
  \item{Counts}{A dataframe of 200 obs and 16 variables,
  caontaining candidate DMR read counts in each
  paired IP and input replicate of all 8 samples.}
  \item{sf}{A numerical vector,
  containing the size factors of all samples
  estimated from the whole transcriptome}
  \item{MeRatio}{
  A data matrix of 200-by-8, containing the
  estimated methylation ratio
  for each sample.}
  ...
}
}
\usage{
data(DMR_SixWeekvsTwoWeek)
}
\description{
A data list which consists of the transcriptome location
and read counts of 200 candidate DMRs from 8 mouse brain samples.
Each sample has two paired IP and input replicates.
It also contains the size factor of each
sample and the
estimated methylation ratio in each sample
for 200 canidate DMRs.
}
\keyword{datasets}
