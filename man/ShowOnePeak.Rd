\name{ShowOnePeak}
\alias{ShowOnePeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visulization of a single peak along the genome.
}
\description{
This function plots the estimated methylation level (as bars) of each bin within a peak for each replicate, and the corresponding normalized input read depth (grey curve).
}
\usage{
ShowOnePeak(onePeak, allBins, binCounts, ext = 500, ylim = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{onePeak}{
A one-row dataframe containing the genomic position of a single peak: chr, start, end, strand.
}
  \item{allBins}{
  A dataframe contaning the genomic position of all bins used to call peaks: chr, start, end, strand.
}
  \item{binCounts}{
  A dataframe contaning the read counts of all bins for each replicate. The sample order is: input1, ip1, input2, ip2, ...
}
  \item{ext}{
  An integer indicating the length of base paris to extend the peak on both sides: (start - ext, end + ext). Default is 500.
}
  \item{ylim}{
  The range of y-axis to plot. Default is c(0, 1)
}
}

\value{
It only generates a plot. No specific output.
}



%% ~Make other sections like Warning with \section{Warning}{....} ~

\seealso{
ShowOneDMR from "DSS" package.
}

\examples{
### read peaks
peaks = read.table(file.path(system.file(package = "TRESS"),
                             "extdata/examplebyBam_peaks.xls"),
                   sep = "\t", header = TRUE)
### load annotation and bin counts
load(file.path(system.file(package = "TRESS"), "extdata/examplebyBam.rda"))
 allBins = as.data.frame(bins$bins)
 colnames(allBins)[1] = "chr"
 allBins$strand = binStrand
 for (i in 1:4) {
   ShowOnePeak(onePeak = peaks[i,], allBins = allBins, binCounts = allCounts)
 }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
