\name{DivideBins}
\alias{DivideBins}
\title{
Obtain genomic bins and bin-level read counts from BAM files.
}
\description{
This function first divides the whole genome into equal-sized bins
and then calculates read counts in each bin for all samples.
The number of bins depends on the input annotation file, bin size and
whether or not including intronic regions.
}
\usage{
DivideBins(IP.file, Input.file, Path_To_AnnoSqlite,
           InputDir,OutputDir, experimentName,
           binsize = 50, filetype = "bam",
           IncludeIntron = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP.file}{
  A vector of characters containing the name of BAM files
  for all IP samples.
  }
  \item{Input.file}{
  A vector of characters containing the name of BAM files
  for all input control samples.
  }

  \item{Path_To_AnnoSqlite}{
  A character to specify the path to a "*.Sqlite" file used
  for genome annotation.
  }

  \item{experimentName}{
  A character to specify a name for output results.
  }

  \item{binsize}{
  A numerical value to specify a size of window to
  bin the genome and get bin-level read counts.
  Default value is 50.
  }


  \item{filetype}{
  A character to specify the format of input data.
  Possible choices are: "bam", "bed" and "GRanges".
  Default is "bam".
  }

  \item{InputDir}{
  A character to specify the input directory of all BAM files.
  }

  \item{OutputDir}{
  A character to specify an output directory to save all results.
  Default is NA, which will not save any results.
  }

  \item{IncludeIntron}{
  A logical value indicating whether to include (TRUE)
  intronic regions or not (False). Default is FALSE.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The value returned by this function is a list containing two components:
  \item{bins}{A dataframe containing the genomic cordinates of all bins.}
  \item{binCount}{A M-by-N matrix containg bin-level read counts in M bins
  and N samples, where N is two times of the
  length of "IP.file" or "Input.file".
  The column order depends on the sample order in
  "Input.file" and "IP.file".
  }
If the "OutputDir" is specified, then both
genomic bins and corresponding bin-level
read counts would be saved as an ".rda" file.
}



\examples{
# use data in pakage datasetTRES
# available on github, which can be installed by
# install_github("https://github.com/ZhenxingGuo0015/datasetTRES")
\dontrun{
library(datasetTRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
Path_sqlit = file.path(system.file(package = "datasetTRES"),
 "extdata/mm9_chr19_knownGene.sqlite")
#OutDir = "/Users/zhenxingguo/Documents/research/m6a/packagetest"
allBins = DivideBins(
    IP.file = IP.file,
    Input.file = Input.file,
    Path_To_AnnoSqlite = Path_sqlit,
    InputDir = BamDir
    )
 }

}

