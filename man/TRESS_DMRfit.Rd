\name{TRESS_DMRfit}
\alias{TRESS_DMRfit}
\title{
Differential m6A methylation analysis for MeRIP-seq data
under general experimental design
}
\description{
This function performs differential m6A analysis through the
following three steps:
\itemize{
  \item Divide the whole genome to obtain bin-level read counts:
  \code{\link{DivideBins}}
  \item Call candidate differential m6A methylation regions (DMRs):
  \code{\link{CallCandidates}}
  \item Model fitting on candidate DMRs based on
  Negative Binomial distribution: \code{\link{CallDMRs.paramEsti}}
}

}
\usage{
TRESS_DMRfit(IP.file, Input.file, Path_To_AnnoSqlite,
             variable = NULL, model = NULL,
             InputDir, OutputDir = NA,
             experimentName = NA,
             binsize = 50,
             filetype = "bam",
             IncludeIntron = FALSE,
             filterRegion = TRUE,
             shrkPhi = TRUE,
             addsuedo = FALSE)
}

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
  A character to specify the path to a "*.sqlite" file used
  for genome annotation.
  }
\item{variable}{
  A dataframe containing condition information of all samples.
  Default is NULL.
  }
\item{model}{
  A formula to specify which factor in "variable"
  will be included into design for model fitting. Default is NULL.
  }
\item{InputDir}{
 A character to specify the input directory of all BA, files.
}
\item{OutputDir}{
 A character to specify an output directory to save
 bin-level and region-level data.
 Default is NA, which will not save any results.
 }

\item{experimentName}{
 A character to specify the name of results if "OutputDir" is provided.
 }

\item{binsize}{
  A numerical value to specify the size of window to
  bin the genome.
  Default value is 50.}

\item{filetype}{
 A character to specify the format of input data.
 Possible choices are: "bam", "bed" and "GRanges".
 Default is "bam".
 }
\item{IncludeIntron}{
  A logical value indicating whether to include (TRUE)
  intronic regions or not (False). Default is FALSE.
  }

\item{filterRegion}{
  A logical value indicating whether to filter out candidate DMRs
  based on their marginal coefficient of variation (CV) in
  methylation ratios.
  If TRUE, then a candidate DMR with CV < 25\% quantile
  would be filtered out. Default value is TRUE.
  }

\item{shrkPhi}{
  A logical value to indicate whether conducting
  shringkage estimate for dispersion parameter. Default is TRUE.}

\item{addsuedo}{
  A logical value to indicate whether or not adding a
  psuedo count 5 on raw read counts. Default is FALSE.
  }
}

\details{
For complete details on each step (especially step 3)
in above "Description" section,
please see the manual pages of respective functions.
}
\value{
This function generates three sets of results: "allBins",
"Candidates" and "DMRfit" returned respectively by function \code{\link{DivideBins}}, \code{\link{CallCandidates}}
and \code{\link{CallDMRs.paramEsti}}.
If "OutputDir" is not specified, only "DMRfit" will be returned.
If "OutputDir" is specified, in addition to returning "DMRfit",
"allBins" and "Candidates" will also be saved
under the provided output directory.

Detailed structure of "DMRfit", "Candidates" and "allBins"
can be found in the manual of respective functions.
}

\author{
Zhenxing Guo <zhenxing.guo@emory.edu>
}

\references{
Zhenxing Guo, Andrew M. Shafik, Peng Jin, Hao Wu.
"Differential RNA Methylation Analysis for MeRIP-seq Data
under General Experimental Design."
}

\seealso{
\code{\link{DivideBins}}, \code{\link{CallCandidates}},
\code{\link{CallDMRs.paramEsti}}
}
\examples{
\dontrun{
Input.file = c("input1.bam", "input2.bam",..., "inputN.bam")
IP.file = c("ip1.bam", "ip2.bam", ..., "ipN.bam")
InputDir = "/directory/to/BAMfile"
OutputDir = "/directory/to/output"
Path_sqlit = "/path/to/xxx.sqlite"
design = "YourDesign"
model = "YourModel"
DMR.fit = TRESS_DMRfit(IP.file = IP.file,
                       Input.file = Input.file,
                       Path_To_AnnoSqlite = Path_sqlit,
                       variable = design,
                       model = model,
                       InputDir = InputDir,
                       OutputDir = OutputDir,
                       experimentName = "example"
                       )
                       }
}
