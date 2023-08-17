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
             Path_To_OrgdbSqlite = NA,
             variable = NULL, model = NULL,
             InputDir, OutputDir = NA,
             experimentName = NA,
             binsize = 50,
             WhichThreshold = "fdr",
             pval.cutoff = 1e-5,
             fdr.cutoff = 0.05,
             lfc.cutoff = 0.4,
             IncludeIntron = TRUE,
             filetype = "bam",
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
\item{Path_To_OrgdbSqlite}{
 A character to specify the path to a "*.sqlite" file containg mappings
 between different version of gene names. If specified, gene SYMBOL for
 each DMR will be provided in the results. Default is NA.
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

\item{WhichThreshold}{
  A character to specify which criterion to select significant bins
  in order to obtain candidate DMRs from the first step.
  It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": The inference is only based on P-values;
  "fdr": The inference is only based on FDR;
  "lfc": The inference is only based on log fold changes
  between normalized IP and normalized input read counts;
  "pval_lfc": The inference is based on both p-values and log fold changes;
  "fdr_lfc": The inference is based on both FDR and log fold changes.
  Default is "fdr".
}
  \item{pval.cutoff}{
  A numerical value to specify a p-value cutoff in the selection of significant bins to form candidate DMRs.
  Default is 1e-5.
  }
  \item{fdr.cutoff}{
  A numerical value to specify a FDR cutoff in the selection of significant bins to form candidate DMRs. Default is 0.05.
  }
 \item{lfc.cutoff}{
  A numerical value to specify a cutoff of log fold change
  between normalized IP and input counts in the selection of significant bins to form candidate DMRs.
  Default is 0.4 for fold change of 1.5.
  }

\item{filetype}{
 A character to specify the format of input data.
 Possible choices are: "bam", "bed" and "GRanges".
 Default is "bam".
 }
\item{IncludeIntron}{
  A logical value indicating whether to include (TRUE) bins overlapping with
  intronic regions or not (False). Default is TRUE.
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
If "OutputDir" is specified, "allBins" and "Candidates" will also be saved under the provided output directory. Detailed structure of  "allBins" and "Candidates" can be found in the manual of \code{\link{DivideBins}} and \code{\link{CallCandidates}}.

In "DMRfit", the elements are
\item{Candidates}{ A list contains:
A data frame containng genomic coordinate
  for each candidate region, and a data matrix containing read counts of
  all samples for each candidate region.
}
\item{Ratio}{A dataframe containing the IP/input ratio from all samples.}
\item{loglik}{A numerical vector containing the log-likelihood of all DMRs.}
\item{Coef}{A matrix containing estimates of coefficients in the design.}
\item{Cov}{A list of variance-covariance matrix estimates for all DMRs.}
}

\author{
Zhenxing Guo <guozhenxing@cuhk.edu.cn>
}

\references{
Zhenxing Guo, Andrew M. Shafik, Peng Jin, Hao Wu. (2022)
Differential RNA Methylation Analysis for MeRIP-seq Data
under General Experimental Design. Bioinformatics, 38 (20), 4705-4712. \url{
https://academic.oup.com/bioinformatics/article/38/20/4705/6692302?login=true}
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
                       WhichThreshold = "fdr",
                       InputDir = InputDir,
                       OutputDir = OutputDir,
                       experimentName = "example"
                       )
                       }
}
