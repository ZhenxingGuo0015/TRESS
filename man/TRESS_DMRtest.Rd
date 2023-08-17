\name{TRESS_DMRtest}
\alias{TRESS_DMRtest}
\title{
Hypothesis testing on candidate DMRs.
}
\description{
This function conducts statistical test for each candidate DMR
based on user specified contrast of coefficients in design.
}
\usage{
TRESS_DMRtest(DMR, contrast, nullModel = "standN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMR}{A list at least containing IP/input ratio,
  the coefficients estimate,
  variance-covariance estimate.
  This can be obtained from the ouput of \code{\link{TRESS_DMRfit}}.}
  \item{contrast}{A contrast for all coefficients in the design.
  It can be either a (p+1) vector or a m-by-(p+1) matrix,
  where p is the number of columns in the design.
  m depends on the number of relationships that users want to test.}
  \item{nullModel}{A character to specify a method to calculate
  p-value based on the statistics. It can be "standN", "2mix" and
  "trunN" for standard normal,
  two-component mixed gaussian and truncated normal respectively.
  Defult is "standN".}
}
\details{
The hypothesis for each of candidate DMR \eqn{i}{i} is of the form:
\deqn{H_0: C^TR_i = 0 vs. H_1: C^TR_i \neq 0}{H_0:
C^TR_i = 0 vs. H_1: C^TR_i \neq 0}
where \eqn{C}{C} is a contrast of all coefficients in model design;
\eqn{R_i}{R_i} is coefficient vector for DMR \eqn{i}{i}.
If the \eqn{C}{C} is a vector, then TRESS performs Wald test;
if the \eqn{C}{C} is a matrix, then TRESS conducts F-test.
}

\value{
This function returns a dataframe containing the testing
results for specified contrast. The columns are
\item{chr}{ Chromosome number of each region}
\item{start}{Start position of each region}
\item{end}{End position of each region}
\item{strand}{Strand information of each region}
\item{gene_Symbol}{If "Path_To_OrgdbSqlite" is specified, this column contains the gene symbol of each region}
\item{baseMean}{Averaged methylation level cross all samples.}
\item{logOR}{Estimated value of contrast: \eqn{C^TR_i}{C^TR_i}. For two-group comparison, it represents the log odds ratio of methylation probabilities.}
\item{lorSE}{Standard error of log odds ratio.}
\item{stat}{Wald statistics.}
\item{pvalue}{P-values from statistical tests.}
\item{padj}{Benjamini-Hochberg procedure adjusted p-values.}
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

\examples{
# A toy example
data(DMR_M3vsWT) # data from TRESS
variable = data.frame(predictor = rep(c("WT", "M3"), c(2, 2)))
model = ~1+predictor
DMR.fit = CallDMRs.paramEsti(
    counts = DMR_M3vsWT$Counts,
    sf = DMR_M3vsWT$sf,
    variable = variable,
    model = model
    )
DMR.fit$Candidates = DMR_M3vsWT
DMR.test = TRESS_DMRtest(DMR = DMR.fit, contrast = c(0, 1))
head(DMR.test, 3)
head(DMR_M3vsWT$Regions[which(DMR.test$padj < 0.05), ], 3)
}

