\name{M6Apeak}
\alias{M6Apeak}
\title{
Parameter estimation in m6A peak calling problem.
}
\description{
This is a wrap function to estimate all involved parameters in the proposed Bayesian hierarchical negative binomial model, which is the second step of our proposed peak calling method.
}
\usage{
M6Apeak(mat, sf = NULL, cutoff = NULL, update = "Joint",
trans = NULL, optM = "L-BFGS-B", myfscale = -1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A count matrix containing both IP and input samples. The order of columns are: input1, IP1, input2, IP2,...
}
  \item{sf}{
  A numerical vector to record the size factor for each sample. It can be mannually provided by the users or estimated automatically from the count data itself. Default is NULL.
}
  \item{cutoff}{
  The background methylation level, which is automatically estimated based on the background read counts in IP and input samples.
}
  \item{update}{
  A logical value indicating whether jointly estimating the nuisance parameter theta with dispersion parameter phi listed in the proposed model. Possible options are "OnlyPhi", "Iterative" and "Joint". "OnlyPhi" means only updating phi_i using R function "optimize" while fixing parameter theta as the plug-in moment estimator.; "Iterative" means iteratively updating theta and phi using R function "optimize"; "Joint" means updating them together using R function "optim". Default is "Joint".
}
  \item{optM}{
A charactor value to specify which optimization algorithm used in the R function ``optim". The options are: "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". Default is "L-BFGS-B". See more details in help of R function "optim".}
  \item{trans}{
  Needed when "optM == Nelder-Mead". It specifies a transformation function used in the estimation of dispersion and/or theta parameter(s) which are subjected to the nonnegative constraints. Possible options are "sin" and "exp".
}
  \item{myfscale}{
The stop criteria in "optim". Default is -1e+06. See parameter "fnscale" used in "control" option of function "optim".
}

}
\details{
This function mainly runs three set of estimation procedures:
1. Estimate methylation levels with the sub-function "EstiMu";
2. Estimate dispersion parameters with the sub-function "EstiPhi";
3. Estimate the variance of the estimated methylation levels with the sub-function "EstiVarofMu";
4. Calculate test statistics and p-values with sub-function "getInfer". Also, it calculates a new score used for peak ranking.
}
\value{
\item{mu}{
Estimation of methylation levels of all sites.
}
\item{mu.var}{
Estimated variance of the estimates of methylation level.
}
\item{shrkPhi}{
Shrinkage estimator for dispersion parameter phi_i.
}
\item{shrkTheta}{
Shrinkage estimator for parameter theta_i if update == "Joint" or "Iterative". Otherwise it would be a plug-in moment estimator.
}

\item{stats}{
Wald-test statisitcs.
}
\item{pvals}{
P-values derived from normal distribution based on the Wald-test statisitcs.
}
\item{p.adj}{Adjusted p-values using FDR}
\item{rSocre}{A score used to ranke each region. The higher the score, the higher the rank would be.}

}

\examples{
 ### A toy example using basal samples from mouse cortex
 data("Basal_regionlevel")
 ix = grepl("rep", colnames(Basal_regionlevel))
 res = M6Apeak(mat = as.matrix(Basal_regionlevel[, ix]), cutoff = 0.5)
}
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
