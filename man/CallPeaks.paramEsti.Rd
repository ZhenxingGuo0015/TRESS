\name{CallPeaks.paramEsti}
\alias{CallPeaks.paramEsti}
\title{
Parameter estimation in m6A peak calling with multiple replicates.
}
\description{
This function estimates all involved parameters
in Bayesian hierarchical negative binomial model,
which is built for read counts from candidate regions
generated from multiple input\& IP replicates.
}
\usage{
CallPeaks.paramEsti(mat, sf = NULL, cutoff = NULL,
                    update = "Joint",
                    trans = NULL,
                    optM = "L-BFGS-B",
                    myfscale = -1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
  A matrix containing read counts from all paired
  input \& input replicates. The order of samples are:
  input1, IP1, input2, IP2,...
}
  \item{sf}{
  A vector of size factors for each sample.
  It can be provided by the users or estimated
  automatically from the data. Default is NULL.
}
  \item{cutoff}{
  Background methylation level, which can be automatically
  estimated based on the background read counts in
  IP and input samples, or provided by users. Defauls is NULL.
}
  \item{update}{
  A logical value indicating whether jointly estimating the
  nuisance parameter theta with dispersion parameter phi
  listed in the proposed model. Possible options are "OnlyPhi",
  "Iterative" and "Joint". "OnlyPhi" means only updating phi_i
  using R function \code{\link{optimize}} while fixing
  parameter theta as the plug-in moment estimator;
  "Iterative" means iteratively
  updating and phi using R function \code{\link{optimize}};
  "Joint" means updating them together using R function
  \code{\link{optim}}. Default is "Joint".
}
  \item{optM}{
  A charactor value to specify which optimization algorithm used
  in the R function \code{\link{optim}}.
  The options are: "Nelder-Mead",
  "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent".
  Default is "L-BFGS-B".
  See more details in help pages of \code{\link{optim}}.}
  \item{trans}{
  Needed when **optM == "Nelder-Mead"**. It specifies which
  transformation function used in the estimation of
  dispersion and/or theta parameter(s) which
  are subjected to the nonnegative constraints.
  Possible options are "sin()" and "exp()". Default is NULL.
}
  \item{myfscale}{
  A stop criteria in \code{\link{optim}}. Default is -1e+06.
}

}
\details{
This function mainly involves three estimation procedures:
\itemize{
    \item Estimate methylation levels
    \item Estimate dispersion parameters and
          the variance of the estimated methylation levels
    \item Calculate test statistics and p-values.
    Also, it calculates a score used for peak ranking.
}

}
\value{
 \item{mu}{
 Estimation of methylation levels of all peaks.
 }
 \item{mu.var}{
 Estimated variance for estimated methylation level.
 }
 \item{shrkPhi}{
 Shrinkage estimator for dispersion parameter phi_i.
 }
 \item{shrkTheta}{
 Shrinkage estimator for parameter theta_i
 if update == "Joint" or "Iterative".
 Otherwise it would be a plug-in moment estimator.
 }

 \item{stats}{
 Wald-test statisitcs.
 }
 \item{pvals}{
 P-values derived from normal distribution based on the
 Wald-test statisitcs.
}
 \item{p.adj}{Adjusted p-values using Benjamini-Hochberg procedure.}
 \item{rSocre}{A score used to ranke each region.
 The higher the score, the higher the rank would be.}
}



\examples{
### A toy example using basal samples from mouse cortex
data("Basal")
res = CallPeaks.paramEsti(
    mat = as.matrix(Basal$Candidates$Counts),
    sf = Basal$Bins$sf,
    cutoff = 0.5
    )
}
