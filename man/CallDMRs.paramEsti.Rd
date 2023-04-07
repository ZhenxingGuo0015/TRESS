\name{CallDMRs.paramEsti}
\alias{CallDMRs.paramEsti}
\title{
Model fitting and parameter estimation by TRESS for each candidate DMR.
}
\description{
TRESS models the read counts in candidate DMR using hierarchical negative
binomial distribution, with methylation level of each
DMR linked to multi-factors in the design by a linear framework.
This function conducts model fitting, parameter estimation,
and the variance-covariance matrix computation.
}
\usage{
CallDMRs.paramEsti(counts, sf,
                   model, variable,
                   shrkPhi = TRUE,
                   addsuedo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
  A dataframe containing read counts in each candidate DMR
  across all samples.
}
\item{sf}{A numerical vector of
size factors for all samples.}

\item{variable}{A dataframe containing condition information
of all samples.}

\item{model}{ A formula to specify which factor in "variable"
to be included in model fitting.}

\item{shrkPhi}{A logical value indicating whether conducting
shringkage estimate for dispersion parameter. Default is TRUE.}

\item{addsuedo}{A logical value indicating whether or not adding a
psuedo count of 5 on raw read counts. Default is FALSE. }
}

\value{
This function returns a list containing:
\item{Ratio}{A dataframe containing the IP/input ratio from all samples.}
\item{loglik}{A numerical vector containing the log-likelihood of all DMRs.}
\item{Coef}{A matrix containing estimates of coefficients in the design.}
\item{Cov}{A list of variance-covariance matrix estimates for all DMRs.}
%\item{shrkPhi}{A vector containing shrinkage estimate of
%methylation dispersion.}
%\item{shrkTheta}{A vector containing shrinkage estimate for
%the scale parameter in gamma model.}
%\item{id.rm}{A vector containing the index of
%candidate DMRs whose methylation level are smaller than 0.5
%across all samples.}
}


\examples{
# A toy example
data(DMR_M3vsWT) # data from TRESS
variable = data.frame(predictor = rep(c("WT", "M3"), c(2, 2)))
model = ~1+predictor
DMRfit = CallDMRs.paramEsti(
    counts = DMR_M3vsWT$Counts,
    sf = DMR_M3vsWT$sf,
    variable = variable,
    model = model
    )
}
