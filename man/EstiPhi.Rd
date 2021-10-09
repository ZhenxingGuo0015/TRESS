\name{EstiPhi}
\alias{EstiPhi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dispersion estimation.
}
\description{
This is a wrapper function to estimate the posterior of
methylation dispersion for each candidate region.
}
\usage{
EstiPhi(counts, sf,
        update = c("OnlyPhi", "Iterative", "Joint"),
        optM = "L-BFGS-B", myfscale = -1e+6, trans = "sin")
}
%- maybe also 'usage' for other objects documented here.
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
  \item{update}{
  A character specifying which strategy used for estimation of
  dispersion. There are three options: "OnlyPhi", "Iterative" and
  "Joint". "OnlyPhi" means only estimate posterior of phi with
  theta fixed as its moment estimate. "Iterative" means iteratively
  update phi and theta based on their posterior distribution. "Joint"
   means jointly estimate phi and theta based on their joint posterior.
}

  \item{optM}{
  A character to specify which maximizing algorithm used for
  optimization.
  Default is "L-BFGS-B". See \code{\link{optim}} for more details.
}
  \item{myfscale}{
An overall scaling to be applied to the value of fn and gr during optimization. If negative, turns the problem into a maximization problem. Optimization is performed on fn(par)/fnscale. Default is -1e+6. See
\code{\link{optim}} for more details.
}
  \item{trans}{
  A character specifying which transformation function used for phi,
  in the process of
  obtaining its posterior mode.
  There are two options: "sin" and "exp" for transformation
  \eqn{phi = sin(s)+1)/2}{phi = sin(s)+1)/2} and
  \eqn{phi = exp(-exp(s))}{phi = exp(-exp(s))}, where maximization is performed
  on s. Default is "sin".
  Note that, this argument only works when the value of
  optM is not "L-BFGS-B".
}
}

\value{
This function returns a dataframe containing: phi and theta
estimates for all candidate regions.
}


\examples{
data("Basal")
## methylatinon level
res = EstiPhi(counts = as.matrix(Basal$Candidates$Counts),
              sf = Basal$Bins$sf,
              update = "Joint")
head(res, 5)
}
