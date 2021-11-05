\name{CoefName}
\alias{CoefName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the name of each coefficient in design matrix.
}
\description{
This functions returns the name of each coefficient in the design,
for the convenience of constructing correct contrast for
hypothesis testing.
}
\usage{
CoefName(DMR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMR}{
  A list which is the output from \code{\link{TRESS_DMRfit}}.
  }
}

\value{
A character vector containing the name of each coefficient
in design matrix.
}


\examples{
data("DMR_SixWeekvsTwoWeek")
design = data.frame(time = rep(c("2wk", "6wk"), each = 4),
                    region = rep(rep(c("Cortex", "Hypothalamus"),
                                       each = 2),2)
                    )
model = ~1 + time + region + time*region
DMRfit = CallDMRs.paramEsti(
    counts = DMR_SixWeekvsTwoWeek$Counts,
    sf = DMR_SixWeekvsTwoWeek$sf,
    variable = design,
    model = model
    )
CoefName(DMRfit)
}
