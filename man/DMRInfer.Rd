\name{DMRInfer}
\alias{DMRInfer}
\title{
P-value calculation given Wald statistics.
}
\description{
This function calculates p-values for candidate DMRs given
Wald statistics, using respectively
two-component mixtures of normal, truncated normal and standard normal distributions.
}
\usage{
DMRInfer(stat, nullModel = "standN")
}
\arguments{
  \item{stat}{
   A vector containing Wald statistics of all candidate DMRs.
}
  \item{nullModel}{A character to specify a method to calculate
  p-value based on the statistics. It can be "standN", "2mix" and
  "trunN" for standard normal,
  two-component mixed gaussian and truncated normal respectively.
  Defult is "standN".}
}
\details{
In addition to standard normal distribution, TRESS provides
another two distributions to calculate p-values given
Wald statistics, in case that statistics are inflated
by potentially underestimated dispersion in data.

One is two-component mixtures of normal distribution,
where TRESS assumes that,
\deqn{T_w \sim p N(0, \sigma_0^2) + (1-p)N(0, \sigma_1^2)}{
T_w \sim p N(0, \sigma_0^2) + (1-p)N(0, \sigma_1^2)}
where \eqn{T_w}{T_w} is Wald statistics, \eqn{\sigma_0}{\sigma_0} and
\eqn{\sigma_1}{\sigma_1} are standard deviations of distribution that
\eqn{T_w}{T_w} follows under null and alternative hypothesis in
\code{\link{TRESS_DMRtest}}. P-values for \eqn{T_w}{T_w}
are calculated using
\eqn{N(0, \hat{\sigma}_0^2)}{N(0, \hat{\sigma}_0^2)}.


The other one is truncated normal distribution,
where TRESS assumes that,
\deqn{ tT_w\sim p N(0, \sigma^2)}{tT_w\sim p N(0, \sigma^2)}
where \eqn{tT_w = T_w \in [-b, b]}{tT_w = T_w \in [-b, b]}.
Here, \eqn{T_w}{T_w} within range [-b, b] is assumed to
sampled from null
distribuion \eqn{N(0, \sigma^2)}{N(0, \sigma^2)}.
For this truncated normal distribution, TRESS explores
different values for boundary \eqn{b}{b} ranging from 1.5 to
2 by step 0.1.
TRESS estimates a \eqn{\hat{\sigma}}{\hat{\sigma}}
for each of 6 boundaries.
If \eqn{\hat{\sigma}_{max} - \hat{\sigma}_{min} > 0.5
}{\hat{\sigma}_{max} - \hat{\sigma}_{min} > 0.5},
TRESS calculates p-values for \eqn{T_w}{T_w} using
\eqn{N(0, \hat{\sigma}_{min}^2)}{N(0, \hat{\sigma}_{min}^2)}.
Otherwise, p-values are obtained using
\eqn{N(0, \hat{\sigma}^2)}{N(0, \hat{\sigma}^2)},
with \eqn{\hat{\sigma}}{\hat{\sigma}} estimated under
\eqn{b=2}{b=2}.

}


\value{
This function returns a dataframe containing p-values and
Benjamini-Hochberg procedure adjusted p-values.
}

\examples{
### use a randomly generated toy data as an
### illustrate of DMRInfer
set.seed(12345)
p = 0.8
nsites = 10000
flag.TP = rep(NA, nsites)
T_w = rep(NA, nsites)
for (i in seq_len(nsites)) {
  u = runif(1, min = 0, max = 1)
  if(u < p){
    flag.TP[i] = FALSE
    T_w[i] = rnorm(1, 0, sd = 1)
  }else{
    flag.TP[i] = TRUE
    T_w[i] = rnorm(1, 0, sd = 5)
  }
}
res = DMRInfer(stat = T_w, nullModel = "standN")
sum(res$padj < 0.05 & !flag.TP)/sum(res$padj < 0.05)

res = DMRInfer(stat = T_w, nullModel = "2mix")
sum(res$padj < 0.05 & !flag.TP)/sum(res$padj < 0.05)

res = DMRInfer(stat = T_w, nullModel = "trunN")
sum(res$padj < 0.05 & !flag.TP)/sum(res$padj < 0.05)
}

