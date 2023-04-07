\name{findBumps}
\alias{findBumps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Bump-finding from transcriptome bins.
 }
\description{
 This function constructs transcriptome m6A bumps for each input \& IP
 replicate, by merging together bins having significant enrichment of
 IP over input control reads.
 }
\usage{
findBumps(chr, pos, strand, x, count,
          use = "pval",
          pval.cutoff,
          fdr.cutoff,
          lfc.cutoff,
          sep = 2000,
          minlen = 100,
          minCount = 3,
          dis.merge = 100,
          scorefun = mean,
          sort = TRUE)
          }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
  Chromosome number of all bins.
  }
  \item{pos}{
  Transcriptome start position of all bins.
  }
  \item{strand}{ Strand of all bins.
  }
  \item{x}{
  A dataframe containing the p-values, fdrs and log fold changes
  of all bins.
  }
  \item{count}{
  Read counts in each bin from paired input and IP sample.
  }
\item{use}{
  A character to specify which criterion to select significant bins.
  It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": The selection is only based on P-values;
  "fdr": The selection is only based on FDR;
  "lfc": The selection is only based on log fold changes
  between normalized IP and normalized input read counts;
  "pval_lfc": The selection is based on both p-values and log fold changes;
  "fdr_lfc": The selection is based on both FDR and log fold changes.
  Default is "pval".}

  \item{pval.cutoff}{
  A numerical value to specify a cutoff for p-value.
  Default is 1e-5.
  }
  \item{fdr.cutoff}{
  A numerical value to specify a cutoff for fdr. Default is 0.05.
  }
 \item{lfc.cutoff}{
  A numerical value to specify a cutoff for log fold change
  between normalized IP and input read counts.
  Default is 0.7 for fold change of 2.
  }

  \item{sep}{
  A constant used divide genome into consecutive sequenced regions.
  Any two bins with distance greater than sep will be grouped into
  different regions. Default is 2000.
  }
  \item{minlen}{
   A constant to select bumps who have minimum length of minlen.
   Default is 100.
   }
  \item{minCount}{
  A constant to select bumps who have
   at least minlen number of bins. Default is 3.
   }
  \item{dis.merge}{
  A constant. Any twp bumps with distance smaller than dis.merge
  would be merged. Default is 100.
  }
  \item{scorefun}{
  A character indicating a function used to assign a score
  for each bump base on p-values of all spanned bins.
  Default is "mean", meaning that the score is an average of
  bin-level p-values.
  }
  \item{sort}{
  A logical value indicating whether rank (TRUE)
  bumps with the score output from scorefun or not (FALSE).
  Default is TRUE.
  }
}


\value{
  This function returns a dataframe containing the
  chromosome, start position, end position, length, strand,
  summit, total read counts (both IP and input) and score of each
  bump.
  }



\examples{
### Use example dataset "Basal" in TRESS
### to illustrate usage of this function
data("Basal")
bins = Basal$Bins$Bins
Counts = Basal$Bins$Counts
sf = Basal$Bins$sf
colnames(Counts)
dat = Counts[, 1:2]
thissf = sf[1:2]
### pvals based on binomial test
idx = rowSums(dat) > 0
Pvals = rep(1, nrow(dat))
Pvals[idx] = 1 - pbinom(dat[idx, 2],
                        rowSums(dat[idx, ]),
                        prob = 0.5)
### lfc
c0 = mean(as.matrix(dat), na.rm = TRUE)  ### pseudocount
lfc = log((dat[, 2]/thissf[2] + c0)/(dat[, 1]/thissf[1] + c0))
x.vals = data.frame(pvals = Pvals,
                    fdr = p.adjust(Pvals, method = "fdr"),
                    lfc = lfc)
### find bumps based on pvals, fdr or lfc
Bumps = findBumps(chr = bins$chr,
                  pos = bins$start,
                  strand = bins$strand,
                  x = x.vals,
                  use = "fdr_lfc",
                  fdr.cutoff = 0.01,
                  lfc.cutoff = 0.5,
                  count = dat)
head(Bumps, 3)
}
