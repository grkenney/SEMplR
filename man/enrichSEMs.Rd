% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichSEMs.R
\name{enrichSEMs}
\alias{enrichSEMs}
\title{Calculate enrichment of motifs broken among variant set}
\usage{
enrichSEMs(semScores, d = "changed", lfc = 0.5)
}
\arguments{
\item{semScores}{a SemplScores object}

\item{d}{direction of test. Options are: "changed", "gained", "lost"}

\item{lfc}{Log2FC cutoff (default = 0.5)}
}
\value{
a data.table with results of the enrichment test
\itemize{
\item \code{n.changed}: number of motifs changed in direction specified
\item \code{odds.ratio}: odds ratio of fisher's test
\item \code{ci.lower}: lower 95\% confidence interval
\item \code{ci.upper}: upper 95\% confidence interval
\item \code{pvalue}: p-value of Fisher's test
\item \code{adj.pvalue}: Benjamini & Hochberg adjusted pvalue
}
}
\description{
Performs a Fisher's exact test to test whether a motif is changed in the
specified direction of change more often than expected within the set of
variants.
}
\examples{
library(VariantAnnotation)

# create an SNP Effect Matrix (SEM)
sem <- matrix(rnorm(12), ncol = 4)
colnames(sem) <- c("A", "C", "G", "T")

# create a VRanges object
vr <- VRanges(seqnames = "chr12",
              ranges = 94136009, 
              ref = "G", alt = "C")

# create a list of SNPEffectMatrix objects
semList <- SNPEffectMatrix(sem, baseline = 0.5, semId = "sem_id")

# calculate binding propensity
s <- scoreVariants(vr, semList)

# calculate enrichment
enrichSEMs(s, d = "changed")

}
\keyword{internal}
