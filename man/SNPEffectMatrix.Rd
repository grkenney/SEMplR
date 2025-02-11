% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{SNPEffectMatrix}
\alias{SNPEffectMatrix}
\title{SNPEffectMatrix object and constructor}
\usage{
SNPEffectMatrix(
  sem,
  baseline,
  semId,
  tf = "",
  ensembl = "",
  uniprot = "",
  cellType = ""
)
}
\arguments{
\item{sem}{A \code{data.table} object to hold one or more variants}

\item{baseline}{A \code{numeric} scrambled baseline, representing the binding
score of randomly scrambled kmers of the same length.}

\item{semId}{basename of the sem file}

\item{tf}{(optional) Name of the transcription factor relevant to the SEM}

\item{ensembl}{(optional) ENSEMBL gene id of the transcription factor}

\item{uniprot}{(optional) Uniprot protein id of the transcription factor}

\item{cellType}{(optional) cell type/line used for ChipSeq experiment}
}
\value{
a SNPEffectMatrix object
}
\description{
Constructs a SNPEffectMatrix class object.
}
