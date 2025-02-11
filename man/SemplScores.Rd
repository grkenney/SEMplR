% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{SemplScores}
\alias{SemplScores}
\title{SemplScores object and constructor}
\usage{
SemplScores(variants = NULL, semList = NULL, semScores = NULL)
}
\arguments{
\item{variants}{A \code{VRanges} object to hold one or more variants}

\item{semList}{A named list of SNPEffectMatrix objects}

\item{semScores}{(optional) A \code{data.table} object for motif information and
binding scores}
}
\value{
a SemplScores object
}
\description{
Constructs a SemplScores class object.
}
