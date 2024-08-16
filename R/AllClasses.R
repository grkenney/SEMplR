
## SNP Effect Matrix (SEM) -----------------------------------------------------
#' https://doi.org/10.1093/bioinformatics/btz612
#' A SNP Effect Matrix (SEM) estimates the binding affinity of every possible
#' mutation in a particular transcription factor (TF) binding motif
#' This class contains three slots: the matrix, the baseline value, and
#' metadata, which contains the TF name and
#' @slot matrix The SEM itself as a data table Rows represent sequence position
#' (variable length), columns represent effects due to each nucleotide base
#' A, C, G, T (fixed length: 4)
#' @slot baseline A scrambled baseline, representing the binding score of
#' randomly scrambled kmers of the same length. This is the binding cutoff
#' for a TF
#' @slot tf_name Name of the TF relevant to the SEM
#' @slot pwm_filename Name of the original PWM file used to generate the SEM
#'
#' @export
setClass(
  Class = "SNPEffectMatrix",
  slots = c(
    matrix = "data.table",
    tf_name = "character",
    baseline = "numeric",
    pwm_filename = "character"
  )
)

#' SNPEffectMatrix object and constructor
#'
#' Constructs a SNPEffectMatrix class object.
#'
#' @param matrix A `data.table` object to hold one or more variants
#' @param tf_name A `character` name for the transcription factor
#' @param baseline A `numeric` scrambled baseline, representing the binding
#' score of randomly scrambled kmers of the same length.
#' @param pwm_filename A `character` name of the original PWM file used to
#' generate the SEM
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @aliases SNPEffectMatrix-class
#' @rdname SNPEffectMatrix
#' @export
SNPEffectMatrix <- function(matrix, tf_name, baseline, pwm_filename = "") {
  # convert matrix to data.table in case given in another format
  matrix <- as.data.table(matrix)

  # convert baseline to numeric
  baseline <- as.numeric(baseline)

  new("SNPEffectMatrix",
      matrix = matrix,
      tf_name = tf_name,
      baseline = baseline,
      pwm_filename = pwm_filename
  )
}

## SemplR class ----------------------------------------------------------------

#' Class for storing SEM motif binding calculations for multiple variants
#'
#' @slot variants A `VRanges` object to hold one or more variants
#' @slot scores A `data.table` object for motif information and binding scores
#'
#' @export
setClass("SemplR",
         slots = c(variants = "VRanges",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VRanges(),
           scores = data.table()
         )
)

#' SemplR object and constructor
#'
#' Constructs a SemplR class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param scores A `data.table` object for motif information and binding scores
#'
#' @return a SemplR object
#' @docType class
#' @aliases SemplR-class
#' @rdname SemplR
#' @export
SemplR <- function(variants, scores) {
  new("SemplR",
      variants = VRanges(),
      scores = data.table()
      )
}
