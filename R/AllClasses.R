
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
#' @importFrom methods new
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @aliases SNPEffectMatrix SNPEffectMatrix-class
#' @rdname SNPEffectMatrix
#' @export
SNPEffectMatrix <- function(matrix, tf_name, baseline, pwm_filename = "") {
  # convert matrix to data.table in case given in another format
  matrix <- data.table::as.data.table(matrix)

  # convert baseline to numeric
  baseline <- as.numeric(baseline)

  new("SNPEffectMatrix",
      matrix = matrix,
      tf_name = tf_name,
      baseline = baseline,
      pwm_filename = pwm_filename
  )
}

## SEMCollection ----------------------------------------------------------------

setClass("SEMCollection",
         representation("VIRTUAL"),
         prototype = prototype(elementType = "SNPEffectMatrix"),
         contains = "list")


## SemplScores class ----------------------------------------------------------------

#' Class for storing SEM motif binding calculations for multiple variants
#'
#' @slot variants A `VRanges` object to hold one or more variants
#' @slot scores A `data.table` object for motif information and binding scores
#'
#' @importFrom VariantAnnotation VRanges
#' @importFrom data.table data.table
#'
#' @export
setClass("SemplScores",
         slots = c(variants = "VRanges",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VRanges(),
           scores = data.table()
         )
)

#' SemplScores object and constructor
#'
#' Constructs a SemplScores class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param scores A `data.table` object for motif information and binding scores
#'
#' @importFrom methods new
#'
#' @return a SemplScores object
#' @docType class
#' @aliases SemplScores SemplScores-class
#' @rdname SemplScores
#' @export
SemplScores <- function(variants=NA, scores=NA) {
  if (is.na(scores)){
    scores_table <- data.table(seqnames=character(),
                               ranges=numeric(),
                               sem=character(),
                               nonRiskSeq=numeric(), riskSeq=numeric(),
                               nonRiskScore=numeric(), riskScore=numeric(),
                               nonRiskNorm=numeric(), riskNorm=numeric())
  }
  else {
    scores_table = scores
  }

  if (all(is.na(variants))) {
    vr = VRanges()
  } else {
    vr = variants
  }

  new("SemplScores",
      variants = vr,
      scores = scores_table
      )
}

setValidity("SemplScores", function(object) {
  expected_column_names <- c("seqnames", "ranges", "sem",
                             "nonRiskSeq", "riskSeq",
                             "nonRiskScore", "riskScore",
                             "nonRiskNorm", "riskNorm")
  actual_column_names <- colnames(object@scores)
  if (sum(expected_column_names %in%
          actual_column_names) != length(expected_column_names)) {
    "@scores must contain columns with names: seqnames, ranges, sem, nonRiskSeq, riskSeq, nonRiskScore, riskScore, nonRiskNorm, riskNorm"
  } else {
    TRUE
  }
})


setGeneric("getMotif", function(x, motif) standardGeneric("getMotif"))
setMethod("getMotif", "SemplScores", function(x, motif) x@scores[x@scores$sem == motif])

getMotif(sem_scores, "BHLHB2_GM12878")
