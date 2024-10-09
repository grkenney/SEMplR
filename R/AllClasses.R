
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
#' @slot sem_id basename of the sem file
#' @slot tf_name Name of the TF relevant to the SEM
#' @slot ensembl_id ENSEMBL gene id of the transcription factor
#' @slot uniprot_id Uniprot protein id of the transcription factor
#' @slot cell_type cell type/line used for ChipSeq experiment
#'
#' @export
setClass(
  Class = "SNPEffectMatrix",
  slots = c(
    matrix = "data.table",
    baseline = "numeric",
    sem_id = "character",
    tf_name = "character",
    ensembl_id = "character",
    uniprot_id = "character",
    cell_type = "character"
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
#' @slot sem_id basename of the sem file
#' @slot tf_name (optional) Name of the TF relevant to the SEM
#' @slot ensembl_id (optional) ENSEMBL gene id of the transcription factor
#' @slot uniprot_id (optional) Uniprot protein id of the transcription factor
#' @slot cell_type (optional) cell type/line used for ChipSeq experiment
#'
#' @importFrom methods new
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @aliases SNPEffectMatrix SNPEffectMatrix-class
#' @rdname SNPEffectMatrix
#' @export
SNPEffectMatrix <- function(matrix, baseline, sem_id, tf_name = "", 
                            ensembl_id = "", uniprot_id = "", cell_type = "") {
  # convert matrix to data.table in case given in another format
  matrix <- data.table::as.data.table(matrix)

  # convert baseline to numeric
  baseline <- as.numeric(baseline)

  new("SNPEffectMatrix",
      matrix = matrix,
      baseline = baseline,
      sem_id = sem_id,
      tf_name = tf_name,
      ensembl_id = ensembl_id,
      uniprot_id = uniprot_id,
      cell_type = cell_type
  )
}


setValidity("SNPEffectMatrix", function(object) {
  if (1 != 1) {
    "length(unique(object$sem_id)) != length(object): sem_id must be unique"
  } else {
    TRUE
  }
})

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
                   sem_metadata = "data.table",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VRanges(),
           sem_metadata = data.table(),
           scores = data.table()
         )
)

#' SemplScores object and constructor
#'
#' Constructs a SemplScores class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param scores A `data.table` object for motif information and binding scores
#' @param sem_metadata A `data.table` object of metadata for the SEMs with one 
#' row for each SEM
#'
#' @importFrom methods new
#'
#' @return a SemplScores object
#' @docType class
#' @aliases SemplScores SemplScores-class
#' @rdname SemplScores
#' @export
SemplScores <- function(variants=NA, scores=NA, sem_metadata=NA) {
  if (is.na(scores)){
    scores_table <- data.table(seqnames=character(),
                               ranges=numeric(),
                               sem=character(),
                               nonRiskSeq=numeric(), riskSeq=numeric(),
                               nonRiskScore=numeric(), riskScore=numeric(),
                               nonRiskNorm=numeric(), riskNorm=numeric())
  }
  else {
    scores_table <- scores
  }
  
  if (is.na(sem_metadata)) {
    sem_metadata <- data.table(tf_name=character(),
                               ensembl_id=character(),
                               uniprot_id=character(),
                               cell_type=character())
  } else {
    sem_metadata <- sem_metadata
  }

  if (all(is.na(variants))) {
    vr = VRanges()
  } else {
    vr = variants
  }

  new("SemplScores",
      variants = vr,
      sem_metadata = sem_metadata,
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


setGeneric("motifScores", function(x, motif) standardGeneric("motifScores"))
setMethod("motifScores", "SemplScores", 
          function(x, motif) x@scores[x@scores$sem %in% motif])

setGeneric("getMotif", function(x, motif) standardGeneric("getMotif"))
setMethod("getMotif", "SemplScores", 
          function(x, motif) x@scores[x@scores$sem == motif])

