
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
#' @aliases SNPEffectMatrix-class
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
#' @param sem_id basename of the sem file
#' @param tf_name (optional) Name of the TF relevant to the SEM
#' @param ensembl_id (optional) ENSEMBL gene id of the transcription factor
#' @param uniprot_id (optional) Uniprot protein id of the transcription factor
#' @param cell_type (optional) cell type/line used for ChipSeq experiment
#'
#' @importFrom methods new
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @aliases SNPEffectMatrix-class
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
#' @slot sem_metadata A `data.table` object of metadata for the SEMs with one 
#' row for each SEM
#'
#' @importFrom VariantAnnotation VRanges
#' @importFrom data.table data.table
#' 
#' @aliases SemplScores-class
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
#' @aliases SemplScores-class
#' @rdname SemplScores
#' @export
SemplScores <- function(variants=NULL, scores=NULL, sem_metadata=NULL) {
  if (is.null(scores)){
    scores_table <- data.table(varId=character(),
                               sem=character(),
                               nonRiskSeq=numeric(), riskSeq=numeric(),
                               nonRiskScore=numeric(), riskScore=numeric(),
                               nonRiskNorm=numeric(), riskNorm=numeric())
  }
  else {
    scores_table <- scores
  }
  
  if (is.null(sem_metadata)) {
    sem_metadata <- data.table(tf_name=character(),
                               ensembl_id=character(),
                               uniprot_id=character(),
                               cell_type=character())
  } else {
    sem_metadata <- sem_metadata
  }

  if (all(is.null(variants))) {
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
  expected_column_names <- c("varId", "sem",
                             "nonRiskSeq", "riskSeq",
                             "nonRiskScore", "riskScore",
                             "nonRiskNorm", "riskNorm")
  actual_column_names <- colnames(object@scores)
  if (sum(expected_column_names %in%
          actual_column_names) != length(expected_column_names)) {
    "@scores must contain columns with names: varId, sem, nonRiskSeq, riskSeq, nonRiskScore, riskScore, nonRiskNorm, riskNorm"
  } else {
    TRUE
  }
})

#' @rdname motifSub
#' @export
setGeneric("motifSub", function(x, motif) standardGeneric("motifSub"))

#' Accessor function to subset a SemplScores object by a motif
#' 
#' @param x a SemplScores object
#' @param motif a sem id represented in x
#' 
#' @rdname motifSub
#' 
#' @export
setMethod("motifSub", "SemplScores", 
          function(x, motif) {
            SemplScores(variants = variants(x),
                        sem_metadata = sem_metadata(x)[sem_id %in% motif],
                        scores = x@scores[sem %in% motif])
                        })

#' @rdname motifScores
#' @export
setGeneric("motifScores", function(x, motif) standardGeneric("motifScores"))

#' Accessor for sem binding scores for a motif
#' 
#' @param x a SemplScores object
#' @param motif a sem id represented in x
#' 
#' @rdname motifScores
#' 
#' @export
setMethod("motifScores", "SemplScores", 
          function(x, motif) x@scores[sem %in% motif])

#' @rdname variantSub
#' @export
setGeneric("variantSub", function(x, v) standardGeneric("variantSub"))

#' Accessor function to subset a SemplScores object by a variant
#' 
#' @param x a SemplScores object
#' @param v a variant id represented in x
#' 
#' @rdname variantSub
#' 
#' @export
setMethod("variantSub", "SemplScores", 
          function(x, v) {
            SemplScores(variants = variants(x)[variants(x)$id %in% v],
                        sem_metadata = sem_metadata(x),
                        scores = x@scores[varId %in% v])
          })

#' @rdname variantScores
#' @export
setGeneric("variantScores", function(x, v) standardGeneric("variantScores"))

#' Accessor for sem binding scores for a variant
#' 
#' @param x a SemplScores object
#' @param v a variant id represented in x
#' 
#' @rdname variantScores
#' 
#' @export
setMethod("variantScores", "SemplScores", 
          function(x, v) x@scores[varId %in% v])

#' @rdname scores
#' @export
setGeneric("scores", function(x) standardGeneric("scores"))

#' Accessor scores slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' 
#' @rdname scores
#' 
#' @export
setMethod("scores", "SemplScores", 
          function(x) x@scores)

#' @rdname variants
#' @export
setGeneric("variants", function(x) standardGeneric("variants"))

#' Accessor variants slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' 
#' @rdname variants
#' 
#' @export
setMethod("variants", "SemplScores", 
          function(x) x@variants)

#' @rdname sem_metadata
#' @export
setGeneric("sem_metadata", function(x) standardGeneric("sem_metadata"))

#' Accessor sem_metadata slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' 
#' @rdname sem_metadata
#' 
#' @export
setMethod("sem_metadata", "SemplScores", 
          function(x) x@sem_metadata)

#' @rdname tf_name
#' @export
setGeneric("tf_name", function(x, s) standardGeneric("tf_name"))

#' Accessor function to pull tf_name for a given sem id
#' 
#' @param x a SemplScores object
#' @param s a sem id represented in x
#'  
#' @rdname tf_name
#'  
#' @export
setMethod("tf_name", "SemplScores", 
          function(x, s) x@sem_metadata[sem_id %in% s]$tf_name)

#' @rdname ensembl_id
#' @export
setGeneric("ensembl_id", function(x, s) standardGeneric("ensembl_id"))

#' Accessor function to pull ensembl_id for a given sem id
#' 
#' @param x a SemplScores object
#' @param s a sem id represented in x
#'  
#' @rdname ensembl_id
#'  
#' @export
setMethod("ensembl_id", "SemplScores", 
          function(x, s) x@sem_metadata[sem_id==s]$ensembl_id)

#' @rdname uniprot_id
#' @export
setGeneric("uniprot_id", function(x, s) standardGeneric("uniprot_id"))

#' Accessor function to pull uniprot_id for a given sem id
#' 
#' @param x a SemplScores object
#' @param s a sem id represented in x
#' 
#' @rdname uniprot_id
#' 
#' @export
setMethod("uniprot_id", "SemplScores", 
          function(x, s) x@sem_metadata[sem_id==s]$uniprot_id)

#' @rdname cell_type
#' @export
setGeneric("cell_type", function(x, s) standardGeneric("cell_type"))

#' Accessor function to pull cell_type for a given sem id
#' 
#' @param x a SemplScores object
#' @param s a sem id represented in x
#' 
#' @rdname cell_type
#' 
#' @export
setMethod("cell_type", "SemplScores", 
          function(x, s) x@sem_metadata[sem_id==s]$cell_type)

#' @rdname changed_motif
#' @export
setGeneric("changed_motif", 
           function(scores_table, direction="changed") 
             standardGeneric("changed_motif"))

#' Accessor function to subset the scores slot to changed motifs
#' 
#' @param scores_table the scores slot of a SemplScores object
#' @param direction direction of binding change. options are: 
#' 'changed', 'gained', 'lost', 'maintained'
#' 
#' @rdname changed_motif
#' 
#' @export
setMethod("changed_motif", "data.table", 
          function(scores_table, direction="changed") {
            if (direction == "gained"){
              scores_table[(scores_table$nonRiskNorm < 0) & (scores_table$riskNorm > 0), ]
            } else if (direction == "lost") {
              scores_table[(scores_table$nonRiskNorm > 0) & (scores_table$riskNorm < 0), ]
            } else if (direction == "maintained") {
              scores_table[(scores_table$nonRiskNorm > 0) & (scores_table$riskNorm > 0), ]
            } else if (direction == "changed") {
              scores_table[((scores_table$nonRiskNorm < 0) & (scores_table$riskNorm > 0)) |
                         ((scores_table$nonRiskNorm > 0) & (scores_table$riskNorm < 0)), ]
            } else {
              stop("direction is not valid. Options are 'gained', 'lost', 'maintained', or 'changed'")
            }
          }
          )

