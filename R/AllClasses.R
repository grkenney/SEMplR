
## SNP Effect Matrix (SEM) -----------------------------------------------------
#' https://doi.org/10.1093/bioinformatics/btz612
#' A SNP Effect Matrix (SEM) estimates the binding affinity of every possible
#' mutation in a particular transcription factor (TF) binding motif
#' This class contains three slots: the matrix, the baseline value, and
#' a unique id
#' @slot sem The SEM itself as a data table Rows represent sequence position
#' (variable length), columns represent effects due to each nucleotide base
#' A, C, G, T (fixed length: 4)
#' @slot baseline A scrambled baseline, representing the binding score of
#' randomly scrambled kmers of the same length. This is the binding cutoff
#' for a TF
#' @slot semId basename of the sem file
#'
#' @export
setClass(
  Class = "SNPEffectMatrix",
  slots = c(
    sem = "data.table",
    baseline = "numeric",
    semId = "character"
  )
)


## SEM Collection Class --------------------------------------------------------
#' A collection of SNPEffectMatrix objects and their corresponding meta data
#' @slot sems A list of SNPEffectMatrix objects
#' @slot semData A `data.table` with meta data on each SEM
#' @slot semKey The name of a column in semData that matches the semIds in the 
#' sems list as a `character`. If column entries have a .sem suffix, a new 
#' column named SEM_KEY will be created without the .sem suffixes.
#'
#' @export
setClass(
  Class = "SNPEffectMatrixCollection",
  slots = c(
    sems = "list",
    semData = "data.table",
    semKey = "character"
  )
)


setValidity("SNPEffectMatrixCollection", function(object) {
  # sem_ids must be unique
  sem_ids <- names(object@sems)
  if (length(sem_ids) != length(unique(sem_ids))) {
    "semIds for objects in list must be unique"
  } else {
    TRUE
  }
  
  # semKey column must be unique
  semKeys <- object@semData[, object@semKey]
  if (length(semKeys) != length(unique(semKeys))) {
    paste0("Column designated as semKey must be unique. Number of rows (", 
           length(semKeys), ") != unique keys (", length(unique(semKeys)), ")")
  } else {
    TRUE
  }
  
  # number of sems and semData rows must be same length
  if ((length(object@sems) != nrow(object@semData)) & 
      nrow(object@semData) != 0) {
    paste0("Length of sems (", length(object@sems) ,
           ") and semData (", nrow(object@semData), ") are not equal")
  } else {
    TRUE
  }
  
  # all sem_ids must be in semData keys
  keys <- object@semData[, data.table::key(object@semData), with = FALSE] |>
    unlist()
  id_diffs <- setdiff(sem_ids, keys)
  # if there are sem_ids not in keys and semData is not empty
  if (length(id_diffs) > 0 & length(object@semData) > 0) {
    paste0("not all semIds are in the specified semKey column of semData.\n",
           "\tsemId(s) not in keys: ", paste0(id_diffs, collapse = ", "))
  }
})


## SemplScores class -----------------------------------------------------------

#' Class for storing SEM motif binding calculations for multiple variants
#'
#' @slot variants A `VRanges` object to hold one or more variants
#' @slot scores A `data.table` object for motif information and binding scores
#' @slot semData A `data.table` object of metadata for the SEMs with one 
#' row for each SEM
#'
#' @importFrom VariantAnnotation VRanges
#' @importFrom data.table data.table
#'
#' @export
setClass("SemplScores",
         slots = c(variants = "VRanges",
                   semData = "data.table",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VariantAnnotation::VRanges(),
           semData = data.table::data.table(),
           scores = data.table::data.table()
         )
)


setValidity("SemplScores", function(object) {
  if (nrow(object@scores) > 0) {
    expected_column_names <- c("varId", "semId",
                               "refSeq", "altSeq",
                               "refScore", "altScore",
                               "refNorm", "altNorm")
    actual_column_names <- colnames(object@scores)
    if (sum(expected_column_names %in%
            actual_column_names) != length(expected_column_names)) {
      "@scores must contain columns with names: varId, semId, refSeq, altSeq, refScore, altScore, refNorm, altNorm"
    } else {
      TRUE
    }
  } else {
    TRUE
  }
})


## SequenceFrame -----------------------------------------------------
#' Class to store information about sequencing frame and variant location
#'
#' @slot seqName name to display alongside sequence
#' @slot frameStart index of frame start
#' @slot frameStop index of frame stop
#' @slot sequence character sequence
#' @slot variantIndex index of variant base pair(s)
#' @slot motif sem id to display
#' @slot variantName name of the variant to display
#' 
#' @export
setClass(
  Class = "SequenceFrame",
  slots = c(
    seqName = "character",
    frameStart = "numeric",
    frameStop = "numeric",
    sequence = "character",
    variantIndex = "numeric",
    motif = "character",
    variantName = "character"
  )
)


#' SequenceFrame object and constructor
#'
#' Constructs a SequenceFrame class object.
#'
#' @param seqName name to display alongside sequence
#' @param frameStart index of frame start
#' @param frameStop index of frame stop
#' @param sequence character sequence
#' @param variantIndex index of variant base pair(s)
#' @param motif sem id to display
#' @param variantName name of the variant to display
#' 
#' @importFrom methods new
#'
#' @return a SequenceFrame object
#' @docType class
#' 
#' @export
SequenceFrame <- function(seqName, frameStart, frameStop, 
                          sequence, variantIndex, motif, variantName) {
  new("SequenceFrame",
      seqName = seqName,
      frameStart = frameStart,
      frameStop = frameStop,
      sequence = sequence,
      variantIndex = variantIndex,
      motif = as.character(motif),
      variantName = as.character(variantName))
}
