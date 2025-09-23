## Class Unions ----------------------------------------------------------------
#' Class union for "GRanges-like" objects
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom VariantAnnotation VRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom VariantAnnotation VRanges

#' @import methods
#' @noRd
setClassUnion(
    "GRangesOrVRanges",
    c("GRanges", "VRanges")
)


## SNP Effect Matrix (SEM) -----------------------------------------------------
#' SNP Effect Matrix (SEM)
#'
#' A SNP Effect Matrix (SEM) estimates the binding affinity of every possible
#' mutation in a particular transcription factor (TF) binding motif.
#' Read more about SEMs here: https://doi.org/10.1093/bioinformatics/btz612
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
#' @examples
#' # build a SEM example matrix
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#' 
#' SNPEffectMatrix(sem = m, baseline = 1, semId = "sem_id")
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
#' @examples
#' # build a SEM example matrix
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#' 
#' # build a SNPEffectMatrix object
#' sm <- SNPEffectMatrix(sem = m, baseline = 1, semId = "sem_name")
#' 
#' # create a meta data table
#' md <- data.table::data.table(tf = "tf_name", sem = "sem_name")
#' 
#' # build a collection with 1 SEM
#' SNPEffectMatrixCollection(sems = list(sm), semData = md, semKey = "sem")
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
    if (object@semKey != "") {
        semKeys <- object@semData[, object@semKey]
        if (length(semKeys) != length(unique(semKeys))) {
            paste0(
                "Column designated as semKey must be unique. Number of rows (",
                length(semKeys), ") != unique keys (",
                length(unique(semKeys)), ")"
            )
        } else {
            TRUE
        }
    } else {
        TRUE
    }

    # number of sems and semData rows must be same length
    if ((length(object@sems) != nrow(object@semData)) &
        nrow(object@semData) != 0) {
        paste0(
            "Length of sems (", length(object@sems),
            ") and semData (", nrow(object@semData), ") are not equal"
        )
    } else {
        TRUE
    }

    # all sem_ids must be in semData keys
    keys <- object@semData[, data.table::key(object@semData), with = FALSE] |>
        unlist()
    id_diffs <- setdiff(sem_ids, keys)
    # if there are sem_ids not in keys and semData is not empty
    if (length(id_diffs) > 0 & length(object@semData) > 0) {
        paste0(
            "not all semIds are in the specified semKey column of semData.\n",
            "\tsemId(s) not in keys: ", paste0(id_diffs, collapse = ", ")
        )
    }
})


## SEMplScores class -----------------------------------------------------------

#' Class for storing SEM motif binding calculations for multiple genomic ranges
#' or variants
#'
#' @slot ranges A `GRanges` or `VRanges` object to hold one or more genomic
#' ranges
#' @slot semData A `data.table` object of metadata for the SEMs with one
#' row for each SEM
#' @slot scores A `data.table` object for motif information and binding scores
#'
#' @importFrom VariantAnnotation VRanges
#' @importFrom data.table data.table
#'
setClass("SEMplScores",
    slots = c(
        ranges = "GRangesOrVRanges",
        semData = "data.table",
        scores = "data.table"
    ),
    prototype = list(
        ranges = VariantAnnotation::VRanges(),
        semData = data.table::data.table(),
        scores = data.table::data.table()
    )
)


setValidity("SEMplScores", function(object) {
    # if (nrow(object@scores) > 0) {
    #   expected_column_names <- c("varId", "semId",
    #                              "refSeq", "altSeq",
    #                              "refScore", "altScore",
    #                              "refNorm", "altNorm")
    #   actual_column_names <- colnames(object@scores)
    #   if (sum(expected_column_names %in%
    #           actual_column_names) != length(expected_column_names)) {
    #     "@scores must contain columns with names: varId, semId, refSeq,
    # altSeq, refScore, altScore, refNorm, altNorm"
    #   } else {
    #     TRUE
    #   }
    # } else {
    #   TRUE
    # }
    TRUE
})
