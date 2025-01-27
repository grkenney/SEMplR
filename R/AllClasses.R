
## SNP Effect Matrix (SEM) -----------------------------------------------------
#' https://doi.org/10.1093/bioinformatics/btz612
#' A SNP Effect Matrix (SEM) estimates the binding affinity of every possible
#' mutation in a particular transcription factor (TF) binding motif
#' This class contains three slots: the matrix, the baseline value, and
#' metadata, which contains the TF name and
#' @slot sem The SEM itself as a data table Rows represent sequence position
#' (variable length), columns represent effects due to each nucleotide base
#' A, C, G, T (fixed length: 4)
#' @slot baseline A scrambled baseline, representing the binding score of
#' randomly scrambled kmers of the same length. This is the binding cutoff
#' for a TF
#' @slot semId basename of the sem file
#' @slot tf Name of the TF relevant to the SEM
#' @slot ensembl ENSEMBL gene id of the transcription factor
#' @slot uniprot Uniprot protein id of the transcription factor
#' @slot cellType cell type/line used for ChipSeq experiment
#'
#' @export
setClass(
  Class = "SNPEffectMatrix",
  slots = c(
    sem = "data.table",
    baseline = "numeric",
    semId = "character",
    tf = "character",
    ensembl = "character",
    uniprot = "character",
    cellType = "character"
  )
)

#' SNPEffectMatrix object and constructor
#'
#' Constructs a SNPEffectMatrix class object.
#'
#' @param sem A `data.table` object to hold one or more variants
#' @param baseline A `numeric` scrambled baseline, representing the binding
#' score of randomly scrambled kmers of the same length.
#' @param semId basename of the sem file
#' @param tf (optional) Name of the transcription factor relevant to the SEM
#' @param ensembl (optional) ENSEMBL gene id of the transcription factor
#' @param uniprot (optional) Uniprot protein id of the transcription factor
#' @param cellType (optional) cell type/line used for ChipSeq experiment
#'
#' @importFrom methods new
#'
#' @return a SNPEffectMatrix object
#' @docType class
#' @export
SNPEffectMatrix <- function(sem, baseline, semId, tf = "", 
                            ensembl = "", uniprot = "", cellType = "") {
  # convert matrix to data.table in case given in another format
  sem <- data.table::as.data.table(sem)

  # convert baseline to numeric
  baseline <- as.numeric(baseline)

  new("SNPEffectMatrix",
      sem = sem,
      baseline = baseline,
      semId = semId,
      tf = tf,
      ensembl = ensembl,
      uniprot = uniprot,
      cellType = cellType
  )
}


#' @rdname sem
#' @export
setGeneric("sem", function(x) standardGeneric("sem"))

#' Accessor function to pull the matrix slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname sem
#'  
#' @export
setMethod("sem", "SNPEffectMatrix", 
          function(x) x@sem)

#' @rdname baseline
#' @export
setGeneric("baseline", function(x) standardGeneric("baseline"))

#' Accessor function to pull the matrix slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname baseline
#'  
#' @export
setMethod("baseline", "SNPEffectMatrix", 
          function(x) x@baseline)

#' @rdname semId
#' @export
setGeneric("semId", function(x) standardGeneric("semId"))

#' Accessor function to pull the semId slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname semId
#'  
#' @export
setMethod("semId", "SNPEffectMatrix", 
          function(x) x@semId)

#' @rdname tf
#' @export
setGeneric("tf", function(x) standardGeneric("tf"))

#' Accessor function to pull the tf slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname tf
#'  
#' @export
setMethod("tf", signature("SNPEffectMatrix"), 
          function(x) x@tf)

#' @rdname ensembl
#' @export
setGeneric("ensembl", function(x) standardGeneric("ensembl"))

#' Accessor function to pull the ensembl_id slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname ensembl
#'  
#' @export
setMethod("ensembl", "SNPEffectMatrix", 
          function(x) x@ensembl)

#' @rdname uniprot
#' @export
setGeneric("uniprot", function(x) standardGeneric("uniprot"))

#' Accessor function to pull the uniprot slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname uniprot
#'  
#' @export
setMethod("uniprot", "SNPEffectMatrix", 
          function(x) x@uniprot_id)

#' @rdname cellType
#' @export
setGeneric("cellType", function(x) standardGeneric("cellType"))

#' Accessor function to pull the cell_type slot from a SNPEffectMatrix object
#' 
#' @param x a SNPEffectMatrix object
#'  
#' @rdname cellType
#'  
#' @export
setMethod("cellType", "SNPEffectMatrix", 
          function(x) x@cellType)


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
#' @slot metadata A `data.table` object of metadata for the SEMs with one 
#' row for each SEM
#'
#' @importFrom VariantAnnotation VRanges
#' @importFrom data.table data.table
#'
#' @export
setClass("SemplScores",
         slots = c(variants = "VRanges",
                   metadata = "data.table",
                   scores = "data.table"
                   ),
         prototype = list(
           variants = VariantAnnotation::VRanges(),
           metadata = data.table::data.table(),
           scores = data.table::data.table()
         )
)

#' SemplScores object and constructor
#'
#' Constructs a SemplScores class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param semList A named list of SNPEffectMatrix objects
#' @param semScores (optional) A `data.table` object for motif information and 
#' binding scores
#'
#' @importFrom methods new
#' @importFrom VariantAnnotation VRanges
#' @importFrom S4Vectors mcols
#'
#' @return a SemplScores object
#' @docType class
#' @export
SemplScores <- function(variants=NULL, semList=NULL, semScores=NULL) {
  if (all(is.null(variants))) {
    vr <-  VariantAnnotation::VRanges()
  } else {
    vr <- variants
  }
  
  if (length(vr) == 0) {
    S4Vectors::mcols(vr) <- data.frame(id=NA)
  } else if (!("id" %in% names(S4Vectors::mcols(vr)))) {
    S4Vectors::mcols(vr)$id <- 1:length(vr)
  }
  
  if (length(vr) != 0 | length(semList) != 0) {
    vr_sem_combos <- expand.grid(names(semList), S4Vectors::mcols(vr)$id)
    sem_col <- vr_sem_combos[, 1]
    vr_col <- vr_sem_combos[, 2]
  } else {
    sem_col <- character()
    vr_col <- character()
  }

  if (is.null(semScores)){
    empty_scores <- data.table(nonRiskSeq=numeric(), riskSeq=numeric(),
                               nonRiskScore=numeric(), riskScore=numeric(),
                               nonRiskNorm=numeric(), riskNorm=numeric())
    scores_table <- cbind(data.table(varId=vr_col, semId=sem_col), empty_scores)
  } else {
    scores_table <- semScores
  }
  
  if (is.null(semList)) {
    sem_metadata <- data.table(tf=character(),
                               ensembl=character(),
                               uniprot=character(),
                               cellType=character())
  } else {
    sem_metadata <- lapply(semList, 
                           function(x) 
                             data.table::data.table(semId=x@semId,
                                                    tf=x@tf,
                                                    ensembl=x@ensembl,
                                                    uniprot=x@uniprot,
                                                    cellType=x@cellType)) |>
      data.table::rbindlist()
  }

  new("SemplScores",
      variants = vr,
      metadata = sem_metadata,
      scores = scores_table
      )
}


setValidity("SemplScores", function(object) {
  expected_column_names <- c("varId", "semId",
                             "nonRiskSeq", "riskSeq",
                             "nonRiskScore", "riskScore",
                             "nonRiskNorm", "riskNorm")
  actual_column_names <- colnames(object@scores)
  if (sum(expected_column_names %in%
          actual_column_names) != length(expected_column_names)) {
    "@scores must contain columns with names: varId, semId, nonRiskSeq, riskSeq, nonRiskScore, riskScore, nonRiskNorm, riskNorm"
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
                        semList = metadata(x)[semId %in% motif],
                        semScores = x@scores[semId %in% motif])
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
          function(x, motif) x@scores[semId %in% motif])

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
            varId <- NULL
            SemplScores(variants = variants(x)[variants(x)$id %in% v],
                        semList = metadata(x),
                        semScores = x@scores[varId %in% v])
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
#' @export
setMethod("variantScores", "SemplScores", 
          function(x, v) {
            varId <- NULL
            x@scores[varId %in% v]
            })

#' @rdname scores
#' @export
setGeneric("scores", function(x) standardGeneric("scores"))
setGeneric("scores<-", 
           function(x, value) standardGeneric("scores<-"))


#' Accessor scores slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' 
#' @rdname scores
#' @keywords internal
#' @export
setMethod("scores", "SemplScores", 
          function(x) x@scores)

setMethod("scores<-", "SemplScores", function(x, value) {
  x@scores <- value
  x
})

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

#' @rdname metadata
#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))
setGeneric("metadata<-", 
           function(x, value) standardGeneric("metadata<-"))

#' Accessor metadata slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' 
#' @rdname metadata
#' 
#' @export
setMethod("metadata", "SemplScores", 
          function(x) x@metadata)
setMethod("metadata<-", "SemplScores", function(x, value) {
  x@metadata <- value
  x
})

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


## SequenceFrame -----------------------------------------------------
#' Class to store information about sequencing frame and variant location
#'
#' @slot seqName 
#' @slot frameStart
#' @slot frameStop
#' @slot sequence
#'
#' @export
setClass(
  Class = "SequenceFrame",
  slots = c(
    seqName = "character",
    frameStart = "numeric",
    frameStop = "numeric",
    sequence = "character",
    variantIndex = "numeric"
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
#' 
#' @importFrom methods new
#'
#' @return a SequenceFrame object
#' @docType class
#' @export
SequenceFrame <- function(seqName, frameStart, frameStop, 
                          sequence, variantIndex) {
  new("SequenceFrame",
      seqName = seqName,
      frameStart = frameStart,
      frameStop = frameStop,
      sequence = sequence,
      variantIndex = variantIndex)
}


#' Show function for SequenceFrame
#' 
#' @param object a SequenceFrame object
#' 
#' @importFrom crayon bgRed
#' 
#' @rdname show
#' 
#' @export
setMethod("show", "SequenceFrame",
          function(object) {
            for (i in 1:length(object@sequence)) {
              prefix <- substr(object@sequence[i], 
                               start = 1, 
                               stop = object@frameStart[i] - 1)
              
              fp <- substr(object@sequence[i], 
                           object@frameStart[i], 
                           min(object@variantIndex)-1)
              
              v <- substr(object@sequence[i], 
                          min(object@variantIndex),
                          max(object@variantIndex))
              
              fs <- substr(object@sequence[i], 
                           max(object@variantIndex)+1, 
                           object@frameStop[i])
              
              suffix <- substr(object@sequence[i], 
                               start = object@frameStop[i]+1, 
                               stop = nchar(object@sequence[i]))
              
              grey7 <- crayon::make_style("#777777", bg=TRUE)

              cat(object@seqName[i], ": ", prefix, grey7(fp), bgRed(v), grey7(fs), suffix, "\n", sep="")
            }
          }
)
