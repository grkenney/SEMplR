# ---- helper functions ----

# Given a single variant in a VRanges object, construct a unique id from the 
# position information
.makeDefaultId <- \(vr) {
  start_pos <- IRanges::start(IRanges::ranges(vr))
  end_pos <- IRanges::end(IRanges::ranges(vr))
  sn <- GenomeInfoDb::seqnames(vr)
  
  ref_allele <- as.character(VariantAnnotation::ref(vr))
  alt_allele <- as.character(VariantAnnotation::alt(vr))
  
  if (ref_allele == "") {
    allele_str <- paste0("ins", alt_allele)
  } else if (alt_allele == "") {
    allele_str <- paste0("del", ref_allele)
  } else {
    allele_str <- paste0(ref_allele, ">", alt_allele)
  }
  
  pos_str <- ifelse(start_pos == end_pos, 
                    start_pos, paste0(start_pos, "-", end_pos))
  
  vid <- paste0(sn, ":", pos_str, ":", allele_str)
  return(vid)
}

# ---- constructor ----

#' SemplScores object and constructor
#'
#' Constructs a SemplScores class object.
#'
#' @param variants A `VRanges` object to hold one or more variants
#' @param semData A named list of SNPEffectMatrix objects
#' @param scores (optional) A `data.table` object for motif information and 
#' binding scores
#'
#' @importFrom methods new
#' @importFrom VariantAnnotation VRanges
#' @importFrom S4Vectors mcols
#'
#' @return a SemplScores object
#' @docType class
#' @export
SemplScores <- function(variants=NULL, semData=NULL, scores=NULL) {
  # if no variants provided, make an empty VRanges object
  if (all(is.null(variants))) {
    vr <-  VariantAnnotation::VRanges()
  } else {
    vr <- variants
  }
  
  # if more than one VRange, and there is no id column, make a unique id from
  # position information
  if (length(vr) == 0) {
    S4Vectors::mcols(vr) <- data.frame(id=NA)
  } else if (!("id" %in% names(S4Vectors::mcols(vr)))) {
    vid <- lapply(1:length(vr), \(i) .makeDefaultId(vr = vr[i])) |>
      unlist()
    S4Vectors::mcols(vr)$id <- vid
  }
  
  if (is.null(scores)){
    scores_table <- data.table()
  } else {
    scores_table <- scores
  }
  
  new("SemplScores",
      variants = vr,
      semData = semData,
      scores = scores_table
  )
}


# ---- accessors ----

#' Accessor variants slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' @rdname variants
#' @export
setMethod("variants", "SemplScores", 
          function(x) x@variants)

#' Accessor semData slot in a SemplScores object
#' 
#' @param x a SemplScores object
#' @rdname semData
#' @export
setMethod("semData", "SemplScores", 
          function(x) x@semData)


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
              scores_table[(scores_table$refNorm < 0) & (scores_table$altNorm > 0), ]
            } else if (direction == "lost") {
              scores_table[(scores_table$refNorm > 0) & (scores_table$altNorm < 0), ]
            } else if (direction == "maintained") {
              scores_table[(scores_table$refNorm > 0) & (scores_table$altNorm > 0), ]
            } else if (direction == "changed") {
              scores_table[((scores_table$refNorm < 0) & (scores_table$altNorm > 0)) |
                             ((scores_table$refNorm > 0) & (scores_table$altNorm < 0)), ]
            } else {
              stop("direction is not valid. Options are 'gained', 'lost', 'maintained', or 'changed'")
            }
          }
)

# ---- show ----


#' Show method for SemplScores objects
#'
#' Prints information about the number of variants, SEM meta data columns, and
#' the scoring table if scoreVariants has been run.
#'
#' @param object a SemplScores object
#' 
#' @rdname show
#' 
#' @export
setMethod("show", "SemplScores",
          function(object) {
            cat("An object of class SemplScores\n")
            
            # show variants
            num_vars <- length(object@variants)
            
            if (num_vars > 5) {
              first2 <- object@variants$id[1:2]
              last2 <- object@variants$id[(num_vars-1):num_vars]
              vars_id_list <- paste0(paste(first2, collapse = ", "), " ... ", 
                                     paste(last2, collapse = ", "))
            } else {
              vars_id_list <- object@variants$id |>
                paste0(collapse = ", ")
            }
            
            cat("variants(",num_vars,"): ", sep = "")
            cat(paste(vars_id_list, collapse = " "))
            
            # show semData
            meta_cols <- names(object@semData)
            n_meta_cols <- length(meta_cols)
            if (length(meta_cols) > 5) {
              meta_cols_list <- c(meta_cols[1:2], " ... ",
                                  meta_cols[(n_meta_cols-1):n_meta_cols]) |>
                paste0(collapse = ", ")
            } else {
              meta_cols_list <- paste0(meta_cols, collapse = ", ")
            }

            cat("\nsemData(", n_meta_cols, "): ", meta_cols_list, sep = "")
            
            # show scores
            n_scores <- nrow(object@scores)
            cat("\nscores(", n_scores, "):\n", sep="")
            if(n_scores > 0) { print(object@scores) }
          }
)


