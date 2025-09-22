# Strip Ensembl/RefSeq versions from ids
.stripIds <- \(ids) {
  pat <- "^(ENS[GTPS]\\d+|N[MRP]_\\d+)\\.(\\d+)$"
  idx <- grepl(pat, ids, ignore.case = TRUE)
  ids[idx] <- sub(pat, "\\1", ids[idx], ignore.case = TRUE)
  return(ids)
}


.mapToEntrezIds <- \(id_type, orgdb, ids, threshold, 
                     print_msg = TRUE, transcript = FALSE) {
  if(! id_type %in% keytypes(orgdb)) {
    rlang::abort(c(
      paste0("'", id_type,
             "' is not an available key for in this mapping object."),
      "i" = "Run AnnotationDbi::keytypes(orgdb) to determine valid keytypes."
    ))
  }
  
  mapped_ids <- AnnotationDbi::select(x = orgdb, keys = ids, 
                                      columns = c("ENTREZID", "GENETYPE"), 
                                      keytype = id_type)
  
  pct_mapped <- nrow(mapped_ids) / length(ids)
  
  if (pct_mapped < threshold) {
    rlang::abort(c(paste0("Unable to map >=", threshold*100, "% of your IDs."),
                   "i" = "Ensure the idType provided is correct.",
                   "i" = paste0("Run AnnotationDbi::columns(orgdb) ",
                                "to see available id types.")
    ))
  } else {
    if (print_msg) {
      rlang::inform(paste0( 
        "Successfully mapped ", pct_mapped*100, 
        "% of the provided foreground ids."))
    }
  }
  return(mapped_ids)
}


# What if user specified gene level analysis (transcript = FALSE) but
# foreground_ids provided are for transcripts.
# Conceivable but rare.
#
# Specifically detect most common version of that case (best = ensembltrans &
# transcript = FALSE and warn they aey are losing transcript level coordinate
# specificity with this flag.
#
# Additionally reverse the mapping (entrez -> mapped id) to detect other
# transcript-style ids by  1->many gene to mapped_id inflations and warn
# again.
.checkForInflation <- \(orgdb, ids, id_type, inflateThresh) {
  if (id_type == "ENSEMBLTRANS") {
    rlang::warn(c(paste0(
      "It looks like you provided Ensembl transcript IDs ('", id_type,
      "') but requested gene-level analysis (transcript = FALSE)."),
      "i" = paste0("Any downstream coordinate lookup will use gene (Entrez) ",
                   "IDs, so you'll lose the per-transcript specificity of ",
                   "your input."),
      "i" = paste0("If you really want transcript-level coordinates, set",
                   "transcript = TRUE or supply Ensembl gene IDs instead.")
    ))
    
  } else {
    # test for other instances of 1:many gene:foreground_id inflations
    # reverse the mapping from the fg data
    rev_df <- AnnotationDbi::select(orgdb, 
                                    keytype = "ENTREZID", 
                                    columns = id_type, 
                                    keys = ids$ENTREZID)
    
    # count genes, and mapped_ids and look for excessive inflation
    nGenes   <- unique(rev_df[, "ENTREZID"]) |> length()
    nMapped  <- unique(rev_df[, id_type]) |> length()
    if (nGenes > 0) {
      inflation <- nMapped / nGenes - 1
      if (inflation > inflateThresh) {
        rlang::warn(c(
          "Your IDs appear transcript-like: ",
          paste0(
            "Reverse-mapping shows ", inflation * 100, "% more ",
            "(", nMapped, " unique transcript IDs for ", nGenes, " genes). "),
          "Downstream, only gene-level coordinates will be used.",
          paste0("If you need transcript-level analyses, set, ",
                 "transcript = TRUE and use Ensembl transcript IDs.")
        ))
      }
    }
  }
}


#' Map User IDs to Entrez and Determine Best Keytype
#'
#' Given a vector of user-supplied gene/transcript IDs, finds the AnnotationDbi
#' keytype (e.g. ens, refseq, symbol, etc.) that maps the highest fraction of
#' inputs, and returns both foreground and full background sets as Entrez IDs.
#' Optionally collapses transcript-style inputs to genes when requested or when
#' reverse-mapping inflation exceeds a threshold.
#'
#' @param foreground_ids Character vector of input IDs (e.g. ENSG, ENST, NM_*)
#' @param mapping        A list from \code{buildMappingObject()}, containing:
#'   \itemize{
#'     \item \code{so_obj}: the \code{src_organism} object
#'     \item \code{orgdb}:  the loaded OrgDb package object
#'     \item \code{organism}, \code{genomeBuild}, \code{txdb}: parameters used
#'   }
#' @param threshold      Fraction in range 0 to 1; minimum mapping rate to 
#' accept a keytype without falling back (default 0.9).
#' @param transcript     Logical; if \code{TRUE}, analyze as transcript-level 
#'   IDs (default \code{FALSE}).
#' @param stripVersions  Logical; strip trailing ".1", ".2" from IDs 
#' (default \code{TRUE}).
#' @param inflateThresh  Fraction in range 0 to 1; if reverse-mapping shows 
#' excessive, inflation automatically collapse transcripts to genes 
#' (default 1 ie. 100%).
#'
#' @return
#' A named list combining the original \code{mapping} components with:
#' \describe{
#'  \item{\code{fg_ids}}{data.frame(entrez, mappedID) for your foreground set}
#'  \item{\code{bg_ids}}{data.frame(entrez, mappedID) for the full background}
#'  \item{\code{userIDtype}}{the chosen keytype (e.g. "ensembl")}
#'  \item{\code{transcript}}{logical, whether transcript-level mapping was used}
#' }
#'
#' @examples
#'   # Gene Ids
#'   mapping <- buildMappingObject("Homo sapiens")
#'   ids <- mapIDs(
#'     mapping = mapping,
#'     foreground_ids = c("ENSG00000139618","ENSG00000157764")
#'   )
#'   ids
#'   
#'   # Transcript Ids
#'   mapping <- buildMappingObject("Homo sapiens")
#'   ids <- mapIDs(
#'     mapping = mapping,
#'     foreground_ids = c("ENST00000245479", "ENST00000633194"),
#'     transcript = TRUE
#'   )
#'   ids
#'
#' @importFrom AnnotationDbi columns keytypes
#' @export
mapIDs <- \(orgdb,                
            foreground_ids,
            background_ids = NULL,
            id_type,
            threshold = 0.9,     
            transcript = FALSE,   
            stripVersions = TRUE,    
            inflateThresh = 1) {
  
  if (stripVersions) {
    foreground_ids <- .stripIds(foreground_ids)
    background_ids <- .stripIds(background_ids)
  }
  
  rlang::inform("Mapping foreground ids to ENTREZIDs...")
  fg_id_map <- .mapToEntrezIds(id_type = id_type, 
                               orgdb = orgdb, 
                               ids = foreground_ids, 
                               threshold = threshold)

  # Background pool can be restricted according to background_ids (if provided)
  if(is.null(background_ids)) {
    rlang::inform("Building background id set...")
    # Otherwise background pool is all records in orgdb
    background_ids <- keys(orgdb, keytype = id_type)
    print_msg <- FALSE
  } else {
    rlang::inform("Mapping background ids to ENTREZIDs...")
    print_msg <- TRUE
  }
  
  bg_id_map <- .mapToEntrezIds(id_type = id_type, orgdb = orgdb, 
                               ids = background_ids, 
                               threshold = threshold,
                               print_msg = print_msg)

  # Ensure background pool is the same universe as foreground by dropping any
  # rows with no available value for best mappedID type
  bg_na_filter <- is.na(bg_id_map[, id_type])
  bg_id_map <- bg_id_map[!bg_na_filter, ]
  
  if (!transcript) {
    rlang::inform("Checking for inflation...")
    .checkForInflation(orgdb = orgdb, 
                       ids = fg_id_map, 
                       id_type = id_type, 
                       inflateThresh = inflateThresh)
  }

  # Report mapping statistics:
  # Make output list
  mapped <- c(
    orgdb = orgdb,
    list(
      fg_ids = fg_id_map,
      bg_ids = bg_id_map,
      userIDtype = id_type,
      transcript = transcript
    )
  )

  return(mapped)
}
