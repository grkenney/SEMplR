.validateGeneType <- \(geneType, orgdb) {
  # fetch every GENETYPE in the orgdb
  valid_types <- AnnotationDbi::keys(orgdb, "GENETYPE")
  
  # if the user's geneType isn't in that set, stop and list the valid ones
  if (!(geneType %in% valid_types)) {
    rlang::abort(c(
      paste0("Invalid geneType '", geneType, "'."),
      "i" = "Run keys(orgdb, 'GENETYPE') to see accepted gene types."
    ))
  }
}


#' Filter Foreground and Background ID Sets by Gene Type
#'
#' @description
#' `poolFilter()` takes the mapped foreground and background ID data frames
#' (as produced by `mapIDs()`) and, if requested, filters both sets
#' to only include genes (or their transcripts) of a specified biotype
#' (e.g. “protein-coding”).
#'
#' @param mapped    A list returned by `mapIDs()`, containing at least:
#'   \itemize{
#'     \item `fg_ids`: data.frame with columns `entrez` and `mappedID` 
#'     (foreground).
#'     \item `bg_ids`: data.frame with columns `entrez` and `mappedID` 
#'     (background).
#'     \item `so_obj`: a `src_organism` object for transcript lookups.
#'     \item `orgdb`:  the loaded OrgDb package object.
#'     \item `transcript`: logical, whether IDs are transcripts.
#'   }
#' @param geneType  Optional character scalar; if not NULL, only genes of this
#'   biotype (`GENETYPE` in the OrgDb) will be retained.  Valid values vary by
#'   organism (e.g. “protein-coding”, “lncRNA”, etc.).
#'
#' @return
#' The original `mapped` list, but with `fg_ids` and `bg_ids` replaced by
#' filtered versions (only rows matching `geneType`, if provided).
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' 
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' orgdb <- org.Hs.eg.db
#' 
#' my_genes <- c("ENSG00000139618", "ENSG00000157764")
#' ids <- mapIDs(orgdb = orgdb, 
#'               foreground_ids = my_genes, 
#'               id_type = "ENSEMBL")
#' filtered <- poolFilter(ids, geneType="protein-coding")
#'
#' @keywords internal
#' @export
poolFilter <- function(mapped,
                       geneType = NULL) {
  # unpack
  fg_df      <- mapped$fg_ids
  bg_df      <- mapped$bg_ids
  transcript <- mapped$transcript
  orgdb      <- mapped$orgdb

  # geneType filtering, only if requested
  if (!is.null(geneType)) {
    .validateGeneType(geneType, orgdb)

    if (transcript) {
      # transcripts that match geneType
      ensembl_trans <- AnnotationDbi::keys(orgdb, "ENSEMBLTRANS")
      ensembl_trans_df <- AnnotationDbi::select(orgdb, 
                                                keys = ensembl_trans,
                                                columns = c("ENTREZID", 
                                                            "GENETYPE"), 
                                                keytype = "ENSEMBLTRANS")
      keep_trans <- ensembl_trans_df[ensembl_trans_df$GENETYPE == geneType, ]
      
      # restrict both bg and fg to transcripts whose gene is in geneType
      bg_df <- bg_df[bg_df$ENTREZID %in% keep_trans$ENTREZID, ]
      fg_df <- fg_df[fg_df$ENTREZID %in% keep_trans$ENTREZID, ]

    } else {
      # gene‐mode: intersect by entrez
      bg_df <- bg_df[ bg_df$GENETYPE == geneType, ]
      fg_df <- fg_df[ fg_df$GENETYPE == geneType, ]
    }
  }

  # reinject filtered dfs and return all components
  mapped$fg_ids <- fg_df
  mapped$bg_ids <- bg_df

  return(mapped)
}


