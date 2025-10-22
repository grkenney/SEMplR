#' Run the Full Differential‐Expression Motif Enrichment Pipeline
#'
#' @description
#' \code{enrichmentSets()} is a one‐stop wrapper that takes a vector of
#' user‐provided IDs from a differential expression analysis and returns
#' promoter regions plus matched background sets for downstream motif
#' enrichment.  It performs:
#' 1. ID mapping via \code{buildMappingObject()} and
#'    \code{mapIDs()}.
#' 2. Optional biotype filtering via \code{poolFilter()}.
#' 3. Promoter coordinate extraction via \code{getCoordinates()}.
#' 4. Background sampling (pool, random, or matched) within that same call.
#'
#' @param txdb A TxDb object. 
#' (e.g. \code{TxDb.Hsapiens.UCSC.hg38.knownGene}).
#' @param orgdb An OrgDb object. (e.g. \code{org.Hs.eg.db}).
#' @param id_type Type of identifier supplied in foreground and background IDs.
#' @param foreground_ids Character vector of gene or
#' transcript IDs (e.g. Ensembl, RefSeq, gene symbols) to analyze.
#' @param background_ids Character vector of gene or transcript 
#' IDs to use as background set.
#' @param transcript Logical; \code{TRUE} to treat inputs as transcript‐level,
#'   \code{FALSE} for gene‐level.
#' @param threshold Numeric in range 0 to 1. Min fraction of IDs that must map 
#' to pick a keytype (default 0.9).
#' @param stripVersions Logical; strip version suffixes (e.g. ".1") from
#' Ensembl/RefSeq IDs.
#' @param inflateThresh Numeric in range 0 to 1; max allowed transcript:gene 
#' inflation before auto‐collapsing (default 1).
#' @param geneType Optional character; biotype filter
#' (e.g. \code{"protein-coding"}).
#' @param overlapMinGap Numeric; minimum gap when reducing
#' overlapping promoters.
#' @param onePromoterPerGene Logical; if
#' \code{TRUE}, pick only one promoter per gene.
#' @param n_ratio Numeric; ratio of ranges to retain in the background set.
#' The number of background ranges will be equal to \code{n_ratio} multiplied
#' by the number of foreground ranges.
#' @param promoterWindow Numeric named vector \code{c(upstream, downstream)};
#'   promoter flank widths in bp (default \code{c(300,50)}).
#' @param standardChroms Logical; restrict to standard chromosomes.
#' @param reduceOverlaps Logical; merge overlapping promoter windows.
#'
#' @return A \code{list} containing the output of
#'   \code{\link{getCoordinates}()}, namely:
#'   \describe{
#'     \item{backgroundElements}{\code{GRanges} of sampled background promoters}
#'     \item{foregroundElements}{\code{GRanges} of foreground promoters}
#'     \item{backgroundUniverse}{\code{GRanges} of the pruned promoter pool}
#'     \item{matchObject}{\code{MatchedGRanges} if \code{bgMethod="matched"},
#'       else \code{NULL}}
#'   }
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' 
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' orgdb <- org.Hs.eg.db
#' 
#' my_genes <- c("ENSG00000139618", "ENSG00000157764")
#' # Minimal run with defaults:
#' results <- enrichmentSets(my_genes, 
#'                           txdb = txdb, 
#'                           orgdb = orgdb, 
#'                           id_type = "ENSEMBL")
#'
#' @export
enrichmentSets <- \(txdb,
                    orgdb,
                    id_type,
                    foreground_ids,
                    background_ids = NULL,
                    transcript            = FALSE,
                    threshold             = 0.9,
                    stripVersions         = TRUE,
                    inflateThresh         = 1,
                    geneType              = NULL,
                    overlapMinGap         = 0,
                    onePromoterPerGene    = FALSE,
                    n_ratio               = 1,
                    promoterWindow        = c(upstream=300,
                                              downstream=50),
                    standardChroms        = TRUE,
                    reduceOverlaps        = TRUE ) {
  
  # Require that orgdb has a GENETYPE column if using geneType param
  if (!is.null(geneType)) {
    od_cols <- AnnotationDbi::columns(orgdb)
    if (!"GENETYPE" %in% od_cols) {
      rlang::abort(paste0(
        "Your OrgDb (", orgdb,
        ") does not contain a 'GENETYPE' column;\n",
        "cannot apply geneType filter '", geneType, "'.\n",
        "Please omit geneType or choose a supported species OrgDb."
      ))
    }
  }
  
  # Map the user's IDs — a bit heavier, but now we know geneType is valid
  mapped <- mapIDs(
    orgdb = orgdb,
    foreground_ids = foreground_ids,
    background_ids = background_ids,
    id_type = id_type,
    threshold = threshold,
    transcript = transcript,
    stripVersions = stripVersions,
    inflateThresh = inflateThresh
  )

  # Pool‐level filtering
  filtered <- poolFilter(
    mapped    = mapped,
    geneType  = geneType
  )

  # Coordinate extraction (lazy until collect, then quick)
  coords <- getCoordinates(
    mapped = filtered,
    txdb = txdb,
    transcript = transcript,
    n_ratio              = n_ratio,
    promoterWindow       = promoterWindow,
    standardChroms       = standardChroms,
    reduceOverlaps       = reduceOverlaps,
    overlapMinGap        = overlapMinGap,
    onePromoterPerGene   = onePromoterPerGene)

  return(coords)

}



