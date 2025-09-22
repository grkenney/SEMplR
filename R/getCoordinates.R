.getTranscriptRanges <- \(txdb, fg_df, bg_df, standardChroms) {
  columns_to_keep <- c("GENEID", "TXID", "TXNAME")
  # get all transcripts
  txs <- transcripts(txdb, columns = columns_to_keep)
  
  # get ranges for foreground
  fg_ix <- lapply(fg_df$ENSEMBLTRANS, 
                  function(x) grep(pattern = x, x = txs$tx_name)) |>
    unlist()
  fg_ranges <- txs[fg_ix, ]
  
  # get ranges for background
  bg_ix <- lapply(bg_df$ENSEMBLTRANS, 
                  function(x) grep(pattern = x, x = txs$tx_name)) |>
    unlist()
  bg_ranges <- txs[bg_ix, ]
  
  organism <- metadata(txdb)[metadata(txdb)$name == "Organism", "value"]
  # Restrict to standard chromosomes? Optional but default yes
  if(standardChroms){
    bg_ranges <- GenomeInfoDb::keepStandardChromosomes(bg_ranges,
                                                       species = organism, 
                                                       pruning.mode = "coarse")
  }
  # Generate foreground elements granges by subsetting bg_gr by mappedID
  return(c(fg = fg_ranges, bg = bg_ranges))
}


.getGeneRanges <- \(txdb, fg_df, bg_df, standardChroms) {
  columns_to_keep <- "GENEID"
  all_genes <- GenomicFeatures::genes(txdb, columns = columns_to_keep) |> 
    suppressMessages()
  
  # find indices of ENTREZIDs in all_genes
  fg_ix <- match(fg_df$ENTREZID, unlist(all_genes$GENEID))
  bg_ix <- match(bg_df$ENTREZID, unlist(all_genes$GENEID))
  
  # remove non-matching genes from list
  fg_df_in_genes <- fg_df[!is.na(fg_ix), ]
  bg_df_in_genes <- bg_df[!is.na(bg_ix), ]
  
  fg_ranges <- all_genes[na.omit(fg_ix)] # ignore NAs in indices now
  bg_ranges <- all_genes[na.omit(bg_ix)] # ignore NAs in indices now
  
  # add metadata
  S4Vectors::mcols(fg_ranges) <- fg_df_in_genes
  S4Vectors::mcols(bg_ranges) <- bg_df_in_genes
  
  return(c(fg = fg_ranges, bg = bg_ranges))
}


.extractPromoters <- \(gr, promoterWindow, transcript,
                       reduceOverlaps, onePromoterPerGene) {
  upstream <- promoterWindow[["upstream"]]
  downstream <- promoterWindow[["downstream"]]
  
  prom_bg <- GenomicRanges::promoters(
    gr$bg,
    upstream   = upstream,
    downstream = downstream,
    use.names  = TRUE
  )
  prom_fg <- GenomicRanges::promoters(
    gr$fg,
    upstream   = upstream,
    downstream = downstream,
    use.names  = TRUE
  )

  # (optional) merge any overlapping promoter windows within each set
  if (reduceOverlaps) {
    rlang::inform(paste0("Combining overlapping promoter ranges within genes.",
    " (This step may take 1-2 minutes)..."))
    prom_fg <- .reduceOverlapsWithinGenes(prom_fg)
    prom_bg <- .reduceOverlapsWithinGenes(prom_bg)
  }

  # (optional) enforce one promoter per gene: choose the widest window
  if (onePromoterPerGene) {
    rlang::inform("Restricting to one promoter range per gene...")
    prom_bg <- .subsetToOneTxPerGene(prom_bg)
    prom_fg <- .subsetToOneTxPerGene(prom_fg)
  }
  return(c(fg = prom_fg, bg = prom_bg))
}


# Reduce overlapping ranges within each gene
.reduceOverlapsWithinGenes <- function(gr) {
  # reduce within gene
  gr <- lapply(unique(gr$ENTREZID), 
               function(x) 
                 .reduceGene(gr[gr$ENTREZID == x])) |>
    GenomicRanges::GRangesList() |>
    unlist()
  
  return(gr)
}


.reduceGene <- \(foo_gene) {
  if (length(foo_gene) > 1) {
    mcol_foo_gene <- S4Vectors::mcols(foo_gene)
    foo_gene_red <- GenomicRanges::reduce(foo_gene, with.revmap = TRUE)
    for (cn in colnames(mcol_foo_gene)) {
      S4Vectors::mcols(foo_gene_red)[cn] <-
        sapply(foo_gene_red$revmap,
               function(i) paste(unique(mcol_foo_gene[i, cn]), collapse = ", "))
    }
    foo_gene_red$revmap <- NULL
    return(foo_gene_red)
  } else {
    return(foo_gene)
  }
}


# Select one transcript per gene (widest, then 5'-most)
.subsetToOneTxPerGene <- function(gr) {
  gr <- lapply(unique(gr$ENTREZID), 
                    function(x) 
                      .findWidestOrMost5PrimeRange(gr[gr$ENTREZID == x])) |>
    GenomicRanges::GRangesList() |>
    unlist()
  return(gr)
}


.findWidestOrMost5PrimeRange <- \(gr) {
  # if only one grange, don't need to subset
  if (length(gr) == 1) {
    return(gr)
  } else {
    # if multiple grange, first try to subset by width
    gr_widths <- GenomicRanges::width(gr)
    max_width_grs <- gr[gr_widths == max(gr_widths)]
    
    # if multiple with largest width, subset to most 5'
    if (length(max_width_grs) > 1) {
      gr_strand <- GenomicRanges::strand(max_width_grs)[1]
      if (as.character(gr_strand) == "+") {
        ix <- which.min(GenomicRanges::start(max_width_grs))
        return(max_width_grs[gr_i])
      } else {
        ix <- which.max(GenomicRanges::end(max_width_grs))
        return(max_width_grs[gr_i])
      }
    } else {
      return(max_width_grs)
    }
  }
}


# Define and Sample Background Elements for Motif Enrichment
.defineBackgroundElements <- function(
    background_universe,
    foreground_elements,
    n_ratio) {

  ## Pruning steps
  # Remove all background ranges from pool if any overlap with
  # foreground range
  background_universe <- IRanges::subsetByOverlaps(
    background_universe, foreground_elements, invert = TRUE
    )

  # Remove all background genes from pool if they appear in foreground
  fg_genes <- unique(S4Vectors::mcols(foreground_elements)["ENTREZID"])
  background_universe <- background_universe[
    ! S4Vectors::mcols(background_universe)$ENTREZID %in% fg_genes$ENTREZID
  ]
  
  selectedBg <- .randomBackground(pool = background_universe,
                                  focal = foreground_elements,
                                  n_ratio = n_ratio)

  bg_gr       <- selectedBg
  fg_gr       <- foreground_elements
  universe    <- background_universe

  # Return a consistent list object
  out <- list(
    backgroundElements = bg_gr,
    foregroundElements = fg_gr,
    backgroundUniverse = universe
  )
  return(out)
}


# Randomly sample background regions
.randomBackground <- \(pool, focal, n_ratio) {
  n_fg   <- length(focal)
  n_pool <- length(pool)
  n_bg   <- n_ratio * n_fg
  
  if (n_bg > n_pool) {
    stop(
      "Requested ", n_bg, " background regions (",
      n_ratio, "x", n_fg, ") but only ",
      n_pool, " available in the pool."
    )
  }
  
  idx <- sample(
    seq_len(n_pool),
    size    = n_bg
  )
  out <- pool[idx]
  return(out)
}


#' Retrieve Promoter Regions and Sample Background Elements
#'
#' @description
#' Given a mapped set of foreground IDs (with Entrez and mappedID),
#' this function:
#' 1. Fetches genomic coordinates for each feature (gene or transcript),
#'    according to the chosen \code{TSS.method}.
#' 2. Extracts promoter windows around those coordinates.
#' 3. Optionally reduces overlaps and selects one promoter per gene.
#' 4. Samples background promoter elements matching the foreground
#'
#' @param mapped A list as returned by \code{\link{mapIDs}()}, containing
#'   at least \code{fg_ids}, \code{bg_ids}, \code{so_obj}, \code{orgdb},
#'   \code{genomeBuild}, and \code{organism}.
#' @param transcript Logical; \code{TRUE} for transcript-level coordinates,
#'   \code{FALSE} for gene-level.
#' @param TSS.method Character(1). Method to pick TSS when
#' \code{transcript=FALSE}:
#'   \code{"UCSCgene"}, \code{"Ensembl_canonical"}, \code{"commonTSS"},
#'   \code{"uniqueTSS"}, \code{"fivePrimeTSS"}, or \code{"allTSS"}.
#' @param bgMethod Character; one of \code{"pool"}, \code{"random"}, or
#'   \code{"matched"}.
#' @param n_ratio Numeric; for \code{bgMethod="random"}, number of backgrounds =
#'   \code{n_ratio * length(foreground)}.
#' @param bgExcludeFgOverlaps Logical; if \code{TRUE}, drop any background whose
#'   promoter overlaps a foreground.
#' @param bgExcludeFgGenes Logical; if \code{TRUE}, drop any background gene
#'   whose EntrezID is in the foreground.
#' @param covariates Character vector of covariate names for matching:
#'   \code{"width"}, \code{"gc"}, \code{"seqnames"} (or \code{"chromosome"}), or
#'   any user-provided \code{mcols()} column.
#' @param bgReplace Logical; allow replacement when sampling backgrounds.
#' @param nrMethod Character; \code{"rejection"}, \code{"nearest"}, or
#'   \code{"stratified"}.
#' @param genome A \pkg{genome} object. If \code{NULL} and \code{"gc"} is in
#'   \code{covariates}, will be auto-resolved.
#' @param promoterWindow Numeric named vector of lengths:
#'   \code{c(upstream, downstream)} (default \code{c(300,50)}).
#' @param standardChroms Logical; restrict to standard chromosomes.
#' @param reduceOverlaps Logical; merge any overlapping promoter windows.
#' @param overlapMinGap Numeric; minimum gap when reducing overlaps.
#' @param onePromoterPerGene Logical; if \code{TRUE}, choose one
#' promoter per gene.
#' @param ensdb Optional \code{EnsDb} object; required for
#'   \code{TSS.method="Ensembl_canonical"}.
#'
#' @return A named \code{list} from the final
#' \code{.defineBackgroundElements()}:
#'   \item{backgroundElements}{\code{GRanges} of sampled background promoters}
#'   \item{foregroundElements}{\code{GRanges} of foreground promoters}
#'   \item{backgroundUniverse}{\code{GRanges} of the pruned universe}
#'   \item{matchObject}{\code{MatchedGRanges} if \code{bgMethod="matched"},
#'   else \code{NULL}}
#'
#' @examples
#' my_genes <- c("ENSG00000139618", "ENSG00000157764")
#' mapping <- buildMappingObject("Homo sapiens")
#' ids <- mapIDs(mapping = mapping, foreground_ids = my_genes)
#' filtered <- poolFilter(mapped, geneType="protein-coding")
#' coords <- getCoordinates(
#'   filtered,
#'   transcript     = FALSE,
#'   TSS.method     = "UCSCgene",
#'   bgMethod       = "matched",
#'   n_ratio        = 1,
#'   covariates     = c("width","gc","chromosome"),
#'   nrMethod       = "stratified",
#'   genome         = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'   promoterWindow = c(upstream=300, downstream=50)
#' )
#'
#' @export
getCoordinates <- function(mapped,
                           transcript = FALSE,
                           n_ratio = 1,
                           promoterWindow = c(upstream=300,
                                              downstream=50),
                           standardChroms = TRUE,
                           reduceOverlaps = TRUE,
                           overlapMinGap = 0,
                           onePromoterPerGene = FALSE) {
  
  if (transcript && (reduceOverlaps || onePromoterPerGene)) {
    rlang::abort(paste0(
      "reduceOverlaps=TRUE and/or onePromoterPerGene=TRUE only ",
      "makes sense for gene-level analysis (transcript=FALSE). ",
      "Please rerun with transcript=FALSE to use those options."
    ))
  }

  # unpack
  fg_df       <- mapped$fg_ids
  bg_df       <- mapped$bg_ids
  transcript  <- mapped$transcript

  if (transcript) {
    gr <- .getTranscriptRanges(txdb, fg_df, bg_df, standardChroms)
    } else {
    gr <- .getGeneRanges(txdb, fg_df, bg_df, standardChroms)
    }
  
  #### ---- promoter extraction ----
  proms <- .extractPromoters(gr = gr,
                             promoterWindow = promoterWindow,
                             transcript = transcript,
                             reduceOverlaps = reduceOverlaps,
                             onePromoterPerGene = onePromoterPerGene)

  #### ---- background element selection ----
  rlang::inform("Defining background elements...")
  out <- .defineBackgroundElements(background_universe = proms$bg,
                                  foreground_elements = proms$fg,
                                  n_ratio = n_ratio)

  return(out)

}
