# format a list of items to print well in show function
.formatList <- function(x) {
    num_items <- length(x)
    if (num_items > 5) {
        first2 <- x[seq_len(2)]
        last2 <- x[(num_items - 1):num_items]
        formatted_list <- paste0(
            paste(first2, collapse = ", "), " ... ",
            paste(last2, collapse = ", ")
        )
    } else {
        formatted_list <- x |>
            paste0(collapse = ", ")
    }
    return(formatted_list)
}


# most functions in this package expect a SNPEffectMatrixCollection
# for convenience, this functions converts lists of SNPEffectMatrix's,
# or single SNPEffectMatrix into SNPEffectMatrixCollections
.convertToSNPEffectMatrixCollection <- function(x) {
    # if it's already a collection, return as is
    # if given a list of SNPEffectMatrices or a single SNPEffectMatrix,
    # make a collection
    # else, fail if unable to convert
    if (is(x, "SNPEffectMatrixCollection")) {
        return(x)
    } else if (is(x, "list") | is(x, "vector")) {
        # check that all elements are SNPEffectMatrices
        class_check <- lapply(x, function(y) {
            is(y, class2 = "SNPEffectMatrix")
        }) |>
            unlist()
        if (all(class_check)) {
            return(SNPEffectMatrixCollection(x))
        } else {
            invalid_class <- is(x[!class_check][1])[1]
            rlang::abort(paste0(
                "unable to convert object of class ", invalid_class,
                " to class SNPEffectMatrixCollection.\n",
                "See ?SNPEffectMatrixCollection or use the provided ",
                "default 'sc'"
            ))
        }
    } else if (is(x, "SNPEffectMatrix")) {
        return(SNPEffectMatrixCollection(sems = x))
    } else {
        rlang::abort(paste0(
            "unable to convert object of class ", is(x)[1],
            " to class SNPEffectMatrixCollection.\n",
            "See ?SNPEffectMatrixCollection or use the provided default 'sc'"
        ))
    }
}


# Given a single VRange or GRange, construct a unique id from the
# position and allele information
.makeVariantId <- function(x, refCol = NULL, altCol = NULL) {
    start_pos <- IRanges::start(IRanges::ranges(x))
    end_pos <- IRanges::end(IRanges::ranges(x))
    sn <- GenomeInfoDb::seqnames(x)
    
    if (is(x, "VRanges")) {
        ref_allele <- as.character(VariantAnnotation::ref(x))
        alt_allele <- as.character(VariantAnnotation::alt(x))
    } else {
        if (is.null(refCol) | is.null(altCol)) {
            rlang::abort(paste0(
                "If providing a GRanges object, ",
                "both refCol and altCol must be defined"
            ))
        }
        ref_allele <- as.character(S4Vectors::mcols(x)[, refCol])
        alt_allele <- as.character(S4Vectors::mcols(x)[, altCol])
    }
    
    if (ref_allele == "") {
        allele_str <- paste0("ins", alt_allele)
    } else if (alt_allele == "") {
        allele_str <- paste0("del", ref_allele)
    } else {
        allele_str <- paste0(ref_allele, ">", alt_allele)
    }
    
    pos_str <- ifelse(start_pos == end_pos,
                      start_pos, paste0(start_pos, "-", end_pos)
    )
    
    vid <- paste0(sn, ":", pos_str, ":", allele_str)
    return(vid)
}


# convert a SNP Effect Matrix to a Position Probability Matrix
.semToPpm <- function(s) {
    # normalize matrix
    norm_score <- apply(
        getSEM(s), 1,
        function(x) {
            (2^x - 2^getBaseline(s)) /
                abs(2^getBaseline(s))
        }
    )
    # replace negative scores with zero
    norm_score[norm_score < 0] <- 0
    # make all rows sum to 1
    ppm <- apply(norm_score, 2, function(x) x / sum(x))
    return(ppm)
}


.getTranscriptRanges <- function(txdb, fg_df, bg_df, standardChroms) {
    columns_to_keep <- c("GENEID", "TXID", "TXNAME")
    # get all transcripts
    txs <- GenomicFeatures::transcripts(txdb, columns = columns_to_keep)
    
    # get ranges for foreground
    fg_ix <- lapply(
        fg_df$ENSEMBLTRANS,
        function(x) grep(pattern = x, x = txs$tx_name)
    ) |>
        unlist()
    fg_ranges <- txs[fg_ix, ]
    
    # get ranges for background
    bg_ix <- lapply(
        bg_df$ENSEMBLTRANS,
        function(x) grep(pattern = x, x = txs$tx_name)
    ) |>
        unlist()
    bg_ranges <- txs[bg_ix, ]
    
    organism <- S4Vectors::metadata(txdb)[S4Vectors::metadata(txdb)$name ==
                                              "Organism", "value"]
    # Restrict to standard chromosomes? Optional but default yes
    if (standardChroms) {
        bg_ranges <- GenomeInfoDb::keepStandardChromosomes(bg_ranges,
                                                           species = organism,
                                                           pruning.mode = 
                                                               "coarse"
        )
    }
    # Generate foreground elements granges by subsetting bg_gr by mappedID
    return(c(fg = fg_ranges, bg = bg_ranges))
}


.getGeneRanges <- function(txdb, fg_df, bg_df, standardChroms) {
    columns_to_keep <- "GENEID"
    all_genes <- GenomicFeatures::genes(txdb, columns = columns_to_keep) |>
        suppressMessages()
    
    # find indices of ENTREZIDs in all_genes
    fg_ix <- match(fg_df$ENTREZID, unlist(all_genes$GENEID))
    bg_ix <- match(bg_df$ENTREZID, unlist(all_genes$GENEID))
    
    # remove non-matching genes from list
    fg_df_in_genes <- fg_df[!is.na(fg_ix), ]
    bg_df_in_genes <- bg_df[!is.na(bg_ix), ]
    
    fg_ranges <- all_genes[stats::na.omit(fg_ix)] # ignore NAs in indices now
    bg_ranges <- all_genes[stats::na.omit(bg_ix)] # ignore NAs in indices now
    
    # add metadata
    S4Vectors::mcols(fg_ranges) <- fg_df_in_genes
    S4Vectors::mcols(bg_ranges) <- bg_df_in_genes
    
    return(c(fg = fg_ranges, bg = bg_ranges))
}


.extractPromoters <- function(gr, promoterWindow, transcript,
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
        rlang::inform(paste0(
            "Combining overlapping promoter ranges within genes.",
            " (This step may take 1-2 minutes)..."
        ))
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
    gr <- lapply(
        unique(gr$ENTREZID),
        function(x) {
            .reduceGene(gr[gr$ENTREZID == x])
        }
    ) |>
        GenomicRanges::GRangesList() |>
        unlist()
    
    return(gr)
}


# collapse overlapping ranges while preseving meta data
.reduceGene <- function(gene_ranges) {
    if (length(gene_ranges) > 1) {
        mcol_foo_gene <- S4Vectors::mcols(gene_ranges)
        foo_gene_red <- GenomicRanges::reduce(gene_ranges, with.revmap = TRUE)
        for (cn in colnames(mcol_foo_gene)) {
            S4Vectors::mcols(foo_gene_red)[cn] <-
                vapply(
                    foo_gene_red$revmap,
                    function(i) {
                        paste(unique(mcol_foo_gene[i, cn]), collapse = ", ")
                    }, ""
                )
        }
        foo_gene_red$revmap <- NULL
        return(foo_gene_red)
    } else {
        return(gene_ranges)
    }
}


# Select one transcript per gene (widest, then 5'-most)
.subsetToOneTxPerGene <- function(gr) {
    gr <- lapply(
        unique(gr$ENTREZID),
        function(x) {
            .findWidestOrMost5PrimeRange(gr[gr$ENTREZID == x])
        }
    ) |>
        GenomicRanges::GRangesList() |>
        unlist()
    return(gr)
}


.findWidestOrMost5PrimeRange <- function(gr) {
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
                return(max_width_grs[ix])
            } else {
                ix <- which.max(GenomicRanges::end(max_width_grs))
                return(max_width_grs[ix])
            }
        } else {
            return(max_width_grs)
        }
    }
}


# Define and Sample Background Elements for Motif Enrichment
.defineBackgroundElements <- function(background_universe,
                                      foreground_elements,
                                      n_ratio) {
    ## Pruning steps
    # Remove all background ranges from pool if any overlap with
    # foreground range
    background_universe <- IRanges::subsetByOverlaps(
        background_universe, foreground_elements,
        invert = TRUE
    )
    
    # Remove all background genes from pool if they appear in foreground
    fg_genes <- unique(S4Vectors::mcols(foreground_elements)["ENTREZID"])
    background_universe <- background_universe[
        !S4Vectors::mcols(background_universe)$ENTREZID %in% fg_genes$ENTREZID
    ]
    
    selectedBg <- .randomBackground(
        pool = background_universe,
        focal = foreground_elements,
        n_ratio = n_ratio
    )
    
    bg_gr <- selectedBg
    fg_gr <- foreground_elements
    universe <- background_universe
    
    # Return a consistent list object
    out <- list(
        backgroundElements = bg_gr,
        foregroundElements = fg_gr,
        backgroundUniverse = universe
    )
    return(out)
}


# Randomly sample background regions
.randomBackground <- function(pool, focal, n_ratio) {
    n_fg <- length(focal)
    n_pool <- length(pool)
    n_bg <- n_ratio * n_fg
    
    if (n_bg > n_pool) {
        rlang::abort(paste0(
            "Requested ", n_bg, " background regions (",
            n_ratio, "x", n_fg, ") but only ",
            n_pool, " available in the pool."
        ))
    }
    
    idx <- sample(
        seq_len(n_pool),
        size = n_bg
    )
    out <- pool[idx]
    return(out)
}


# Retrieve Promoter Regions and Sample Background Elements
#
# @description
# Given a mapped set of foreground IDs (with Entrez and mappedID),
# this function:
# 1. Fetches genomic coordinates for each feature (gene or transcript),
#    according to the chosen \code{TSS.method}.
# 2. Extracts promoter windows around those coordinates.
# 3. Optionally reduces overlaps and selects one promoter per gene.
# 4. Samples background promoter elements matching the foreground
#
# @param mapped A list as returned by \code{\link{mapIDs}()}, containing
#   at least \code{fg_ids} and \code{bg_ids}.
# @param txdb A TxDb object.
# (e.g. \code{TxDb.Hsapiens.UCSC.hg38.knownGene}).
# @param transcript Logical; \code{TRUE} for transcript-level coordinates,
#   \code{FALSE} for gene-level.
# @param n_ratio Numeric; ratio of ranges to retain in the background set.
# The number of background ranges will be equal to \code{n_ratio} multiplied
# by the number of foreground ranges.
# @param promoterWindow Numeric named vector of lengths:
#   \code{c(upstream, downstream)} (default \code{c(300,50)}).
# @param standardChroms Logical; restrict to standard chromosomes.
# @param reduceOverlaps Logical; merge any overlapping promoter windows.
# @param overlapMinGap Numeric; minimum gap when reducing overlaps.
# @param onePromoterPerGene Logical; if \code{TRUE}, choose one
# promoter per gene.
#
# @return A named \code{list} from the final
# \code{.defineBackgroundElements()}:
#   \item{backgroundElements}{\code{GRanges} of sampled background promoters}
#   \item{foregroundElements}{\code{GRanges} of foreground promoters}
#   \item{backgroundUniverse}{\code{GRanges} of the pruned universe}
#   \item{matchObject}{\code{MatchedGRanges} if \code{bgMethod="matched"},
#   else \code{NULL}}
#
# @examples
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
#
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# orgdb <- org.Hs.eg.db
#
# my_genes <- c("ENSG00000139618", "ENSG00000157764")
# ids <- mapIDs(
#     orgdb = orgdb,
#     foreground_ids = my_genes,
#     id_type = "ENSEMBL"
# )
# filtered <- poolFilter(ids, geneType = "protein-coding")
# coords <- getCoordinates(mapped = filtered, txdb = txdb)
#
# @export
.getCoordinates <- function(mapped,
                            txdb,
                            transcript = FALSE,
                            n_ratio = 1,
                            promoterWindow = c(
                                upstream = 300,
                                downstream = 50
                            ),
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
    fg_df <- mapped$fg_ids
    bg_df <- mapped$bg_ids
    
    if (transcript) {
        gr <- .getTranscriptRanges(txdb, fg_df, bg_df, standardChroms)
    } else {
        gr <- .getGeneRanges(txdb, fg_df, bg_df, standardChroms)
    }
    
    #### ---- promoter extraction ----
    proms <- .extractPromoters(
        gr = gr,
        promoterWindow = promoterWindow,
        transcript = transcript,
        reduceOverlaps = reduceOverlaps,
        onePromoterPerGene = onePromoterPerGene
    )
    
    #### ---- background element selection ----
    rlang::inform("Defining background elements...")
    out <- .defineBackgroundElements(
        background_universe = proms$bg,
        foreground_elements = proms$fg,
        n_ratio = n_ratio
    )
    
    return(out)
}


# Strip Ensembl/RefSeq versions from ids
.stripIds <- function(ids) {
    pat <- "^(ENS[GTPS]\\d+|N[MRP]_\\d+)\\.(\\d+)$"
    idx <- grepl(pat, ids, ignore.case = TRUE)
    ids[idx] <- sub(pat, "\\1", ids[idx], ignore.case = TRUE)
    return(ids)
}


.mapToEntrezIds <- function(id_type, orgdb, ids, threshold,
                            print_msg = TRUE, transcript = FALSE) {
    if (!id_type %in% AnnotationDbi::keytypes(orgdb)) {
        rlang::abort(c(
            paste0(
                "'", id_type,
                "' is not an available key for in this mapping object."
            ),
            "i" = paste0(
                "Run AnnotationDbi::keytypes(orgdb) to determine",
                " valid keytypes."
            )
        ))
    }
    
    mapped_ids <- AnnotationDbi::select(
        x = orgdb, keys = ids,
        columns = c("ENTREZID", "GENETYPE"),
        keytype = id_type
    )
    
    pct_mapped <- nrow(mapped_ids) / length(ids)
    
    if (pct_mapped < threshold) {
        rlang::abort(c(
            paste0(
                "Unable to map >=",
                threshold * 100, "% of your IDs."
            ),
            "i" = "Ensure the idType provided is correct.",
            "i" = paste0(
                "Run AnnotationDbi::columns(orgdb) ",
                "to see available id types."
            )
        ))
    } else {
        if (print_msg) {
            rlang::inform(paste0(
                "Successfully mapped ", pct_mapped * 100,
                "% of the provided foreground ids."
            ))
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
.checkForInflation <- function(orgdb, ids, id_type, inflateThresh) {
    if (id_type == "ENSEMBLTRANS") {
        rlang::warn(c(
            paste0(
                "It looks like you provided Ensembl transcript IDs ('", id_type,
                "') but requested gene-level analysis (transcript = FALSE)."
            ),
            "i" = paste0(
                "Any downstream coordinate lookup will use gene (Entrez) ",
                "IDs, so you'll lose the per-transcript specificity of ",
                "your input."
            ),
            "i" = paste0(
                "If you really want transcript-level coordinates, set",
                "transcript = TRUE or supply Ensembl gene IDs instead."
            )
        ))
    } else {
        # test for other instances of 1:many gene:foreground_id inflations
        # reverse the mapping from the fg data
        rev_df <- AnnotationDbi::select(orgdb,
                                        keytype = "ENTREZID",
                                        columns = id_type,
                                        keys = ids$ENTREZID
        )
        
        # count genes, and mapped_ids and look for excessive inflation
        nGenes <- unique(rev_df[, "ENTREZID"]) |> length()
        nMapped <- unique(rev_df[, id_type]) |> length()
        if (nGenes > 0) {
            inflation <- nMapped / nGenes - 1
            if (inflation > inflateThresh) {
                rlang::warn(c(
                    "Your IDs appear transcript-like: ",
                    paste0(
                        "Reverse-mapping shows ", inflation * 100, "% more ",
                        "(", nMapped, " unique transcript IDs for ",
                        nGenes, " genes). "
                    ),
                    "Downstream, only gene-level coordinates will be used.",
                    paste0(
                        "If you need transcript-level analyses, set, ",
                        "transcript = TRUE and use Ensembl transcript IDs."
                    )
                ))
            }
        }
    }
}


# Map User IDs to Entrez and Determine Best Keytype
#
# Given a vector of user-supplied gene/transcript IDs, finds the AnnotationDbi
# keytype (e.g. ens, refseq, symbol, etc.) that maps the highest fraction of
# inputs, and returns both foreground and full background sets as Entrez IDs.
# Optionally collapses transcript-style inputs to genes when requested or when
# reverse-mapping inflation exceeds a threshold.
#
# @param orgdb An OrgDb object. (e.g. \code{org.Hs.eg.db}).
# @param id_type Type of identifier supplied in foreground and background IDs.
# @param foreground_ids Character vector of gene or
# transcript IDs (e.g. Ensembl, RefSeq, gene symbols) to analyze.
# @param background_ids Character vector of gene or transcript
# IDs to use as background set.
# @param threshold Fraction in range 0 to 1; minimum mapping rate to
# accept a keytype without falling back (default 0.9).
# @param transcript Logical; if \code{TRUE}, analyze as transcript-level
#   IDs (default \code{FALSE}).
# @param stripVersions Logical; strip version suffixes (e.g. ".1") from
# Ensembl/RefSeq IDs.
# @param inflateThresh  Fraction in range 0 to 1; if reverse-mapping shows
# excessive, inflation automatically collapse transcripts to genes
# (default 1 ie. 100%).
#
# @return
# A named list combining the original \code{mapping} components with:
# \describe{
#  \item{\code{fg_ids}}{data.frame(entrez, mappedID) for your foreground set}
#  \item{\code{bg_ids}}{data.frame(entrez, mappedID) for the full background}
#  \item{\code{userIDtype}}{the chosen keytype (e.g. "ensembl")}
#  \item{\code{transcript}}{logical, whether transcript-level mapping was used}
# }
#
# @examples
# library(org.Hs.eg.db)
# orgdb <- org.Hs.eg.db
#
# my_genes <- c("ENSG00000139618", "ENSG00000157764")
# ids <- mapIDs(
#     orgdb = orgdb,
#     foreground_ids = my_genes,
#     id_type = "ENSEMBL"
# )
#
# # Transcript Ids
# my_transcripts <- c("ENST00000245479", "ENST00000633194")
# ids <- mapIDs(
#     orgdb = orgdb,
#     foreground_ids = my_transcripts,
#     id_type = "ENSEMBLTRANS",
#     transcript = TRUE
# )
#
# @importFrom AnnotationDbi columns keytypes
# @export
.mapIDs <- function(orgdb,
                    id_type,
                    foreground_ids,
                    background_ids = NULL,
                    threshold = 0.9,
                    transcript = FALSE,
                    stripVersions = TRUE,
                    inflateThresh = 1) {
    if (stripVersions) {
        foreground_ids <- .stripIds(foreground_ids)
        background_ids <- .stripIds(background_ids)
    }
    
    rlang::inform("Mapping foreground ids to ENTREZIDs...")
    fg_id_map <- .mapToEntrezIds( id_type = id_type, orgdb = orgdb,
                                  ids = foreground_ids, threshold = threshold )
    
    # Restrict background pool according to background_ids (if provided)
    if (is.null(background_ids)) {
        rlang::inform("Building background id set...")
        # Otherwise background pool is all records in orgdb
        background_ids <- AnnotationDbi::keys(orgdb, keytype = id_type)
        print_msg <- FALSE
    } else {
        rlang::inform("Mapping background ids to ENTREZIDs...")
        print_msg <- TRUE
    }
    
    bg_id_map <- .mapToEntrezIds( id_type = id_type, orgdb = orgdb,
                                  ids = background_ids, threshold = threshold,
                                  print_msg = print_msg )
    
    # Ensure background pool is the same universe as foreground by dropping any
    # rows with no available value for best mappedID type
    bg_na_filter <- is.na(bg_id_map[, id_type])
    bg_id_map <- bg_id_map[!bg_na_filter, ]
    
    if (!transcript) {
        rlang::inform("Checking for inflation...")
        .checkForInflation( orgdb = orgdb, ids = fg_id_map,
                            id_type = id_type, inflateThresh = inflateThresh )
    }
    
    # Report mapping statistics: Make output list
    mapped <- c( orgdb = orgdb,
                 list( fg_ids = fg_id_map, bg_ids = bg_id_map,
                       userIDtype = id_type, transcript = transcript ) )
    return(mapped)
}


.validateGeneType <- function(geneType, orgdb) {
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


# Filter Foreground and Background ID Sets by Gene Type
#
# @description
# `poolFilter()` takes the mapped foreground and background ID data frames
# (as produced by `mapIDs()`) and, if requested, filters both sets
# to only include genes (or their transcripts) of a specified biotype
# (e.g. “protein-coding”).
#
# @param mapped    A list returned by `mapIDs()`, containing at least:
#   \itemize{
#     \item `fg_ids`: data.frame with columns `entrez` and `mappedID`
#     (foreground).
#     \item `bg_ids`: data.frame with columns `entrez` and `mappedID`
#     (background).
#     \item `so_obj`: a `src_organism` object for transcript lookups.
#     \item `orgdb`:  the loaded OrgDb package object.
#     \item `transcript`: logical, whether IDs are transcripts.
#   }
# @param geneType  Optional character scalar; if not NULL, only genes of this
#   biotype (`GENETYPE` in the OrgDb) will be retained.  Valid values vary by
#   organism (e.g. “protein-coding”, “lncRNA”, etc.).
#
# @return
# The original `mapped` list, but with `fg_ids` and `bg_ids` replaced by
# filtered versions (only rows matching `geneType`, if provided).
#
# @examples
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
#
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# orgdb <- org.Hs.eg.db
#
# my_genes <- c("ENSG00000139618", "ENSG00000157764")
# ids <- mapIDs(
#     orgdb = orgdb,
#     foreground_ids = my_genes,
#     id_type = "ENSEMBL"
# )
# filtered <- poolFilter(ids, geneType = "protein-coding")
#
# @keywords internal
# @export
.poolFilter <- function(mapped, geneType = NULL) {
    # unpack
    fg_df <- mapped$fg_ids
    bg_df <- mapped$bg_ids
    transcript <- mapped$transcript
    orgdb <- mapped$orgdb
    
    # geneType filtering, only if requested
    if (!is.null(geneType)) {
        .validateGeneType(geneType, orgdb)
        
        if (transcript) {
            # transcripts that match geneType
            ensembl_trans <- AnnotationDbi::keys(orgdb, "ENSEMBLTRANS")
            ensembl_trans_df <- AnnotationDbi::select(orgdb,
                                                      keys = ensembl_trans,
                                                      columns = c(
                                                          "ENTREZID",
                                                          "GENETYPE"
                                                      ),
                                                      keytype = "ENSEMBLTRANS"
            )
            keep_trans <- ensembl_trans_df[ensembl_trans_df$GENETYPE ==
                                               geneType, ]
            
            # restrict both bg and fg to transcripts whose gene is in geneType
            bg_df <- bg_df[bg_df$ENTREZID %in% keep_trans$ENTREZID, ]
            fg_df <- fg_df[fg_df$ENTREZID %in% keep_trans$ENTREZID, ]
        } else {
            # gene‐mode: intersect by entrez
            bg_df <- bg_df[bg_df$GENETYPE == geneType, ]
            fg_df <- fg_df[fg_df$GENETYPE == geneType, ]
        }
    }
    
    # reinject filtered dfs and return all components
    mapped$fg_ids <- fg_df
    mapped$bg_ids <- bg_df
    
    return(mapped)
}
