# convert SEM to Motif
.formatMotifs <- function(sem, label) {
    .SD <- NULL
    ppms <- convertSEMsToPPMs(getSEMs(sem))
    motifs <- lapply(
        seq_along(ppms),
        function(i) {
            universalmotif::create_motif(ppms[[i]],
                name = names(ppms[i])
            )
        }
    )

    for (i in seq_along(motifs)) {
        motifs[[i]]["altname"] <- semData(sem)[motifs[[i]]["name"],
            .SD,
            .SDcols = label
        ] |>
            unlist() |>
            unname()
    }
    return(motifs)
}


# construct a matrix comparing similarity of each motif
.constructComparisons <- function(motifs, labels, method) {
    comparisons <- universalmotif::compare_motifs(motifs,
        method = method,
        min.mean.ic = 0
    )
    colnames(comparisons) <- labels
    rownames(comparisons) <- labels

    comparisons <- 1 - comparisons
    comparisons <- stats::as.dist(comparisons)

    comparisons <- stats::hclust(comparisons)
    return(comparisons)
}


.addTipLabels <- function(circ, sigIds, textCex, textCols) {
    group <- NA
    if (length(sigIds) > 0) {
        withCallingHandlers(
            {circ <- ggtree::groupOTU(circ, sigIds)},
            message = function(w) if (grepl("Invaild edge matrix", 
                                            conditionMessage(w))) 
                invokeRestart("muffleMessage")
        )
        circ <- circ + ggtree::geom_tiplab(ggplot2::aes(color = group), 
                                           align = TRUE, 
                                           size = textCex, 
                                           offset = 0.065, 
                                           linesize = 0) +
            ggplot2::scale_color_manual(values=c(textCols[1], textCols[2]), 
                                        guide = "none")
    } else {
        circ <- circ + ggtree::geom_tiplab(color = textCols[1], 
                                           align = TRUE, 
                                           size = textCex, 
                                           offset = 0.1, 
                                           linesize = 0)
    }
    return(circ)
}


#' Plot the results of `enrichSEMs`
#'
#' Generates a circular dendrogram, clustering SNP Effect Matrices on
#'  similarity and a heatmap representing the -log10 transformed
#'  adjusted p-value of a SEMPLR enrichment.
#'
#' @param e The resulting data.table from `enrichSEMs`
#' @param sem A `SNPEffectMatrixCollection` object
#' @param label Column in semData(sem) to use for tree labels
#' @param method Method to use for SEM comparison.
#' See ?universalmotif::compare_motifs for options.
#' @param threshold The adjusted p-value threshold for coloring SEMs
#' @param lineWidth A numeric specifying the dendrogram line width
#' @param textCols A vector of two colors to label non-significant and
#' significant SEMs respectively.
#' @param textCex Text size of SEM labels.
#' @param heatmapCols A vector of two colors to use for the heatmap, ordered
#' low to high `-log10(padj)`.
#' @param pvalRange A vector of 2 numerics to use as the scale range for the
#' heatmap of `-log10(padj)`.
#'
#' @return a `ggtree` object
#'
#' @examples
#' # load SEMs
#'
#' # note that this is a small example for demonstration purposes
#' # in actual enrichment analyses sets of 100+ ranges are recommended
#'
#' # create a GRanges object
#' gr <- GenomicRanges::GRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009
#' )
#'
#' # calculate binding propensity
#' sb <- scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
#' e <- enrichSEMs(sb, SEMC)
#' plotEnrich(e, SEMC)
#'
#' @return NULL
#'
#' @export
plotEnrich <- function(e, sem,
                       label = "transcription_factor",
                       method = "WPCC",
                       threshold = 0.05,
                       lineWidth = 0.5,
                       textCols = c("darkgrey", "black"),
                       textCex = 1,
                       heatmapCols = c("white", "red"),
                       pvalRange = c(0, 20)) {
    .SD <- group <- NULL
    sk <- semData(sem) |> data.table::key()
    em <- merge(semData(sem), e, by.x = sk, by.y = "SEM")
    motifs <- .formatMotifs(sem, label)
    labels <- lapply(motifs, function(x) x["altname"]) |> unlist()

    comparisons <- .constructComparisons(
        motifs = motifs,
        labels = labels,
        method = method )
    
    den <- stats::as.dendrogram(comparisons)
    circ <- ggtree::ggtree(den, layout = "circular") + 
        ggtree::theme_dendrogram(bgcolor = "transparent", 
                                 fgcolor = "transparent")
    em_df <- as.data.frame(-log10(em[, "padj"]))
    rownames(em_df) <- em[, .SD, .SDcols = label] |> unlist()
    colnames(em_df) <- "padj"
    sigIds <- em[, .SD, .SDcols = label][which(em$padj <= threshold)] |>
        unlist() |> unname()
    
    circ <- .addTipLabels(circ = circ, sigIds = sigIds, 
                          textCex = textCex, textCols = textCols)
    
    withCallingHandlers( {
        plt <- ggtree::gheatmap(circ, em_df, width=.08, colnames_angle=0, 
                                offset = -0.02, colnames = FALSE) +
            ggplot2::scale_fill_gradient(name = "-log10(Adj. P-value)",
                                         low = heatmapCols[1], 
                                         high = heatmapCols[2], 
                                         limits = pvalRange, 
                                         oob = scales::squish)
            },
        message = function(w) {
            if (grepl("Invaild edge matrix", conditionMessage(w)) |
                grepl("Missing column: label.", conditionMessage(w))) 
                invokeRestart("muffleMessage") }
    )
    return(plt)
}
