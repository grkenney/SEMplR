# convert SEM to Motif
.formatMotifs <- \(sem, label) {
    .SD <- NULL
    ppms <- convertSEMsToPPMs(sems(sem))
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
.constructComparisons <- \(motifs, labels, method) {
    comparisons <- universalmotif::compare_motifs(motifs,
        method = method,
        min.mean.ic = 0
    )
    colnames(comparisons) <- labels
    rownames(comparisons) <- labels

    comparisons <- 1 - comparisons
    comparisons <- stats::as.dist(comparisons)

    comparisons <- ape::as.phylo(stats::hclust(comparisons))
    return(comparisons)
}


# construct phylogenetic tree of motif similarity
.plotMotifTree <- \(e, sem, threshold, comparisons, sigCol) {
    sm <- NA
    # identify significant motifs
    sig_motifs <- e$SEM[e$padj <= threshold]

    # store TF names for significant motifs for labeling
    sm_tf <- unlist(semData(sem)[sig_motifs, "transcription_factor"])
    comparisons <- ggtree::groupOTU(comparisons,
        sm_tf,
        group_name = "sm"
    )

    tree <- ggtree::ggtree(comparisons, layout = "fan") +
        theme(legend.position = "", legend.title = element_blank()) +
        theme(plot.margin = margin(10, 10, 10, 100)) +
        ggtree::geom_tiplab(aes(colour = sm),
            as_ylab = FALSE,
            size = 1.5,
            offset = .05
        ) +
        scale_color_manual(values = c("grey", sigCol), guide = "none")
    return(tree)
}


# add heatmap to ggtree object
.plotTreeHeatmap <- \(plt, e, labels) {
    heatmapData <- e[, "padj"] |> as.data.frame()
    colnames(heatmapData) <- "-log10(padj)"
    rownames(heatmapData) <- labels

    plt <- ggtree::gheatmap(plt, -log10(heatmapData),
        offset = -0.014,
        width = 0.12, colnames = FALSE,
        font.size = 2.5
    )
    # add the scale bar
    plt <- plt +
        theme(
            legend.position = "bottom",
            legend.title = element_text(color = "black"),
            legend.title.position = "top"
        ) +
        scale_fill_viridis_c(
            na.value = "white", direction = -1,
            name = "-log10(padj)"
        )
    return(plt)
}


#' Plot the results of `enrichSEMs`
#'
#' Generates a cladogram clustered by SNP Effect Matrix similarity and
#' a heatmap representing the adjusted p-value of the enrichment.
#'
#' @param e The resulting data.table from `enrichSEMs`
#' @param sem A `SNPEffectMatrixCollection` object
#' @param label Column in semData(sem) to use for tree labels
#' @param method Method to use for SEM comparison.
#' See ?universalmotif::compare_motifs for options.
#' @param threshold The adjusted p-value threshold for coloring SEMs
#' @param sigCol Color to label significant SEM labels
#'
#' @return a `ggtree` object
#'
#' @examples
#' # load SEMs
#' data(sc)
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
#' sb <- scoreBinding(gr, sc, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' e <- enrichSEMs(sb, sc) 
#' # plotEnrich(e, sc)
#'
#' @export
plotEnrich <- \(e, sem,
    label = "transcription_factor",
    method = "WPCC",
    threshold = 0.05,
    sigCol = "purple4") {
    motifs <- .formatMotifs(sem, label)
    labels <- lapply(motifs, function(x) x["altname"]) |> unlist()

    comparisons <- .constructComparisons(motifs, labels, method = method)

    plt <- .plotMotifTree(
        e = e, sem = sem, threshold = threshold,
        comparisons = comparisons, sigCol = sigCol
    )

    plt <- .plotTreeHeatmap(plt = plt, e = e, labels = labels)

    return(plt)
}
