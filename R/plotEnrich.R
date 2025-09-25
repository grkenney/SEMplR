# convert SEM to Motif
.formatMotifs <- \(sem, label) {
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
.constructComparisons <- \(motifs, labels, method) {
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


.circlizePlot <- \(em, sem, comps, ds, col_fun, label, sigIds, sigCols) {
    .SD <- NULL
    dend <- stats::as.dendrogram(comps)

    padjs <- data.frame(padj = em$padj)
    rownames(padjs) <- em[, .SD, .SDcols = label] |> unlist()
    column_od <- comps$order
    ordered_mtx <- as.matrix(padjs[column_od, 1])
    rownames(ordered_mtx) <- rownames(padjs)[comps$order]

    ordered_text_cols <- lapply(
        rownames(ordered_mtx),
        function(x) {
            ifelse(x %in% sigIds,
                sigCols[2], sigCols[1]
            )
        }
    ) |>
        unlist()

    circlize::circos.clear()

    circlize::circos.heatmap(ordered_mtx,
        col = col_fun,
        rownames.side = "outside",
        cluster = FALSE, clustering.method = NULL,
        track.height = 0.04,
        dend.callback = function(dend, m, si) {
            dendsort::dendsort(dend)
        }, rownames.col = ordered_text_cols
    )

    graphics::par(new = TRUE)
    circlize::circos.trackPlotRegion(
        ylim = c(0, 1),
        bg.border = NA,
        track.height = min(ds) * 0.11 / 1.02,
        panel.fun = function(x, y) {
            circlize::circos.dendrogram(
                dend = dend,
                facing = "outside",
                max_height = 1
            )
        }
    )
    circlize::circos.clear()
    return()
}



.addLegend <- \(em, sem, comps, heatmapCols, label, sigIds, sigCols, textCex) {
    graphics::plot.new()
    circle_size <- grid::unit(1, "snpc") # snpc unit gives you a square region

    grid::pushViewport(grid::viewport(
        x = 0, y = 0.5, width = circle_size,
        height = circle_size,
        just = c("left", "center")
    ))
    graphics::par(
        omi = gridBase::gridOMI(), new = TRUE,
        cex = textCex, mar = c(0, 0, 0, 0)
    )
    ds <- grDevices::dev.size()

    col_fun <- circlize::colorRamp2(c(-2, 2), heatmapCols)
    .circlizePlot(
        em = em, sem = sem, comps = comps,
        ds = ds, col_fun = col_fun,
        label = label,
        sigIds = sigIds, sigCols = sigCols
    )
    grid::upViewport()

    lgd <- ComplexHeatmap::Legend(
        title = "Adj. P-value", col_fun = col_fun,
        title_gp = grid::gpar(fontsize = 8),
        labels_gp = grid::gpar(fontsize = 8),
        direction = "horizontal"
    )

    lgd_list <- ComplexHeatmap::packLegend(lgd,
        max_height = unit(0.9 * ds[2], "inch")
    )
    ComplexHeatmap::draw(lgd_list,
        x = circle_size * 0.95,
        y = circle_size * 0.95,
        just = c("center", "top")
    )
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
#' @param textCols A vector of two colors to label non-significant and
#' significant SEMs respectively.
#' @param textCex Text size of SEM labels.
#' @param heatmapCols A vector of two colors to use for the heatmap, ordered
#' low to high -log10(padj).
#'
#' @return a `ggtree` object
#'
#' @importFrom circlize circos.trackPlotRegion
#'
#' @examples
#' # load SEMs
#' data(SEMC)
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
plotEnrich <- \(e, sem,
    label = "transcription_factor",
    method = "WPCC",
    threshold = 0.05,
    textCols = c("darkgrey", "black"),
    textCex = 0.7,
    heatmapCols = c("white", "red")) {
    .SD <- SEM_KEY <- NULL

    em <- merge(semData(sem), e, by.x = "SEM_KEY", by.y = "SEM")

    motifs <- .formatMotifs(sem, label)
    labels <- lapply(motifs, function(x) x["altname"]) |> unlist()

    comparisons <- .constructComparisons(
        motifs = motifs,
        labels = labels,
        method = method
    )

    sigIds <- em[, .SD, .SDcols = label][which(em$padj <= threshold)] |>
        unlist() |>
        unname()

    .addLegend(
        em = em, sem = sem,
        comps = comparisons, heatmapCols = heatmapCols,
        label = label,
        sigIds = sigIds, sigCols = textCols, textCex = textCex
    )
}
