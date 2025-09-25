# make sure motif is not null and that provided motif id is in the collection
.validateMotifInCollection <- \(sem, motif) {
    # motif required if providing a SNPEffectMatrixCollection
    if (is.null(motif)) {
        rlang::abort(paste0(
            "must specify the motif to plot if supplying a ",
            "SNPEffectMatrixCollection"
        ))
    }

    # check that motif given is in collection
    if (!(motif %in% names(getSEMs(sem)))) {
        rlang::abort(paste0(
            "provided motif not found in SNPEffectMatrixCollection.",
            "See the SEM_KEY column in semData(sem) for available motif ids."
        ))
    }
}


# validate that sem is a SNPEffectMatrix object or a SNPEffectMatrixCollection
# and that the motif id provided is valid
# define the baseline and matrix for the SEM
.definePlotSEMParams <- \(sem, motif) {
    if (is(sem, "SNPEffectMatrix")) {
        if (!is.null(motif)) {
            message(
                "motif ignored because SNPEffectMatrix object given.",
                " See motif argument description in ?plotSEM."
            )
        }
        sem_baseline <- getBaseline(sem)
        sem_mtx <- getSEM(sem)
    } else if (is(sem, "SNPEffectMatrixCollection")) {
        .validateMotifInCollection(sem, motif)
        sem_baseline <- getBaseline(getSEMs(sem, semId = motif))
        sem_mtx <- getSEM(getSEMs(sem, semId = motif))
    } else {
        rlang::abort(paste0(
            "sem must be a SNPEffectMatrix or a ",
            "SNPEffectMatrixCollection."
        ))
    }

    return(list(
        sem_baseline = sem_baseline,
        sem_mtx = sem_mtx
    ))
}


# pivot SEM matrix into a long data.table with score, position, and bp columns
.formatPlotSEMTable <- \(sem_mtx) {
    bp <- rep(c("A", "C", "G", "T"), each = nrow(sem_mtx))
    motif_pos <- rep(seq_len(nrow(sem_mtx)), 4)
    sem_score <- sem_mtx |>
        unlist()

    sem_mtx_long <- data.table(sem_score, motif_pos, bp)
    return(sem_mtx_long)
}


# define colors to plot nucleotides
# adds hseq and text_color columns to sem_mtx_long for the sequence to highlight
# and the color to plot each nucleotide
.defineNucleotideColors <- \(sem_mtx_long, motifSeq) {
    motif_length <- max(sem_mtx_long$motif_pos)
    # if sequence is provided for plotting,
    if (!is.null(motifSeq)) {
        if (nchar(motifSeq) != motif_length) {
            rlang::abort(paste0(
                paste0(
                    "Number of characters in motifSeq (", nchar(motifSeq),
                    ") does not equal the number of rows in the SEM (",
                    motif_length, ")"
                )
            ))
        }
        mseq <- unlist(strsplit(motifSeq, ""))

        # going to plot as white first and then overlay bolded, colored sequence
        mseq_dt <- data.table(mseq, seq_len(motif_length), mseq, "white") |>
            stats::setNames(c("bp", "motif_pos", "mseq", "text_color"))

        sem_mtx_long <- merge(sem_mtx_long, mseq_dt,
            by = c("bp", "motif_pos"),
            all.x = TRUE
        )

        sem_mtx_long$text_color[is.na(sem_mtx_long$text_color)] <- "#d7dbdd"
        sem_mtx_long$mseq[is.na(sem_mtx_long$mseq)] <- ""
    } else {
        sem_mtx_long$text_color <- "black"
    }
    return(sem_mtx_long)
}


.createBasePlotSEM <- \(sem_mtx_long, sem_mtx, highlight, motif,
    hwidth, hcol, sem_baseline, size) {
    motif_pos <- sem_score <- NULL

    motif_plot <- ggplot2::ggplot(
        sem_mtx_long,
        aes(x = motif_pos, y = sem_score)
    ) +
        theme_classic()

    if (!is.null(highlight)) {
        motif_plot <- motif_plot + geom_vline(
            xintercept = highlight,
            linewidth = hwidth,
            col = hcol, alpha = 0.1
        )
    }

    motif_plot <- motif_plot +
        geom_hline(yintercept = 0, color = "#eaeaea", linewidth = 1) +
        geom_hline(
            yintercept = sem_baseline, linetype = "dashed",
            color = "#eaeaea", linewidth = 1
        ) +
        geom_text(aes(label = sem_mtx_long$bp),
            size = size,
            col = sem_mtx_long$text_color
        ) +
        scale_x_continuous(breaks = scales::breaks_pretty(n = nrow(sem_mtx))) +
        theme(legend.position = "none") +
        ggplot2::xlab("Position in Motif") +
        ggplot2::ylab("SNP Effect Score") +
        theme(text = element_text(size = 12)) +
        ggtitle(motif)
    return(motif_plot)
}


#' Plot SEM scores for each nucleic acid in each position of the motif
#'
#' @param sem A SNPEffectMatrix or SNPEffectMatrixCollection object.
#' @param motif sem id to plot. If providing a SNPEffectMatrixCollection object,
#' must match a name in the list. If providing a single SNPEffectMatrix object,
#' this parameter is ignored and the SNPEffectMatrix's semId is used for
#' plotting.
#' @param motifSeq Character sequence to color on plot
#' @param highlight Index of variant location on reference
#' @param hcol Color of motifSeq bases
#' @param hwidth Width of highlight bar
#' @param size Font size of nucleotides in plot
#'
#' @import ggplot2
#'
#' @return a `ggplot` with sem scores for each nucleic acid per position
#'
#' @export
#'
#' @examples
#' library(VariantAnnotation)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#'
#' # create a list of SNPEffectMatrix objects
#' s <- SNPEffectMatrix(sem, baseline = -1, semId = "sem_id")
#'
#' # plot the motif
#' plotSEM(s, getSEMId(s))
#'
#' # color by sequence
#' plotSEM(s, getSEMId(s), motifSeq = "ACT", highlight = 2)
plotSEM <- function(sem, motif = NULL,
                    motifSeq = NULL, highlight = NULL,
                    hcol = "dodgerblue", hwidth = 15,
                    size = 7) {
    alt <- ref <- motif_pos <- mseq <- .SD <-
        semId <- sem_score <- sm <- varId <- NA

    sem_params <- .definePlotSEMParams(sem, motif)
    sem_baseline <- sem_params$sem_baseline
    sem_mtx <- sem_params$sem_mtx

    sem_mtx_long <- .formatPlotSEMTable(sem_mtx)

    sem_mtx_long <- .defineNucleotideColors(
        sem_mtx_long = sem_mtx_long,
        motifSeq = motifSeq
    )

    motif_plot <- .createBasePlotSEM(
        sem_mtx_long, sem_mtx,
        highlight, motif, hwidth, hcol,
        sem_baseline, size
    )


    if (!is.null(motifSeq)) {
        motif_plot <- motif_plot +
            geom_text(aes(label = mseq),
                size = size,
                col = hcol, fontface = "bold"
            )
    }

    return(motif_plot)
}
