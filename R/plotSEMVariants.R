# validate label, variant, and cols parameters
.validatePlotSemVariantsInputs <- \(s, label, semId, cols) {
    # check that sem label is in the sem meta data
    meta_cols <- S4Vectors::mcols(getRanges(s)) |>
        colnames()
    if (!(label %in% meta_cols)) {
        rlang::abort("label not found in S4Vectors::mcols(getRanges(s)).")
    }

    # check that semId is a valid id in s
    if (!(semId %in% scores(s)[, semId])) {
        rlang::abort(paste0(
            "variant not found in SEMplScores object. ",
            semId, " is not in scores(s)[, semId]"
        ))
    }

    # check that cols is length 2
    if (length(cols) != 2) {
        rlang::abort("cols must be a vector of length 2")
    }
}


.createBasePlotSEMVariants <- \(scores_dt, cols, label) {
    refNorm <- altNorm <- NULL
    plt <- ggplot2::ggplot(
        data = scores_dt,
        aes(x = refNorm, y = altNorm)
    ) +
        geom_vline(xintercept = 0, linetype = 2, col = "grey") +
        geom_hline(yintercept = 0, linetype = 2, col = "grey") +
        geom_point(
            data = subset(scores_dt, altNorm < 0 & refNorm < 0),
            size = 1,
            color = "grey"
        ) +
        geom_point(
            data = subset(scores_dt, altNorm > 0 & refNorm < 0),
            size = 1,
            color = cols[1]
        ) +
        geom_point(
            data = subset(scores_dt, altNorm < 0 & refNorm > 0),
            size = 1,
            color = cols[2]
        ) +
        geom_point(
            data = subset(scores_dt, altNorm > 0 & refNorm > 0),
            size = 1,
            color = "grey"
        ) +
        ggrepel::geom_text_repel(
            data = subset(
                scores_dt,
                altNorm > 0 & refNorm < 0
            ),
            mapping = aes(label = .data[[label]]),
            size = 4,
            color = cols[1]
        ) +
        ggrepel::geom_text_repel(
            data = subset(
                scores_dt,
                altNorm < 0 & refNorm > 0
            ),
            mapping = aes(label = .data[[label]]),
            size = 4,
            color = cols[2]
        )
    return(plt)
}


#' Plot non-alt versus alt binding propensity for a single motif
#'
#' @param s a SEMplScores object with scores populated
#' @param sem a single character vector matching a semId in the semplObj
#' @param label column in scores slot of semplObj to use for point labels
#' @param cols vector of length 2 with colors to use for plotting gained
#' and lost motifs respectively
#'
#' @import ggplot2
#'
#' @return a `ggplot` scatter plot of non-alt binding propensity versus alt
#' binding propensity
#'
#' @export
#'
#' @examples
#' library(VariantAnnotation)
#' data(SEMC)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#'
#' # create a VRanges object
#' vr <- VRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009,
#'     ref = "G", alt = "C"
#' )
#'
#' # calculate binding propensity
#' s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
#' plotSEMVariants(s, "IKZF1_HUMAN.GM12878")
#'
plotSEMVariants <- function(s, sem, label = "varId",
                            cols = c("#F8766D", "dodgerblue2")) {
    refNorm <- altNorm <- ix <- semId <- NA
    .validatePlotSemVariantsInputs <- \(s = s, label = label,
        semId = sem, cols = cols)

    rlang::inform(paste0("Plotting ", sem, "..."))
    ix <- sem == scores(s)[, semId]
    scores_dt <- scores(s)[ix, ]

    var_plot <- .createBasePlotSEMVariants(scores_dt, cols, label)
    var_plot <- var_plot +
        scale_x_continuous(
            breaks = scales::pretty_breaks(),
            limits = \(x) ifelse(abs(x) < 1, c(-1, 1), x)
        ) +
        scale_y_continuous(
            breaks = scales::pretty_breaks(),
            limits = \(x) ifelse(abs(x) < 1, c(-1, 1), x)
        ) +
        labs(
            x = "ref binding propensity",
            y = "alt binding propensity"
        ) +
        theme_classic() +
        theme(
            panel.grid = element_blank(),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15)
        )

    return(var_plot)
}
