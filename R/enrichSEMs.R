# given a list of DNA sequences, scramble the entire sequence
.scrambleSeqs <- function(seqs) {
    scrambled <- lapply(
        seqs,
        function(x) stringi::stri_rand_shuffle(x)
    ) |> unlist()
    return(scrambled)
}


# perform a binomial test on sem scores versus a scores from a background for a
# single sem
.binomSEMEnrich <- \(xs, bg, semName) {
    xs_sem <- xs[xs$SEM == semName, ]
    bg_sem <- bg[bg$SEM == semName, ]

    successes <- sum(xs_sem$scoreNorm > 0)
    sample_size <- nrow(xs_sem)
    prob_success <- (sum(bg_sem$scoreNorm > 0) + 1) / nrow(bg_sem)

    if (prob_success > 1) {
        prob_success <- 1
    }

    binom_test_res <- stats::binom.test(
        x = successes,
        n = sample_size,
        p = prob_success,
        alternative = "greater"
    )

    n_bound_bg <- sum(bg$scoreNorm[bg$SEM == semName] > 0)
    result <- data.table::data.table(
        SEM = semName,
        pvalue = binom_test_res$p.value,
        n_bound = successes,
        n_bound_bg = n_bound_bg
    )
    return(result)
}


# define the background set if not provided
.defineBackground <- \(x, sem, background, seqs, nFlank, genome) {
    if (is.null(background)) {
        if (is(x, "SEMplScores")) {
            seqs <- getRanges(x)$sequence
        }

        scramb <- .scrambleSeqs(seqs)
        bg <- lapply(
            getSEMs(sem),
            function(y) {
                scoreSequence(
                    sem = as.matrix(getSEM(y)),
                    dna_sequences = scramb,
                    nFlank = nFlank,
                    bl = getBaseline(y),
                    seqIds = seq_along(seqs)
                )
            }
        ) |> data.table::rbindlist(idcol = "SEM")
    } else {
        bg <- scoreBinding(x = background, sem = sem, genome = genome) |>
            scores()
    }
    return(bg)
}


#' Calculates binding enrichment for motif(s)
#'
#' Perform a binomial test to determine if SNP Effect Matrices are bound more
#' often than expected.
#'
#' @param x The scoring table produced by `scoreBinding`
#' @param sem A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object
#' @param background A list of DNA sequences to use for background. The length
#' of each sequence must match the length of sequences in `x`. By default, will
#' scramble the provided sequences.
#' @param seqs The sequences scored in `scoreBinding`
#' @param nFlank Number of flanking nucleotides added to the sequences. Defaults
#' to the length of the longest motif. If no flanks were added (ie, sequences
#' were scored rather than a `GRanges`), use `nFlank = 0`.
#' @param genome A `BSGenome` object. Requried if background isn't specified.
#'
#' @return a `list` of matrices
#'
#' @examples
#'
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
#' enrichSEMs(sb, SEMC)
#'
#' @export
enrichSEMs <- \(x, sem,
    background = NULL, seqs = NULL, nFlank = 0,
    genome = NULL) {
    sem_names <- getSEMs(sem) |> names()

    bg <- .defineBackground(
        x = x,
        sem = sem,
        background = background,
        seqs = seqs,
        nFlank = nFlank,
        genome = genome
    )

    x_scores <- scores(x)
    result <- lapply(
        seq_along(sem_names),
        \(i) .binomSEMEnrich(
            xs = x_scores,
            bg = bg,
            semName = sem_names[i]
        )
    ) |>
        data.table::rbindlist()

    result$padj <- stats::p.adjust(result$pvalue, method = "BH")
    return(result)
}
