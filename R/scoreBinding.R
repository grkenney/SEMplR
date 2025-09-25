# given a GRanges object, x, build a position identifier string of format
# seqname:position
.makePositionId <- \(x) {
    start_pos <- IRanges::start(IRanges::ranges(x))
    end_pos <- IRanges::end(IRanges::ranges(x))
    sn <- GenomeInfoDb::seqnames(x)
    pos_str <- ifelse(start_pos == end_pos,
        start_pos, paste0(start_pos, "-", end_pos)
    )

    id <- paste0(sn, ":", pos_str)
    return(id)
}


# test if x is a list or vector of characters that resembles DNA
.testIfSequenceList <- \(x) {
    # if is neither a vector or a list, it's not collection of DNA seqs
    if (!(is.vector(x)) & !(is.list(x))) {
        # if it's not a GRanges, error
        if (!is(x, "GRanges")) {
            rlang::abort(paste0(
                is(x)[1], " is not an accepted class for scoreBinding. x ",
                "must be a GRanges or VRanges object or vector of DNA sequences"
            ))
        } else {
            return(FALSE)
        }
    }

    # check that x is a character list
    x_unlist <- x |> unlist()
    if (!(is.character(x_unlist))) {
        rlang::abort(paste0(
            is(x)[1], " is not an accepted class for scoreBinding.
    x must be a GRanges or VRanges object or vector of DNA sequences"
        ))
    }
    # check if all elements only contain characters A, C, G, T
    unique_items <- x_unlist |>
        strsplit(split = "") |>
        unlist() |>
        unique()
    if (!all(unique_items %in% c("A", "C", "G", "T"))) {
        invalid_char <- unique_items[!(unique_items %in% c("A", "C", "G", "T"))]
        rlang::abort(paste0(
            "x contains invalid nucleotides: ",
            paste(invalid_char, collapse = ", "),
            ". All characters in x must ",
            "be 'A', 'C', 'G', or 'T'"
        ))
    }
    return(TRUE)
}


# add sequence and unique id columns to GRanges meta data
# x: GRanges object
# sem: SNPEffectMatrixCollection object
# genome: BSGenomes object
# nFlank: integer for number of flanking nucleotides to pull
# seqId: column name in x to use as unique id.
# allele: allele column name
.prepRangeMetadata <- \(x, sem, genome, nFlank, seqId) {
    # generate a unique sequence id if one is not provided
    if (is.null(seqId)) {
        id <- lapply(
            seq_along(x),
            function(i) .makePositionId(x = x[i])
        ) |>
            unlist()
        S4Vectors::mcols(x)[, "id"] <- id
    } else {
        # check that the seqId column exists in the meta data
        if (!(seqId %in% colnames(S4Vectors::mcols(x)))) {
            rlang::abort(paste0(
                "'", seqId,
                "' is not a meta data column. See mcols(x) for ",
                "the meta data column names or don't specify seqId ",
                "to generate unique ids."
            ))
        }

        # make sure ids are unique
        id <- S4Vectors::mcols(x)[, seqId]
        if (length(unique(id)) != length(id)) {
            rlang::abort(paste0(
                "entries in the ", seqId,
                " column are not unique."
            ))
        }
    }

    # must provide a reference genome if sequences are not given
    if (is.null(genome)) {
        rlang::abort("must specify a genome if providing a GRanges object")
    }

    # pull the sequence for each range
    x <- getRangeSeqs(x,
        genome = genome,
        up = nFlank, down = nFlank
    )

    return(x)
}


#' Calculate binding propensity for all SEM motifs and
#' genomic positions provided
#'
#' @param x `GRanges` object or a vector of DNA sequences
#' @param sem A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object
#' @param genome A `BSgenome` object for the genome build to use. ie.
#' `BSgenome.Hsapiens.UCSC.hg19::Hsapiens`. Required if providing a GRanges
#' object. Ignored if providing a vector of sequences.
#' @param nFlank Number of flanking nucleotides to add to provided range. By
#' default will add flank equal to the length of the longest motif. Ignored
#' if providing a vector of sequences.
#' @param seqId Column in `GRanges` object to use for unique id.
#' By default, ids will be generated from the `seqnames` and `ranges.`
#' Ignored if not providing a `GRanges` object.
#'
#' @return If a `GRanges` object is provided, return a `SEMplScores` object.
#' If a list of sequences is provided, just return the scoring table
#'
#' @export
#'
#' @examples
#' # load SEMs
#' data(SEMC)
#'
#' # create a GRanges object
#' gr <- GenomicRanges::GRanges(
#'     seqnames = "chr12",
#'     ranges = 94136009
#' )
#'
#' # calculate binding propensity
#' scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#'
scoreBinding <- \(x, sem, genome, nFlank = NULL,
    seqId = NULL) {
    # make sure nFlank is an integer, if provided
    if (!is.null(nFlank) & !is.numeric(nFlank)) {
        rlang::abort("nFlank must be an integer.")
    } else if (!is.null(nFlank)) {
        nFlank <- as.integer(nFlank)
    }

    # determine x input object class/type
    is_sequence_list <- .testIfSequenceList(x)
    # convert sem to a collection if it isn't one already
    sem <- .convertToSNPEffectMatrixCollection(x = sem)

    # if sequence list not provided, grab sequences from the reference genome
    if (!is_sequence_list) {
        # if the nFlank is not given, use length of the longest SEM
        if (is.null(nFlank)) {
            nFlank <- lapply(getSEMs(sem), function(x) {
                nrow(getSEM(x))
            }) |>
                unlist() |>
                max()
        }

        x <- .prepRangeMetadata(
            x = x,
            sem = sem,
            genome = genome,
            nFlank = nFlank,
            seqId = seqId
        )
        id <- x$id
        seqs <- x$sequence
    } else {
        # if sequence list is given, use sequences and apply an integer id
        seqs <- x
        id <- seq_along(seqs)
        if (is.null(nFlank)) {
            nFlank <- 0
        }
    }

    s <- lapply(
        getSEMs(sem),
        function(y) {
            scoreSequence(
                sem = as.matrix(getSEM(y)),
                dna_sequences = seqs,
                nFlank = nFlank,
                bl = getBaseline(y),
                seqIds = id
            )
        }
    )

    # combine nested list of tables into a single data.table with SEM column
    s <- s |> data.table::rbindlist(idcol = "SEM")
    s <- s[, c("seqId", "SEM", "score", "scoreNorm", "index", "seq")]

    # if GRanges, return in SEMplScores object, otherwise, return the data.table
    if (!is_sequence_list) {
        ss <- SEMplScores(ranges = x, semData = semData(sem), scores = s)
    } else {
        ss <- s
    }

    return(ss)
}
