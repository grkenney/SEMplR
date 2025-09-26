b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens


test_that(".makePositionId builds single position id", {
    x <- GenomicRanges::GRanges(
        seqnames = "chr12",
        ranges = 94136009
    )
    id_a <- .makePositionId(x)
    expect_equal(id_a, "chr12:94136009")
})


test_that(".makePositionId builds single range id", {
    x <- GenomicRanges::GRanges(
        seqnames = "chr12",
        ranges = IRanges::IRanges(94136009, 94136010)
    )
    id_a <- .makePositionId(x)
    expect_equal(id_a, "chr12:94136009-94136010")
})


test_that(".makePositionId builds multiple ids", {
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12", "chr13"),
        ranges = IRanges::IRanges(
            c(94136009, 50000000),
            c(94136010, 50000000)
        )
    )
    id_a <- .makePositionId(x)
    expect_equal(id_a, c("chr12:94136009-94136010", "chr13:50000000"))
})


test_that(".testIfSequenceList test invalid x inputs", {
    expect_false(.testIfSequenceList(GenomicRanges::GRanges()))

    expect_error(.testIfSequenceList(list(TRUE, FALSE)),
        regexp = "is not an accepted class for scoreBinding."
    )
    expect_error(.testIfSequenceList(list(1, 2, 3)),
        regexp = "is not an accepted class for scoreBinding."
    )
    expect_error(.testIfSequenceList(c()),
        regexp = "is not an accepted class for scoreBinding."
    )
    expect_error(.testIfSequenceList(list()),
        regexp = "is not an accepted class for scoreBinding."
    )
    expect_error(.testIfSequenceList(TRUE),
        regexp = "is not an accepted class for scoreBinding."
    )
    expect_error(.testIfSequenceList(list("1", "2")),
        regexp = "contains invalid nucleotide"
    )
    expect_error(.testIfSequenceList(list("AC", "TGN")),
        regexp = "contains invalid nucleotide"
    )
    expect_error(.testIfSequenceList(list("N", "A")),
        regexp = "contains invalid nucleotide"
    )
})


test_that(".testIfSequenceList test valid x inputs", {
    expect_true(.testIfSequenceList(list("AC", "AT")))
    expect_true(.testIfSequenceList("AC"))
    expect_true(.testIfSequenceList(c("AC", "TAG")))
})


test_that(".prepRangeMetadata builds id and grabs sequence", {
    # no flank
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009))
    )
    x_a <- .prepRangeMetadata(
        x = x, sem = sc,
        genome = b, nFlank = 0,
        seqId = NULL
    )
    x_e <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        id = "chr12:94136009",
        sequence = "G"
    )
    expect_equal(x_a, x_e)

    # with flank
    x_a <- .prepRangeMetadata(
        x = x, sem = sc,
        genome = b, nFlank = 2,
        seqId = NULL
    )
    x_e <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        id = "chr12:94136009",
        sequence = "TTGAG"
    )
    expect_equal(x_a, x_e)
})


test_that(".prepRangeMetadata fails on invalid input", {
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12", "chr13"),
        ranges = IRanges::IRanges(c(94136009, 10000000)),
        ids = c("A", "A")
    )
    # invalid id column
    expect_error(
        .prepRangeMetadata(
            x = x, sem = sc,
            genome = b, nFlank = 0,
            seqId = "foo"
        ),
        regexp = "is not a meta data column."
    )
    expect_error(
        .prepRangeMetadata(
            x = x, sem = sc,
            genome = b, nFlank = 0,
            seqId = "ids"
        ),
        regexp = "column are not unique."
    )

    # genome is NULL
    expect_error(
        .prepRangeMetadata(
            x = x, sem = sc,
            genome = NULL, nFlank = 0,
            seqId = NULL
        ),
        regexp = "must specify a genome if providing a GRanges object"
    )
})


test_that("scoreBinding invalid inputs", {
    # invalid x
    expect_error(scoreBinding(x = TRUE, sem = sc, genome = b),
        regexp = "is not an accepted class for scoreBinding"
    )
    expect_error(scoreBinding(x = NULL, sem = sc, genome = b),
        regexp = "is not an accepted class for scoreBinding"
    )
    expect_error(scoreBinding(x = 1, sem = sc, genome = b),
        regexp = "is not an accepted class for scoreBinding"
    )

    # invalid nFlank
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009))
    )
    expect_error(scoreBinding(x = x, sem = sc, genome = b, nFlank = "A"),
        regexp = "nFlank must be an integer."
    )
})


test_that("scoreBinding range scoring", {
    # invalid nFlank
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009))
    )
    sb_a <- scoreBinding(
        x = x,
        sem = getSEMs(SEMC, "JUN"),
        genome = b
    )
    scores_e <- data.table(
        seqId = "chr12:94136009",
        SEM = "JUN",
        score = -0.96313179,
        scoreNorm = -0.02038131,
        index = 7,
        seq = "TGAGGCA"
    )
    x_e <- x
    S4Vectors::mcols(x_e)[, "id"] <- c("chr12:94136009")
    S4Vectors::mcols(x_e)[, "sequence"] <- c("AGGCTTTGAGGCATC")
    sb_e <- SEMplScores(
        ranges = x_e,
        semData = data.table(),
        scores = scores_e
    )
    expect_equal(sb_a, sb_e, tolerance = 1e-6)
})


test_that("scoreBinding sequence list scoring", {
    # invalid nFlank
    x <- "AGGCTTTGAGGCATC"
    sb_a <- scoreBinding(
        x = x,
        sem = getSEMs(SEMC, "JUN"),
        nFlank = 7, seqId = "A"
    )
    sb_e <- data.table(
        seqId = "1",
        SEM = "JUN",
        score = -0.96313179,
        scoreNorm = -0.02038131,
        index = 7,
        seq = "TGAGGCA"
    )
    expect_equal(sb_a, sb_e, tolerance = 1e-6)
})


test_that("scoreBinding test when flank is shorter than SEM", {
    # invalid nFlank
    x <- "TTGAGTCAA"
    sb_a <- scoreBinding(
        x = x,
        sem = getSEMs(SEMC, "JUN"),
        nFlank = 1, seqId = "A"
    )
    sb_e <- data.table(
        seqId = "1",
        SEM = "JUN",
        score = 0.0100808,
        scoreNorm = 0.9231947,
        index = 2,
        seq = "TGAGTCA"
    )
    expect_equal(sb_a, sb_e, tolerance = 1e-6)
})


test_that("scoreBinding test when flank is zero", {
    # invalid nFlank
    x <- "TTGAGTCAA"
    sb_a <- scoreBinding(
        x = x,
        sem = getSEMs(SEMC, "JUN"),
        nFlank = 0, seqId = "A"
    )
    sb_e <- data.table(
        seqId = "1",
        SEM = "JUN",
        score = 0.0100808,
        scoreNorm = 0.9231947,
        index = 2,
        seq = "TGAGTCA"
    )
    expect_equal(sb_a, sb_e, tolerance = 1e-6)
})


test_that("scoreBinding test when flank null", {
    # invalid nFlank
    x <- "TTGAGTCAA"
    sb_a <- scoreBinding(
        x = x,
        sem = getSEMs(SEMC, "JUN"),
        seqId = "A"
    )
    sb_e <- data.table(
        seqId = "1",
        SEM = "JUN",
        score = 0.0100808,
        scoreNorm = 0.9231947,
        index = 2,
        seq = "TGAGTCA"
    )
    expect_equal(sb_a, sb_e, tolerance = 1e-6)
})
