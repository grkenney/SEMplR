test_that(".formatList different num of items", {
    # more than 5 items
    x_a <- .formatList(c(seq_len(10)))
    expect_equal(x_a, "1, 2 ... 9, 10")

    # 5 items
    x_a <- .formatList(c(seq_len(5)))
    expect_equal(x_a, "1, 2, 3, 4, 5")

    # less than 5 items
    x_a <- .formatList(c(seq_len(3)))
    expect_equal(x_a, "1, 2, 3")
})


test_that(".convertToSNPEffectMatrixCollection given collection", {
    expect_s4_class(
        .convertToSNPEffectMatrixCollection(SEMC),
        "SNPEffectMatrixCollection"
    )
})


test_that(".convertToSNPEffectMatrixCollection given list", {
    x <- getSEMs(SEMC)
    expect_s4_class(
        .convertToSNPEffectMatrixCollection(x),
        "SNPEffectMatrixCollection"
    )
})


test_that(".convertToSNPEffectMatrixCollection given SNPEffectMatrix", {
    x <- getSEMs(SEMC)[[1]]
    expect_s4_class(
        .convertToSNPEffectMatrixCollection(x),
        "SNPEffectMatrixCollection"
    )
})


test_that(".convertToSNPEffectMatrixCollection given invalid input", {
    expect_error(
        .convertToSNPEffectMatrixCollection(TRUE),
        "unable to convert object of class"
    )
    expect_error(
        .convertToSNPEffectMatrixCollection(data.table::data.table()),
        "unable to convert object of class"
    )
})


test_that(".makeVariantId invalid input", {
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009))
    )
    expect_error(.makeVariantId(x),
        regexp = "both refCol and altCol must be defined"
    )
})


test_that(".makeVariantId VRanges input", {
    # SNP
    x <- VariantAnnotation::VRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        ref = "G", alt = "C"
    )
    expect_equal(.makeVariantId(x), "chr12:94136009:G>C")

    # insertions
    x <- VariantAnnotation::VRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        ref = "", alt = "C"
    )
    expect_equal(.makeVariantId(x), "chr12:94136009:insC")

    # deletions
    x <- VariantAnnotation::VRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        ref = "G", alt = ""
    )
    expect_equal(.makeVariantId(x), "chr12:94136009:delG")
})


test_that(".makeVariantId GRanges input", {
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009)),
        ref = "G", alt = "C"
    )
    expect_equal(
        .makeVariantId(x, refCol = "ref", altCol = "alt"),
        "chr12:94136009:G>C"
    )
})


test_that(".makeVariantId GRanges input", {
    stp_a <- .semToPpm(getSEMs(SEMC, "JUN"))

    # make sure rows sum to 1
    cs <- colSums(stp_a) |> round(digits = 5)
    expect_true(all(cs == 1))

    # spot check a few scores
    expect_equal(stp_a[seq_len(4), 2],
        c(
            A = 0.10733802,
            C = 0.08459069,
            G = 0.56368717,
            `T` = 0.24438412
        ),
        tolerance = 1e-6
    )

    expect_equal(stp_a[seq_len(4), 4],
        c(
            A = 0.1785122,
            C = 0.3136177,
            G = 0.3113809,
            `T` = 0.1964891
        ),
        tolerance = 1e-6
    )

    # negative scores are converted to zeros
    expect_equal(stp_a[seq_len(4), 7],
        c(
            A = 1,
            C = 0,
            G = 0,
            `T` = 0
        ),
        tolerance = 1e-6
    )
})
