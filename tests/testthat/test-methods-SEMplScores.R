test_that("SEMScores constructs a blank object and accessors work", {
    x <- SEMScores()
    expect_equal(getRanges(x), VariantAnnotation::VRanges())
    expect_equal(semData(x), data.table::data.table())
    expect_equal(scores(x), data.table::data.table())
})


test_that("scores can be updated", {
    x <- SEMScores()
    scores(x) <- data.table()
    expect_equal(x, x)
})


test_that("show SEMScores prints empty object", {
    x <- SEMScores()
    expect_output(print(x), "An object of class SEMScores")
    expect_output(print(x), "ranges\\(0\\):")
    expect_output(print(x), "semData\\(0\\):")
    expect_output(print(x), "scores\\(0\\):")
})


test_that("show SEMScores prints scored object", {
    x <- GenomicRanges::GRanges(
        seqnames = c("chr12"),
        ranges = IRanges::IRanges(c(94136009))
    )
    sb_a <- scoreBinding(
        x = x,
        sem = SEMC,
        genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
        rc = FALSE
    )
    expect_output(print(sb_a), "An object of class SEMScores")
    expect_output(print(sb_a), "ranges\\(1\\): chr12:94136009")
    sem_data_str <- "semData\\(12\\): transcription_factor, ensembl_id ... "
    expect_output(print(sb_a), sem_data_str)
    expect_output(print(sb_a), "scores\\(223\\):")
})
