# setup test example
x <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C",
    id = "1"
)
s <- scoreVariants(
    x = x, sem = SEMC,
    genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
    varId = "id"
)

test_that(".validatePlotSemVariantsInputs errors on invalid input", {
    # no error
    expect_no_error(.validatePlotSemVariantsInputs(
        s = s,
        label = "id",
        semId = "JUN",
        cols = c("red", "blue")
    ))

    # invalid label
    expect_error(
        .validatePlotSemVariantsInputs(
            s = s,
            label = "foo",
            semId = "JUN",
            cols = c("red", "blue")
        ),
        regexp = "label not found"
    )

    # invalid semId
    expect_error(
        .validatePlotSemVariantsInputs(
            s = s,
            label = "id",
            semId = "foo",
            cols = c("red", "blue")
        ),
        regexp = "variant not found"
    )

    # invalid cols
    expect_error(
        .validatePlotSemVariantsInputs(
            s = s,
            label = "id",
            semId = "JUN",
            cols = c("red")
        ),
        regexp = "cols must be a vector of length 2"
    )
    expect_error(
        .validatePlotSemVariantsInputs(
            s = s,
            label = "id",
            semId = "JUN",
            cols = c("red", "blue", "green")
        ),
        regexp = "cols must be a vector of length 2"
    )
})


test_that("plotSEMVariants returns a ggplot object", {
    plt <- plotSEMVariants(
        s = s,
        sem = "JUN"
    )
    expect_s3_class(plt, "gg")
})
