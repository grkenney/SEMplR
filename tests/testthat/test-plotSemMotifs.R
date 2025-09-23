# setup test example
x <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C",
    id = "1"
)
s <- scoreVariants(
    x = x, sem = sc,
    genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
    varId = "id"
)

test_that(".validatePlotSemMotifsInputs errors on invalid input", {
    expect_no_error(.validatePlotSemMotifsInputs(
        s = s,
        label = "transcription_factor",
        variant = "1",
        cols = c("red", "blue")
    ))
    # invalid label
    expect_error(
        .validatePlotSemMotifsInputs(
            s = s,
            label = "tf",
            variant = "1",
            cols = c("red", "blue")
        ),
        regexp = "label not found"
    )
    # invalid variant
    expect_error(
        .validatePlotSemMotifsInputs(
            s = s,
            label = "transcription_factor",
            variant = "2",
            cols = c("red", "blue")
        ),
        regexp = "variant not found"
    )

    # invalid colors
    expect_error(
        .validatePlotSemMotifsInputs(
            s = s,
            label = "transcription_factor",
            variant = "1",
            cols = c("red")
        ),
        regexp = "cols must be a vector of length 2"
    )
    expect_error(
        .validatePlotSemMotifsInputs(
            s = s,
            label = "transcription_factor",
            variant = "1",
            cols = c("red", "blue", "green")
        ),
        regexp = "cols must be a vector of length 2"
    )
})


test_that("plotSemMotifs generates a ggplot object", {
    plt <- plotSemMotifs(
        s = s, variant = "1",
        label = "transcription_factor"
    )
    expect_s3_class(plt, "gg")
})
