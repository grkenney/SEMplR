# setup enrichment
x <- GenomicRanges::GRanges(
    seqnames = c("chr12", "chr16"),
    ranges = IRanges::IRanges(
        start = c(
            94136009,
            67637901
        ),
        end = c(
            94136009,
            67637901
        )
    )
)
b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
sb <- scoreBinding(x = x, sem = SEMC, genome = b)
e <- enrichSEMs(x = sb, sem = SEMC, genome = b)


test_that(".formatMotifs adds correct altname", {
    motifs <- .formatMotifs(
        sem = SEMC,
        label = "transcription_factor"
    )
    expect_equal(motifs[[1]]["altname"], "TFAP2B")
    expect_equal(motifs[[1]]["name"], "TFAP2B")
    expect_equal(motifs[[100]]["altname"], "GATA1")
    expect_equal(motifs[[100]]["name"], "GATA1")
})


test_that(".constructComparisons returns correctly sized object", {
    motifs <- .formatMotifs(
        sem = SEMC,
        label = "transcription_factor"
    )
    labels <- lapply(motifs, function(x) x["altname"]) |> unlist()

    comps <- .constructComparisons(
        motifs = motifs, labels = labels,
        method = "WPCC"
    )
    expect_s3_class(comps, "hclust")
})
