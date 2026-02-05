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
    x <- getSEMs(SEMC)[["TFAP2B"]]
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


orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that(".poolFilter minimal example", {
    ids <- .mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        id_type = "ENSEMBL")
    pf_a <- .poolFilter(ids, geneType = "protein-coding")
    pf_e <- data.frame(ENSEMBL = my_genes, ENTREZID = c("675", "673"), 
                       GENETYPE = "protein-coding")
    expect_equal(pf_a$fg_ids, pf_e)
})


test_that(".poolFilter invalid input", {
    ids <- .mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        id_type = "ENSEMBL")
    expect_error(.poolFilter(ids, geneType = "foo"), 
                 regexp = "Invalid geneType 'foo'")
})


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that(".getCoordinates minimal example", {
    ids <- .mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        background_ids = c("ENSG00000109321", "ENSG00000078061"),
        id_type = "ENSEMBL")
    pf <- .poolFilter(ids, geneType = "protein-coding")
    coords_a <- .getCoordinates(mapped = pf, txdb = txdb, n_ratio = 1)
    expect_in(as.vector(GenomicRanges::seqnames(coords_a$foregroundElements)), 
              c("chr13", "chr7"))
    expect_in(GenomicRanges::start(coords_a$foregroundElements), 
              c(32314786, 140924880))
    expect_in(as.vector(GenomicRanges::seqnames(coords_a$backgroundElements)), 
              c("chrX", "chr4"))
    expect_in(GenomicRanges::start(coords_a$backgroundElements), 
              c(47560905, 74444836))
})


test_that(".getCoordinates invalid input", {
    ids <- .mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        background_ids = "ENSG00000139610",
        id_type = "ENSEMBL")
    pf <- .poolFilter(ids, geneType = "protein-coding")
    expect_error(.getCoordinates(mapped = pf, txdb = txdb),
                 regexp = "Requested 2 background regions")
})


orgdb <- org.Hs.eg.db::org.Hs.eg.db
my_genes <- c("ENSG00000139618", "ENSG00000157764")


test_that(".mapIds minimal example", {
    ids_a <- .mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        id_type = "ENSEMBL")
    fg_ids_e <- data.frame(ENSEMBL = my_genes, ENTREZID = c("675", "673"), 
                           GENETYPE = "protein-coding")
    expect_equal(ids_a$fg_ids, fg_ids_e)
})


test_that(".mapIds errors and messages", {
    expect_error(.mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        id_type = "foo"), 
        regexp = "not an available key for in this mapping object")
    expect_message(.mapIDs(
        orgdb = orgdb,
        foreground_ids = my_genes,
        background_ids = "ENSG00000139610",
        id_type = "ENSEMBL"), 
        regexp = "Mapping background ids")
})
