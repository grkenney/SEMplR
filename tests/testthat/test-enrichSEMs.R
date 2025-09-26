test_that(".scrambleSeqs retains same frequency of bases", {
    dna_seq <- "AAAATGGGG"
    base_freq_e <- strsplit(dna_seq, split = "")[[1]] |>
        table()

    scrabled_seq <- .scrambleSeqs(dna_seq)
    base_freq_a <- strsplit(scrabled_seq, split = "")[[1]] |>
        table()

    expect_equal(base_freq_a, base_freq_a)
})


test_that(".scrambleSeqs retains same frequency of bases", {
    dna_seq <- "AAAATGGGG"
    base_freq_e <- strsplit(dna_seq, split = "")[[1]] |>
        table()

    scrabled_seq <- .scrambleSeqs(dna_seq)
    base_freq_a <- strsplit(scrabled_seq, split = "")[[1]] |>
        table()

    expect_equal(base_freq_a, base_freq_a)
})


test_that(".binomSEMEnrich produces expected stats", {
    # typical case
    xs <- data.table::data.table(
        seqId = c("chr1:1", "chr1:2"),
        SEM = "sem_id",
        score = c(-2, 0),
        scoreNorm = c(-1, 1),
        index = 1,
        seq = "AAAA"
    )
    bg <- data.table::data.table(
        seqId = c("chr2:1", "chr2:2"),
        SEM = "sem_id",
        score = c(-2, -2),
        scoreNorm = c(-1, -1),
        index = 1,
        seq = "AAAA"
    )
    res_a <- .binomSEMEnrich(xs = xs, bg = bg, semName = "sem_id")
    res_e <- data.table::data.table(
        SEM = "sem_id",
        pvalue = 0.75,
        n_bound = 1,
        n_bound_bg = 0
    )
    expect_equal(res_a, res_e)


    # where prob of success is > 1 because of pseudo counts
    bg <- data.table::data.table(
        seqId = c("chr2:1", "chr2:2"),
        SEM = "sem_id",
        score = c(0, 0),
        scoreNorm = c(1, 1),
        index = 1,
        seq = "AAAA"
    )
    res_a <- .binomSEMEnrich(xs = xs, bg = bg, semName = "sem_id")
    res_e <- data.table::data.table(
        SEM = "sem_id",
        pvalue = 1,
        n_bound = 1,
        n_bound_bg = 2
    )
    expect_equal(res_a, res_e)
})


test_that("enrichSEMs on 2 sites without specifying background", {
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
    sb_a <- scoreBinding(
        x = x,
        sem = SEMC,
        genome = b
    )
    e_a <- enrichSEMs(
        x = sb_a, sem = SEMC,
        genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    )

    # all SEMs are represented
    expect_equal(table(e_a$SEM), table(semData(SEMC)$transcription_factor))
    # n_bound is correct
    expect_equal(sum(e_a$n_bound), 2)
    expect_equal(
        e_a[e_a$n_bound > 0, ]$SEM,
        c("IKZF1", "ZNF18")
    )
    expect_equal(
        e_a[e_a$n_bound > 0, ]$n_bound,
        c(1, 1)
    )
})


test_that("enrichSEMs on 2 sites with specifying background", {
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
    bg <- GenomicRanges::GRanges(
        seqnames = c("chr6", "chr6"),
        ranges = IRanges::IRanges(
            start = c(
                167127770,
                135097880
            ),
            end = c(
                167127770,
                135097880
            )
        )
    )
    b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    # score background seperately to validate results
    sb_bg <- scoreBinding(
        x = bg,
        sem = SEMC,
        genome = b
    )

    sb_a <- scoreBinding(
        x = x,
        sem = SEMC,
        genome = b
    )

    e_a <- enrichSEMs(
        x = sb_a, sem = SEMC, background = bg,
        genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    )

    # all SEMs are represented
    expect_equal(table(e_a$SEM), table(semData(SEMC)$transcription_factor))
    # n_bound is correct
    expect_equal(sum(e_a$n_bound), sum(scores(sb_a)$scoreNorm > 0))
    expect_equal(sum(e_a$n_bound), sum(scores(sb_bg)$scoreNorm > 0))

    expect_equal(
        e_a[e_a$n_bound > 0, ]$SEM,
        c("IKZF1", "ZNF18")
    )
    expect_equal(
        e_a[e_a$n_bound > 0, ]$n_bound,
        c(1, 1)
    )

    expect_equal(
        e_a[e_a$n_bound_bg > 0, ]$SEM,
        c("IKZF1", "HNF4A:NR2F2")
    )
    expect_equal(
        e_a[e_a$n_bound_bg > 0, ]$n_bound_bg,
        c(1, 1)
    )
})
