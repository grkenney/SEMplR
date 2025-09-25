test_that(".validateSEM fails on invalid input", {
    # missing a column
    m <- matrix(rnorm(12), nrow = 4)
    colnames(m) <- c("A", "C", "G")
    expect_error(.validateSEM(m, "someFile.tsv"),
        regexp = "3 columns detected in file someFile.tsv"
    )

    # unexpected column
    m <- matrix(rnorm(16), nrow = 4)
    colnames(m) <- c("A", "C", "G", "N")
    expect_error(.validateSEM(m, "someFile.tsv"),
        regexp = "Unexpected column\\(s\\)"
    )
})


test_that(".extractBaseline fails when baseline not in header", {
    m <- matrix(rnorm(16), nrow = 4)
    colnames(m) <- c("A", "C", "G", "T")
    tf <- tempfile()
    write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)
    expect_error(.extractBaseline(tf),
        regexp = "No baseline given."
    )
})


test_that(".extractBaseline extracts correct baseline", {
    m <- matrix(rnorm(16), nrow = 4)
    colnames(m) <- c("A", "C", "G", "T")
    tf <- tempfile()
    baseline_str <- "#BASELINE:1.25"
    m_str <- paste0(m[1, ], collapse = "\t")
    write(
        x = c(baseline_str, m_str),
        file = tf
    )
    bl_a <- .extractBaseline(tf)
    expect_equal(bl_a, "1.25")
})


test_that(".loadSEM loads sem file", {
    m <- matrix(rnorm(16), nrow = 4)
    colnames(m) <- c("A", "C", "G", "T")
    tf <- tempfile()
    write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)

    # without semId
    sem <- .loadSEM(tf, bl = 1.25)
    expect_s4_class(sem, "SNPEffectMatrix")
    expect_equal(getBaseline(sem), 1.25)
    expect_equal(getSEMId(sem), basename(tf))

    # with semId
    sem <- .loadSEM(tf, semId = "foo", bl = 1.25)
    expect_s4_class(sem, "SNPEffectMatrix")
    expect_equal(getBaseline(sem), 1.25)
    expect_equal(getSEMId(sem), "foo")


    m <- matrix(rnorm(16), nrow = 4)
    m_str <- paste0(m[1, ], collapse = "\t")
    m_cols <- paste0(c("A", "C", "G", "T"), collapse = "\t")
    tf <- tempfile()
    baseline_str <- "#BASELINE:1.25"
    write(
        x = c(baseline_str, m_cols, m_str),
        file = tf
    )

    # without baseline
    sem <- .loadSEM(semFile = tf, semId = "foo")
    expect_s4_class(sem, "SNPEffectMatrix")
    expect_equal(getBaseline(sem), 1.25)
    expect_equal(getSEMId(sem), "foo")
})


test_that("loadSEMCollection loads first SEM correctly", {
    one_sem <- getSEMs(SEMC)[[1]]
    one_sem_meta <- semData(SEMC)[1, ]

    # write data to file
    tf <- tempfile()
    write.table(
        x = getSEM(one_sem), file = tf,
        sep = "\t", row.names = FALSE
    )

    sem_col_e <- loadSEMCollection(
        semFiles = tf,
        semMetaData = one_sem_meta,
        semMetaKey = "SEM_KEY",
        semIds = one_sem_meta$SEM_KEY,
        bls = one_sem_meta$SEM_baseline
    )
    # check each slot
    expect_equal(getSEMs(sem_col_e)[[1]], one_sem)
    expect_equal(semData(sem_col_e), one_sem_meta)
    expect_equal(sem_col_e@semKey, "SEM_KEY")
})


test_that("loadSEMCollection errors on missing key with metadata", {
    one_sem <- getSEMs(SEMC)[[1]]
    one_sem_meta <- semData(SEMC)[1, ]

    # write data to file
    tf <- tempfile()
    write.table(
        x = getSEM(one_sem), file = tf,
        sep = "\t", row.names = FALSE
    )

    expect_error(
        loadSEMCollection(
            semFiles = tf,
            semMetaData = one_sem_meta,
            semMetaKey = "",
            semIds = one_sem_meta$SEM_KEY,
            bls = one_sem_meta$SEM_baseline
        ),
        regexp = "If providing semMetaData, must specify a column name"
    )
})
