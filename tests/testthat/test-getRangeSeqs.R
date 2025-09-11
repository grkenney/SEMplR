b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

test_that(".querySeq pulls correct seq", {
  # VRange
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  seq_a <- .querySeq(vr, b = b)
  seq_e <- "G"
  expect_equal(seq_a, seq_e)
  
  # GRange
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  seq_a <- .querySeq(gr, b = b)
  seq_e <- "G"
  expect_equal(seq_a, seq_e)
})


test_that(".querySeq pulls correct seq with adjusted start/ends", {
  # VRange
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  seq_a <- .querySeq(vr, b = b, sp = 94136009-2, ep = 94136009+2)
  expect_equal(seq_a, "TTGAG")
  
  # GRange
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  seq_a <- .querySeq(gr, b = b, sp = 94136009-2, ep = 94136009+2)
  expect_equal(seq_a, "TTGAG")
})


test_that(".checkAlleles catches matching allele to ref", {
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  # don't expect a warning
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = "A"), 
                 regexp = NA)
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = c("A", "C")), 
                 regexp = NA)
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = c("A", "")), 
                 regexp = NA)
  expect_warning(.checkAlleles(x = gr, qs = "AG", allele = c("A", "AG")), 
                 regexp = NA)
  
  # expect warning
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = "T"), 
                 regexp = "does not match reference sequence")
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = c("T", "C")), 
                 regexp = "does not match reference sequence")
  expect_warning(.checkAlleles(x = gr, qs = "A", allele = ""), 
                 regexp = "does not match reference sequence")
})


test_that(".checkAlleles stops if invalid number of alleles", {
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  # expect error
  expect_error(.checkAlleles(x = gr, qs = "A", allele = c("A", "G", "C")), 
               regexp = "The number of alleles must be 2 or less")
})


test_that(".buildSeq builds sequence from range alone", {
  # GRanges
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  bs_a <- .buildSeq(x = gr, up = 1, down = 1, b = b)
  bs_e <- c(sequence = "TGA")
  expect_equal(bs_a, bs_e)
  
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = IRanges::IRanges(94136009, 94136010))
  bs_a <- .buildSeq(x = gr, up = 1, down = 1, b = b)
  bs_e <- c(sequence = "TGAG")
  expect_equal(bs_a, bs_e)
  
  # VRanges
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  bs_a <- .buildSeq(x = vr, up = 1, down = 1, b = b)
  bs_e <- c(sequence = "TGA")
  expect_equal(bs_a, bs_e)
})


test_that(".buildSeq builds sequence from range with alleles", {
  # GRanges
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009)
  bs_a <- .buildSeq(x = gr, up = 1, down = 1, b = b, allele = "G")
  bs_e <- c(ref_seq = "TGA", alt_seq = NA)
  expect_equal(bs_a, bs_e)
  
  bs_a <- .buildSeq(x = gr, up = 1, down = 1, b = b, allele = c("G", "C"))
  bs_e <- c(ref_seq = "TGA", alt_seq = "TCA")
  expect_equal(bs_a, bs_e)
  
  # VRanges
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  bs_a <- .buildSeq(x = vr, up = 1, down = 1, b = b, allele = c("G", "C"))
  bs_e <- c(ref_seq = "TGA", alt_seq = "TCA")
  expect_equal(bs_a, bs_e)
  
  bs_a <- .buildSeq(x = vr, up = 1, down = 1, b = b, allele = c("G", ""))
  bs_e <- c(ref_seq = "TGA", alt_seq = "TA")
  expect_equal(bs_a, bs_e)
})


test_that(".getAllele gets alleles from GRanges object", {
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009,
                               ref = "G", alt = "C")
  ga_a <- .getAllele(x = gr, refCol = "ref", altCol = "alt")
  ga_e <- list(ref = "G", alt = "C")
  expect_equal(ga_a, ga_e)
  
  ga_a <- .getAllele(x = gr, refCol = "ref")
  ga_e <- list(ref = "G", alt = NULL)
  expect_equal(ga_a, ga_e)
  
  ga_a <- .getAllele(x = gr, altCol = "alt")
  ga_e <- list(ref = NULL, alt = "C")
  expect_equal(ga_a, ga_e)
})


test_that(".getAllele gets alleles from VRanges object", {
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  ga_a <- .getAllele(x = vr)
  ga_e <- list(ref = "G", alt = "C")
  expect_equal(ga_a, ga_e)
  
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "")
  ga_a <- .getAllele(x = vr)
  ga_e <- list(ref = "G", alt = "")
  expect_equal(ga_a, ga_e)
})


test_that(".getAllele check errors and messages", {
  gr <- GenomicRanges::GRanges(seqnames = "chr12",
                               ranges = 94136009,
                               ref = "G", alt = "C")
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "C")
  expect_error(.getAllele(x = gr, refCol = "ref", altCol = "foo"))
  expect_error(.getAllele(x = gr, refCol = "foo", altCol = "alt"))
  expect_message(.getAllele(x = vr, refCol = "A", altCol = "T"),
                 regexp = "ignoring allele column")
})


test_that("getRangeSeqs populates GRange mcols", {
  # 2bp range without alleles
  gr <- GenomicRanges::GRanges("chr12", 
                               IRanges::IRanges(94136009, 94136010))
  gfs_a <- getRangeSeqs(gr, up = 1, down = 1, genome = b)
  gfs_e <- GenomicRanges::GRanges("chr12", 
                                  IRanges::IRanges(94136009, 94136010),
                                  sequence = "TGAG")
  expect_equal(gfs_a, gfs_e)
  
  # 1bp range with alleles, but not defined in function
  gr <- GenomicRanges::GRanges("chr12", IRanges::IRanges(94136009),
                               ref = "G", alt = "C")
  gfs_a <- getRangeSeqs(gr, up = 1, down = 1, genome = b)
  gfs_e <- GenomicRanges::GRanges("chr12", IRanges::IRanges(94136009),
                                  ref = "G", alt = "C", sequence = "TGA")
  expect_equal(gfs_a, gfs_e)
  
  # 1bp range with alleles and defined in function
  gfs_a <- getRangeSeqs(gr, up = 1, down = 1, genome = b,
                           refCol = "ref", altCol = "alt")
  gfs_e <- GenomicRanges::GRanges("chr12", IRanges::IRanges(94136009),
                                  ref = "G", alt = "C", 
                                  ref_seq = "TGA", alt_seq = "TCA")
  expect_equal(gfs_a, gfs_e)
  
  # 1bp range with alleles with only ref defined in function
  gfs_a <- getRangeSeqs(gr, up = 1, down = 1, genome = b,
                           refCol = "ref")
  alleles <- c("TGA", NA)
  gfs_e <- GenomicRanges::GRanges("chr12", IRanges::IRanges(94136009),
                                  ref = "G", alt = "C", 
                                  ref_seq = alleles[1], alt_seq = alleles[2])
  expect_equal(gfs_a, gfs_e)
})


test_that("getRangeSeqs populates VRange mcols", {
  vr <- VariantAnnotation::VRanges(seqnames = "chr12",
                                   ranges = 94136009,
                                   ref = "G", alt = "")
  
  # 1bp range with alleles
  gfs_a <- getRangeSeqs(vr, up = 1, down = 1, genome = b)
  gfs_e <- VariantAnnotation::VRanges(seqnames = "chr12",
                                      ranges = 94136009,
                                      ref = "G", alt = "", 
                                      ref_seq = "TGA", alt_seq = "TA")
  expect_equal(gfs_a, gfs_e)
})


test_that("getRangeSeqs populates multiple ranges", {
  range_oi <- IRanges::IRanges(start = c(94136009, 94136010),
                      end = c(c(94136009, 94136010)))
  # GRanges
  gr <- GenomicRanges::GRanges(seqnames = c("chr12", "chr12"),
                               ranges = range_oi,
                               ref = c("G", "A"), 
                               alt = c("C", ""))
  gfs_a <- getRangeSeqs(gr, up = 1, down = 1, genome = b,
                           refCol = "ref", altCol = "alt")
  gfs_e <- GenomicRanges::GRanges(seqnames = c("chr12", "chr12"),
                                  ranges = range_oi,
                                  ref = c("G", "A"), 
                                  alt = c("C", ""),
                                  ref_seq = c("TGA", "GAG"), 
                                  alt_seq = c("TCA", "GG"))
  expect_equal(gfs_a, gfs_e)
  
  # VRanges
  vr <- VariantAnnotation::VRanges(seqnames = c("chr12", "chr12"),
                                   ranges = range_oi,
                                   ref = c("G", "A"), 
                                   alt = c("C", "T"))
  
  # 1bp range with alleles
  gfs_a <- getRangeSeqs(vr, up = 1, down = 1, genome = b)
  gfs_e <- VariantAnnotation::VRanges(seqnames = c("chr12", "chr12"),
                                      ranges = range_oi,
                                      ref = c("G", "A"), 
                                      alt = c("C", "T"),
                                      ref_seq = c("TGA", "GAG"), 
                                      alt_seq = c("TCA", "GTG"))
  expect_equal(gfs_a, gfs_e)
})


test_that("getRangeSeqs test invalid ranges", {
  # Empty GRanges
  gr <- GenomicRanges::GRanges()
  expect_error(getRangeSeqs(x = gr, genome = b), 
               regexp = "x must contain at least one range")
  
  # Not GRanges or VRanges
  expect_error(getRangeSeqs(x = "foo", genome = b), 
               regexp = "x must be of class VRanges or GRanges")
  
})
