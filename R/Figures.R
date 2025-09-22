# library(SEMplR)
# library(cowplot)
# 
# data(sc)
# 
# # ---- Figure 1A ----
# base_motifs <- c("TGACTCA", "TGAGTCA", "TGATTCA", "TGAATCA", 
#                  "TTAGTCA")
# base_motif <- base_motifs[1]
# scoreBinding(base_motif, sems(sc, "MA0099.2_HeLa"), nFlank = 0)
# bps <- c("A", "C", "G", "T")
# 
# mutated_seqs <- matrix("", nrow = 4, ncol = 7)
# rownames(mutated_seqs) <- bps
# for (i in 1:nchar(base_motif)) {
#   char <- substr(base_motif, i, i)
#   mutations <- setdiff(bps, char)
#   for (m in mutations) {
#     prefix <- substr(base_motif, 0, i-1)
#     suffix <- substr(base_motif, i+1, nchar(base_motif))
#     mutated_seqs[m, i] <- paste0(prefix, m, suffix)
#   }
# }
# 
# mutated_seqs <- as.vector(mutated_seqs) |> unique()
# mutated_seqs <- setdiff(mutated_seqs, "")
# 
# # generate random sequence to flank motif and bring it to 39 bps
# set.seed(1)
# prefix <- simulateRandSeqs(38/2)
# suffix <- simulateRandSeqs(38/2)
# prefix <- "ACACTTTGACCCCCGGATA"
# suffix <- "AACGCACTTTGAGTTCAC"
# flanked_seqs <- paste0(prefix, mutated_seqs, suffix)
# write.csv(flanked_seqs, paste0("~/Desktop/flanked_seqs_", base_motif,".csv"))
# 
# pdf(paste0("~/Desktop/mutated_JUN/",base_motif ,".pdf"), height = 7, width = 5)
# for (i in 1:length(mutated_seqs)) {
#   mutation_score <- scoreBinding(x = flanked_seqs[i],
#                                  semList = sc,
#                                  nFlank = 38/2)
#   base_score <- scoreBinding(x = paste0(prefix, base_motif, suffix),
#                              semList = sc,
#                              nFlank = 38/2)
#   scores_m <- merge(base_score, mutation_score, by = "SEM")
# 
#   dt <- merge(scores_m, semData(sc),
#               by.x = "SEM", by.y = data.table::key(semData(sc)))
#   
#   changedCols <- c("#F8766D", "dodgerblue2")
#   label <- "transcription_factor"
# 
#   plt1 <- plotSEM(sc, motif = "MA0099.2_HeLa",
#                   motifSeq = mutated_seqs[i]) + ggtitle("JUN")
# 
#   plt2 <- ggplot2::ggplot(data = dt,
#                   aes(x = scoreNorm.x, y = scoreNorm.y)) +
#     geom_vline(xintercept = 0, linetype = 2, col = 'grey') +
#     geom_hline(yintercept = 0, linetype = 2, col  = 'grey') +
#     geom_point(data = subset(dt, scoreNorm.y < 0 & scoreNorm.x < 0),
#                size = 1,
#                color = 'grey') +
#     geom_point(data = subset(dt, scoreNorm.y > 0 & scoreNorm.x < 0),
#                size = 1,
#                color = changedCols[1]) +
#     geom_point(data = subset(dt, scoreNorm.y < 0 & scoreNorm.x > 0),
#                size = 1,
#                color = changedCols[2]) +
#     geom_point(data = subset(dt, scoreNorm.y > 0 & scoreNorm.x > 0),
#                size = 1,
#                color = 'grey') +
#     geom_point(data = subset(dt, SEM == "MA0099.2_HeLa"),
#                size = 1,
#                color = "purple") +
#     ggrepel::geom_text_repel(data = subset(dt, scoreNorm.y > 0 & scoreNorm.x < 0),
#                              mapping = aes(label = .data[[label]]),
#                              size = 4,
#                              color = changedCols[1]) +
#     ggrepel::geom_text_repel(data = subset(dt, scoreNorm.y < 0 & scoreNorm.x > 0),
#                              mapping = aes(label = .data[[label]]),
#                              size = 4,
#                              color = changedCols[2]) +
#     scale_x_continuous(breaks = scales::pretty_breaks(),
#                        limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
#     scale_y_continuous(breaks = scales::pretty_breaks(),
#                        limits = \(x) ifelse(abs(x) < 1, c(-1,1), x)) +
#     labs(x = "ref binding propensity",
#          y = "alt binding propensity") +
#     theme_classic() +
#     theme(panel.grid = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text = element_text(size=15)) +
#     ggtitle(flanked_seqs[i])
# 
#   plt <- cowplot::plot_grid(plt1, plt2, ncol = 1,
#                             rel_heights = c(1, 2))
#   print(plt)
# }
# dev.off()
# 
# 
# 
# fig1a <- plotSEM(sc,
#         motif = "MA0099.2_HeLa",
#         motifSeq = "TGAGTAA") +
#   ggtitle("JUN")
# 
# pdf(file = "man/figures/Fig1A.pdf",
#     width = 5, height = 5)
# fig1a
# dev.off()
# 
# # ---- Figure 1B ----
# 
# gr <- GenomicRanges::GRanges(seqnames = chr,
#               ranges = pos,
#               allele = ref_allele,
#               id = caqtls$proxy_rsID)[1000:2000]
# # gr <- gr[gr$id == "rs6665973"]
# sb <- scoreBinding(x = gr, semList = sc,
#                    bs_genome_obj = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
#                    allele = "allele")
# 
# # sb[seqId == "rs6665973" & SEM == "MA0099.2_HeLa"]
# # sb[seqId == "chr1:1543500" & SEM == "MA0099.2_HeLa"]
# sb_oi <- sb[scoreNorm > 0 & SEM == "MA0099.2_HeLa"]
# 
# # 9, 15
# pis <- proxy_ID_split[pos == strsplit(sb_oi$seqId, ":")[[8]][2]]
# vr <- VRanges(seqnames = "chr1",
#               ranges = as.integer(pis[[1]][2]),
#               ref = pis[[1]][3],
#               alt = pis[[1]][4], id = "1")
# sempl_obj <- scoreVariants(vr = vr,
#                            semList = sc,
#                            bs_genome_obj=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
# scores(sempl_obj)[varId == "1" & semId == "MA0099.2_HeLa"]
# fig1b <- plotSemMotifs(sempl_obj,
#                        variant = "1",
#                        label = "transcription_factor")
# caqtls[caqtls$proxy_ID == "chr1:85403115:A:G", ] # rs1357635
# fig1b
# 
# pdf(file = "man/figures/Fig1B.pdf",
#     width = 5, height = 5)
# fig1b
# dev.off()
# 
# 
# 
# # scoreBinding("CTGGAGTCCTACCCACTCTCAGTCAGCTTCCTCCTAACA",
# #              semList = sems(sc, "MA0099.2_HeLa"))
# #
# # vr <- VRanges(seqnames = "chr1",
# #               ranges = roi,
# #               ref = ref_allele[pos == roi],
# #               alt = alt_allele[pos == roi])
# # scoreVariants(vr = vr,
# #               semList = sems(sc, "MA0099.2_HeLa"),
# #               bs_genome_obj=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
# #
# # caqtls[caqtls$proxy_ID == "chr1:1543500:T:G", ]
# # scoreBinding(gr,
# #              semList = sc,
# #              bs_genome_obj=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[SEM == "MA0099.2_HeLa" & seqId == "chr1:1005429"]
# 
# 
# # ---- Figure 1C caQTL ----
# #
# # caqtl_file <- "~/Downloads/caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed"
# # caqtls <- data.table::fread(file = caqtl_file, sep = "\t")
# # proxy_ID_split <- strsplit(caqtls$proxy_ID, ":")
# #
# # # formatting position and alleles
# # chr <- lapply(proxy_ID_split, `[[`, 1) |> unlist()
# # pos <- lapply(proxy_ID_split, `[[`, 2) |> unlist()
# # ref_allele <- lapply(proxy_ID_split, `[[`, 3) |> unlist()
# # alt_allele <- lapply(proxy_ID_split, `[[`, 4) |> unlist()
# 
# # gr <- GRanges(seqnames = chr,
# #               ranges = pos,
# #               allele = ref_allele,
# #               id = caqtls$proxy_rsID)
# #
# # e <- enrichSEMs(gr, semList = sc, bs_genome_obj = bs_genome_obj)
# # fig1c <- plotEnrich(e, sc, threshold = 0.01)
# #
# # pdf(file = "man/figures/Fig1C.pdf",
# #     width = 10, height = 10)
# # fig1c
# # dev.off()
# #
# # saveRDS(e, "~/Documents/UNC/Research/SEMplR/enrich_troubleshooting/results/e_all_caqtls.rds")
# 
# 
# # ---- Figure 1C simulated ----
# 
# ppms_from_sems <- list()
# for (s in sems(sc)) {
#   norm_score <- apply(sem(s), 1,
#                       function(x) (2^x - 2^baseline(s)) / abs(2^baseline(s)))
#   norm_score[norm_score < 0] <- 0
#   ppm <- apply(norm_score, 2, function(x) x / sum(x))
#   ppms_from_sems[[semId(s)]] <- ppm
# }
# 
# 
# # create random sequences weighted by ppm probabilities
# simulatePPMSeqs <- function(ppm, nSeqs) {
#   # make a nested list with the top level being the position in the sequence
#   # and the bottom level being the base for each replicate at that position
#   # ie, when generating 100 seqs for a ppm with 6 positions it will make a
#   # list with length 6 where each element is a list of 100 bases
#   ppm_t <- t(ppm)
#   position_samples <- lapply(1:nrow(ppm_t),
#                              function(i) sample(colnames(ppm_t),
#                                                 size = nSeqs,
#                                                 replace = TRUE,
#                                                 prob = ppm_t[i, ]))
# 
#   rand_pwm_seq_mtx <- position_samples |>
#     unlist() |>
#     matrix(ncol = nrow(ppm_t), nrow = nSeqs)
#   rand_pwm_seqs <- apply(rand_pwm_seq_mtx, 1,
#                          function(x) paste0(x, collapse = ""))
# 
#   return(rand_pwm_seqs)
# }
# 
# 
# simulateRandSeqs <- function(seqLength, nSeqs = 1) {
#   bps <- c("A", "C", "G", "T")
#   seqs <- sample(bps, replace = TRUE,
#                  size = seqLength) |>
#     stringi::stri_c(collapse = "") |>
#     replicate(n = nSeqs)
#   return(seqs)
# }
# 
# ppm <- ppms_from_sems[["MA0099.2_HeLa"]]
# random_jun_seqs <- simulatePPMSeqs(ppm = ppm, nSeqs = 200)
# random_jun_seqs_w_flanks <- lapply(random_jun_seqs,
#                                    function(x) paste0(simulateRandSeqs(seqLength = ncol(ppm)),
#                                                       x,
#                                                       simulateRandSeqs(seqLength = ncol(ppm)))) |>
#   unlist()
# random_seqs <- simulateRandSeqs(seqLength = ncol(ppm)*3, nSeqs = 800)
# 
# sb <- scoreBinding(c(random_jun_seqs_w_flanks, random_seqs),
#                    semList = sc)
# 
# e <- enrichSEMs(sb, semList = sc,
#                 seqs = c(random_jun_seqs_w_flanks, random_seqs))
# 
# fig1c <- plotEnrich(e, sc, threshold = 0.01)
# 
# # ---- Build Figure 1 ----
# 
# pdf(file = "man/figures/Fig1.pdf",
#     width = 13, height = 8)
# left_col <- plot_grid(plt1, plt2, ncol = 1, labels = c('A', 'B'))
# plot_grid(left_col, fig1c, ncol = 2, rel_widths = c(1, 2), labels = c('', 'C'))
# dev.off()
# 
