# Build a contingency table for gained motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
gainedContingencyTable <- function(s, si, lfc=0.5) {
  refNorm <- altNorm <- absLog2FC <- NULL
  
  s$absLog2FC <- log2(abs(s$altNorm) / abs(s$refNorm))
  # Gained in SEM
  a <- s[(refNorm < 0 & altNorm > 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Gained in SEM
  b <- s[!(refNorm < 0 & altNorm > 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Gained not in SEM
  c <- s[(refNorm < 0 & altNorm > 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Gained not in SEM
  d <- s[!(refNorm < 0 & altNorm > 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

# Build a contingency table for lost motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
lostContingencyTable <- function(s, si, lfc=0.5) {
  refNorm <- altNorm <- absLog2FC <- NULL
  
  s$absLog2FC <- log2(abs(s$altNorm) / abs(s$refNorm))
  # Lost in SEM
  a <- s[(refNorm > 0 & altNorm < 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!(refNorm > 0 & altNorm < 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[(refNorm > 0 & altNorm < 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!(refNorm > 0 & altNorm < 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

# Build a contingency table for changed motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
changedContingencyTable <- function(s, si, lfc=0.5) {
  refNorm <- altNorm <- absLog2FC <- NULL
  
  s$absLog2FC <- log2(abs(s$altNorm) / abs(s$refNorm))
  # Lost in SEM
  a <- s[((refNorm > 0 & altNorm < 0) |
           (refNorm < 0 & altNorm > 0) & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!((refNorm > 0 & altNorm < 0) |
             (refNorm < 0 & altNorm > 0) & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[((refNorm > 0 & altNorm < 0) |
            (refNorm < 0 & altNorm > 0) & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!((refNorm > 0 & altNorm < 0) |
             (refNorm < 0 & altNorm > 0) & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

# Return a contingency table according to direction of test
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param d direction of test. Options are: "changed", "gained", "lost"
# @param lfc Log2FC cutoff (default = 0.5)
# 
#' @keywords internal
getContingencyTable <- function(s, si, d, lfc) {
  if (d == "changed") {
    ct <- changedContingencyTable(s = s, si = si, lfc=0.5)
  } else if (d == "gained") {
    ct <- gainedContingencyTable(s = s, si = si, lfc=0.5)
  } else if (d == "lost") {
    ct <- lostContingencyTable(s = s, si = si, lfc=0.5)
  } else {
    stop("Invalid input for argument d. Must be 'changed', 'gained', or 'lost'")
  }
  return(ct)
}

#' Calculate enrichment of motifs broken among variant set
#' 
#' Performs a Fisher's exact test to test whether a motif is changed in the
#' specified direction of change more often than expected within the set of 
#' variants.
#'
#' @param semScores a SemplScores object
#' @param d direction of test. Options are: "changed", "gained", "lost"
#' @param lfc Log2FC cutoff (default = 0.5)
#' 
#' @return a data.table with results of the enrichment test
#' - `n.changed`: number of motifs changed in direction specified
#' - `odds.ratio`: odds ratio of fisher's test
#' - `ci.lower`: lower 95% confidence interval
#' - `ci.upper`: upper 95% confidence interval
#' - `pvalue`: p-value of Fisher's test
#' - `adj.pvalue`: Benjamini & Hochberg adjusted pvalue
#' 
#' @keywords internal
#' @export
#' 
#' @examples
#' library(VariantAnnotation)
#'
#' # create an SNP Effect Matrix (SEM)
#' sem <- matrix(rnorm(12), ncol = 4)
#' colnames(sem) <- c("A", "C", "G", "T")
#' 
#' # create a VRanges object
#' vr <- VRanges(seqnames = "chr12",
#'               ranges = 94136009, 
#'               ref = "G", alt = "C")
#' 
#' # create a list of SNPEffectMatrix objects
#' semList <- SNPEffectMatrix(sem, baseline = 0.5, semId = "sem_id")
#' 
#' # calculate binding propensity
#' s <- scoreVariants(vr, semList, 
#'                    bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#' 
#' # calculate enrichment
#' enrichSEMs(s, d = "changed")
#' 
enrichSEMs <- function(semScores, d="changed", lfc=0.5) {
  adj.pvalue <- NULL
  
  s <- scores(semScores)
  
  sem_data_cols <- semData(semScores) |> colnames()
  fisher_scores <- matrix(numeric(), 
                          ncol=7+length(sem_data_cols), 
                          nrow=length(unique(s$semId)))
  fisher_scores <- data.table(fisher_scores)
  colnames(fisher_scores) <- c("semId", sem_data_cols, 
                               "n.changed", "odds.ratio",
                               "ci.lower", "ci.upper",
                               "pvalue", "adj.pvalue")
  data.table::setkey(fisher_scores, semId)
  fisher_scores$semId <- unique(s$semId)
  fisher_scores[, sem_data_cols] <- semData(semScores)
  
  for (i in 1:nrow(fisher_scores)) {
    si <- fisher_scores$semId[i]
    ct <- getContingencyTable(s, si, d, lfc)
    f <- stats::fisher.test(ct)
    
    fisher_scores[i, "n.changed"] <- ct[1,1]
    fisher_scores[i, "pvalue"] <- f$p.value
    fisher_scores[i, "odds.ratio"] <- f$estimate
    fisher_scores[i, "ci.lower"] <- f$conf.int[1]
    fisher_scores[i, "ci.upper"] <- f$conf.int[2]
  }
  fisher_scores$adj.pvalue <- stats::p.adjust(fisher_scores$pvalue, 
                                              method = "BH")
  fisher_scores <- fisher_scores[order(adj.pvalue)]
  return(fisher_scores)
}

