# Build a contingency table for gained motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
# 
# @keywords internal
gainedContingencyTable <- function(s, si, lfc=0.5) {
  s$absLog2FC <- log2(abs(s$riskNorm) / abs(s$nonRiskNorm))
  # Gained in SEM
  a <- s[(nonRiskNorm < 0 & riskNorm > 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Gained in SEM
  b <- s[!(nonRiskNorm < 0 & riskNorm > 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Gained not in SEM
  c <- s[(nonRiskNorm < 0 & riskNorm > 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Gained not in SEM
  d <- s[!(nonRiskNorm < 0 & riskNorm > 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

# Build a contingency table for lost motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
# 
# @keywords internal
lostContingencyTable <- function(s, si, lfc=0.5) {
  s$absLog2FC <- log2(abs(s$riskNorm) / abs(s$nonRiskNorm))
  # Lost in SEM
  a <- s[(nonRiskNorm > 0 & riskNorm < 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!(nonRiskNorm > 0 & riskNorm < 0 & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[(nonRiskNorm > 0 & riskNorm < 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!(nonRiskNorm > 0 & riskNorm < 0 & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

# Build a contingency table for changed motifs
#
# @param s score table from a SemplScores object
# @param si SEM id
# @param lfc Log2FC cutoff (default = 0.5)
# 
# @keywords internal
changedContingencyTable <- function(s, si, lfc=0.5) {
  s$absLog2FC <- log2(abs(s$riskNorm) / abs(s$nonRiskNorm))
  # Lost in SEM
  a <- s[((nonRiskNorm > 0 & riskNorm < 0) |
           (nonRiskNorm < 0 & riskNorm > 0) & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!((nonRiskNorm > 0 & riskNorm < 0) |
             (nonRiskNorm < 0 & riskNorm > 0) & abs(absLog2FC) > lfc) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[((nonRiskNorm > 0 & riskNorm < 0) |
            (nonRiskNorm < 0 & riskNorm > 0) & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!((nonRiskNorm > 0 & riskNorm < 0) |
             (nonRiskNorm < 0 & riskNorm > 0) & abs(absLog2FC) > lfc) & semId != si] |> nrow()
  
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
# @keywords internal
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
#' s <- semMotifBinding(vr, semList)
#' 
#' # calculate enrichment
#' enrichment(s, d = "changed")
#' 
enrichment <- function(semScores, d="changed", lfc=0.5) {
  s <- scores(semScores)
  
  fisher_scores <- matrix(numeric(), ncol=8, nrow=length(unique(s$semId)))
  fisher_scores <- data.table(fisher_scores)
  colnames(fisher_scores) <- c("semId", "tf", "n.changed", "odds.ratio",
                               "ci.lower", "ci.upper",
                               "pvalue", "adj.pvalue")
  fisher_scores$semId <- unique(s$semId)
  fisher_scores$tf <- metadata(semScores)$tf[match(fisher_scores$semId, 
                                                   metadata(semScores)$semId)]
  
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

