#' Build a contingency table for gained motifs
#'
#' @param s score table from a SemplScores object
#' @param si SEM id
gainedContingencyTable <- function(s, si) {
  # Gained in SEM
  a <- s[(nonRiskNorm < 0 & riskNorm > 0) & semId == si] |> nrow()
  # Not Gained in SEM
  b <- s[!(nonRiskNorm < 0 & riskNorm > 0) & semId == si] |> nrow()
  # Gained not in SEM
  c <- s[(nonRiskNorm < 0 & riskNorm > 0) & semId != si] |> nrow()
  # Not Gained not in SEM
  d <- s[!(nonRiskNorm < 0 & riskNorm > 0) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

#' Build a contingency table for lost motifs
#'
#' @param s score table from a SemplScores object
#' @param si SEM id
lostContingencyTable <- function(s, si) {
  # Lost in SEM
  a <- s[(nonRiskNorm > 0 & riskNorm < 0) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!(nonRiskNorm > 0 & riskNorm < 0) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[(nonRiskNorm > 0 & riskNorm < 0) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!(nonRiskNorm > 0 & riskNorm < 0) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

#' Build a contingency table for changed motifs
#'
#' @param s score table from a SemplScores object
#' @param si SEM id
changedContingencyTable <- function(s, si) {
  # Lost in SEM
  a <- s[((nonRiskNorm > 0 & riskNorm < 0) |
           (nonRiskNorm < 0 & riskNorm > 0)) & semId == si] |> nrow()
  # Not Lost in SEM
  b <- s[!((nonRiskNorm > 0 & riskNorm < 0) |
             (nonRiskNorm < 0 & riskNorm > 0)) & semId == si] |> nrow()
  # Lost not in SEM
  c <- s[((nonRiskNorm > 0 & riskNorm < 0) |
            (nonRiskNorm < 0 & riskNorm > 0)) & semId != si] |> nrow()
  # Not Lost not in SEM
  d <- s[!((nonRiskNorm > 0 & riskNorm < 0) |
             (nonRiskNorm < 0 & riskNorm > 0)) & semId != si] |> nrow()
  
  contingencyTable <- matrix(c(a, c, 
                               b, d), ncol=2)
  return(contingencyTable)
}

#' Return a contingency table according to direction of test
#'
#' @param s score table from a SemplScores object
#' @param si SEM id
#' @param d direction of test. Options are: "changed", "gained", "lost"
getContingencyTable <- function(s, si, d) {
  if (d == "changed") {
    ct <- changedContingencyTable(s = s, si = si)
  } else if (d == "gained") {
    ct <- gainedContingencyTable(s = s, si = si)
  } else if (d == "lost") {
    ct <- lostContingencyTable(s = s, si = si)
  } else {
    stop("Invalid input for argument d. Must be 'changed', 'gained', or 'lost'")
  }
  return(ct)
}

#' Return a contingency table according to direction of test
#'
#' @param semScores a SemplScores object
#' @param d direction of test. Options are: "changed", "gained", "lost"
#' 
#' @export
enrichment <- function(semScores, d="changed") {
  s <- scores(semScores)
  
  fisher_scores <- matrix(numeric(), ncol=7, nrow=length(unique(s$semId)))
  fisher_scores <- data.table(fisher_scores)
  colnames(fisher_scores) <- c("semId", "tf", "odds.ratio",
                               "ci.lower", "ci.upper",
                               "pvalue", "adj.pvalue")
  fisher_scores$semId <- unique(s$semId)
  fisher_scores$tf <- metadata(semScores)$tf[match(fisher_scores$semId, 
                                                   metadata(semScores)$semId)]
  
  for (i in 1:nrow(fisher_scores)) {
    si <- fisher_scores$semId[i]
    ct <- getContingencyTable(s, si, d)
    f <- fisher.test(ct)
    
    fisher_scores[i, "pvalue"] <- f$p.value
    fisher_scores[i, "odds.ratio"] <- f$estimate
    fisher_scores[i, "ci.lower"] <- f$conf.int[1]
    fisher_scores[i, "ci.upper"] <- f$conf.int[2]
  }
  fisher_scores$adj.pvalue <- p.adjust(fisher_scores$pvalue, 
                                       method = "BH")
  return(fisher_scores)
}

