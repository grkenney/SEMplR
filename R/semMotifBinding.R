#' Define a function to calculate risk/non-risk binding propensity
#' for all SEM motifs provided
#'
#' @param chr character of chromosome
#' @param pos integer of SNP position
#' @param ref character of reference nucleotide (can be upper or lower)
#' @param alt character of alternate nulceotide (can be upper or lower)
#' @param semFiles character vector of paths to each `.sem` file
#' @param semBaselines file path to sem baselines
#'
semMotifBinding <- \(chr = "chr12",
                     pos = 94136009,
                     ref = "G",
                     alt = "C",
                     semFiles,
                     semBaselines) {

  ## Load SEM files ##
  ## Load libraries
  library(data.table)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg19)

  ## Read in SEM files
  sems <- lapply(semFiles, fread, header = TRUE)

  ## Add names to sems
  names(sems) <-
    lapply(sems, \(x) colnames(x)[[1]]) |>
    unlist() |>
    paste0("_", 1:length(sems))

  ## Remove name column and convert to matrix
  sems <- lapply(sems, \(x) as.matrix(x[,-1]))

  ## Get sequences for each variant ##
  offset <- max(unlist(lapply(sems, nrow)))

  ## Collect up/downstream sequences
  upstream <- getSeq(Hsapiens, chr, pos-offset, pos-1)
  dnstream <- getSeq(Hsapiens, chr, pos+1, pos+offset)

  ## Assemble risk and non-risk sequences
  nonRiskSeq <- xscat(upstream, ref, dnstream)
  riskSeq <- xscat(upstream, alt, dnstream)

  ## Initialize variables for holding sem scores
  nonRisk <- list()
  risk <- list()

  for(i in 1:length(sems)) {

    ## Score risk and non-risk sequences against each sem
    nonRiskScores <-
      PWMscoreStartingAt(pwm = t(sems[[i]]),
                         subject = nonRiskSeq,
                         starting.at = (offset-nrow(sems[[i]])+2):(offset+1))
    riskScores <-
      PWMscoreStartingAt(pwm = t(sems[[i]]),
                         subject = riskSeq,
                         starting.at = (offset-nrow(sems[[i]])+2):(offset+1))

    ## Find which frame has the best binding score (and record this)
    if (max(nonRiskScores) > max(riskScores)){
      frame <- which.max(nonRiskScores)
    } else {
      frame <- which.max(riskScores)
    }

    ## Select the score from the correct frame
    nonRisk[[i]] <- nonRiskScores[frame]
    risk[[i]] <- riskScores[frame]
  }

  ## Add names
  names(nonRisk) <- names(sems)
  names(risk) <- names(sems)

  ## Compile into data table
  dt <-
    data.table(motif = stack(nonRisk)$ind,
               nonrisk = do.call(rbind, nonRisk)[,1],
               risk = do.call(rbind, risk)[,1])

  ## Normalize to threshold
  # threshold <- 2^(unlist(lapply(sems, min)))
  bl <- read.table(semBaselines)
  threshold <- 2^bl[match(gsub(".sem", "", basename(semFiles)), bl[,1]),2]
  dt$nonrisk <- (2^dt$nonrisk - threshold) / abs(threshold)
  dt$risk <- (2^dt$risk - threshold) / abs(threshold)

  ## Return
  return(dt)

}
