#' Permute sems for a given variant
#' 
#' @param sem_scores SemplScores object
#' @param motif_count matrix with a row for each motif and a column for each
#' motif change category (ie, gained, lost, changed, maintained)
#' @param variant_id unique id for variant to be permuted
#'
permuteVariant <- function(sem_scores, motif_count, variant_id) {
  # randomize sem ids
  sem_ids <- sem_scores@sem_metadata$sem_id
  randomized_sem_ids <- sample(sem_ids, replace = F)
  
  # subset scores to variant
  variant_scores <- variantScores(sem_scores, variant_id)
  
  # apply randomized labels
  variant_scores$sem <- randomized_sem_ids
  
  # count number of variants where each motif is gained/lost/maintained
  gained_motifs <- changed_motif(variant_scores, "gained")$sem
  motif_count[gained_motifs, "gained"] <- motif_count[gained_motifs, "gained"] + 1
  
  lost_motifs <- changed_motif(variant_scores, "lost")$sem
  motif_count[lost_motifs, "lost"] <- motif_count[lost_motifs, "lost"] + 1
  
  changed_motifs <- changed_motif(variant_scores, "changed")$sem
  motif_count[changed_motifs, "changed"] <- motif_count[changed_motifs, "changed"] + 1
  
  maintained_motifs <- changed_motif(variant_scores, "maintained")$sem
  motif_count[maintained_motifs, "maintained"] <- motif_count[maintained_motifs, "maintained"] + 1
  
  return(motif_count)
}

#' Permute sems for a given variant
#' 
#' @param permutation_counts matrix storing permutation counts
#' @param motif_count matrix with a row for each motif and a column for each
#' motif change category (ie, gained, lost, changed, maintained)
#' @param sem unique sem id
#'
countAlteredMotifs <- function(permutation_counts, motif_count, sem) {
  if (motif_count[sem, "gained"] >= permutation_counts[sem, "actual_gained"]) {
    permutation_counts[sem, "perm_gt_gained"] <- permutation_counts[sem, "perm_gt_gained"] + 1
  }
  if (motif_count[sem, "lost"] >= permutation_counts[sem, "actual_lost"]) {
    permutation_counts[sem, "perm_gt_lost"] <- permutation_counts[sem, "perm_gt_lost"] + 1
  }
  if (motif_count[sem, "changed"] >= permutation_counts[sem, "actual_changed"]) {
    permutation_counts[sem, "perm_gt_changed"] <- permutation_counts[sem, "perm_gt_changed"] + 1
  }
  if (motif_count[sem, "maintained"] >= permutation_counts[sem, "actual_maintained"]) {
    permutation_counts[sem, "perm_gt_maintained"] <- permutation_counts[sem, "perm_gt_maintained"] + 1
  }
  return(permutation_counts)
}

#' Permute sems for a given variant
#' 
#' @param permutation_counts matrix storing permutation counts
#' @param num_permutations number of permutations to perform
#'
calculatePvalues <- function(permutation_counts, num_permutations) {
  # calculate pvalue
  permutation_counts[, "pval_gained"] <- unlist(lapply(permutation_counts[, "perm_gt_gained"], 
                                                       function(x) x / num_permutations))
  permutation_counts[, "pval_lost"] <- unlist(lapply(permutation_counts[, "perm_gt_lost"], 
                                                     function(x) x / num_permutations))
  permutation_counts[, "pval_changed"] <- unlist(lapply(permutation_counts[, "perm_gt_changed"], 
                                                     function(x) x / num_permutations))
  permutation_counts[, "pval_maintained"] <- unlist(lapply(permutation_counts[, "perm_gt_maintained"], 
                                                           function(x) x / num_permutations))
  
  # calculated BH adjusted pvalue
  permutation_counts[, "adjpval_gained"] <- unlist(lapply(permutation_counts[, "pval_gained"],
                                                          function(x) p.adjust(x, method="BH")))
  permutation_counts[, "adjpval_lost"] <- unlist(lapply(permutation_counts[, "pval_lost"],
                                                        function(x) p.adjust(x, method="BH")))
  permutation_counts[, "adjpval_changed"] <- unlist(lapply(permutation_counts[, "pval_changed"],
                                                        function(x) p.adjust(x, method="BH")))
  permutation_counts[, "adjpval_maintained"] <- unlist(lapply(permutation_counts[, "pval_maintained"],
                                                              function(x) p.adjust(x, method="BH")))
  return(permutation_counts)
}

#' Count the number of variants with altered motifs for each motif
#' 
#' @param sem_scores SemplScores object
#' @param pc matrix storing permutation counts
#' 
#' @return pc matrix with actual counts populated
tallyAlteredMotifs <- function(sem_scores, pc) {
  for (s in rownames(pc)) {
    m <- motifScores(sem_scores, s)
    pc[s, "actual_gained"] <- nrow(changed_motif(m, "gained"))
    pc[s, "actual_lost"] <- nrow(changed_motif(m, "lost"))
    pc[s, "actual_changed"] <- nrow(changed_motif(m, "changed"))
    pc[s, "actual_maintained"]  <- nrow(changed_motif(m, "maintained"))
  }
  return(pc)
}

#' Perform a permuation test to find changed motifs among all variants analyzed
#' 
#' @param sem_scores SemplScores object that has already been scored with 
#' semMotifBinding
#' @param n number of permutations to run
#' 
#' @return a data.frame with numbers of variants with altered motifs and 
#' adjusted pvalues
#' 
#' @export
semEnrichment <- function(sem_scores, n) {
  permutation_counts <- matrix(0, nrow=nrow(sem_scores@sem_metadata), ncol=16)
  colnames(permutation_counts) <- c("actual_gained", "perm_gt_gained", 
                                    "pval_gained", "adjpval_gained",
                                    "actual_lost", "perm_gt_lost", 
                                    "pval_lost", "adjpval_lost",
                                    "actual_changed", "perm_gt_changed",
                                    "pval_changed", "adjpval_changed",
                                    "actual_maintained", "perm_gt_maintained", 
                                    "pval_maintained", "adjpval_maintained")
  rownames(permutation_counts) <- sem_scores@sem_metadata$sem_id
  
  permutation_counts <- tallyAlteredMotifs(sem_scores, permutation_counts)
  
  motif_count_init <- matrix(0, nrow=nrow(sem_scores@sem_metadata), ncol=4)
  rownames(motif_count_init) <- sem_scores@sem_metadata$sem_id
  colnames(motif_count_init) <- c("gained", "lost", "changed", "maintained")
  
  num_permutations <- n
  for (j in 1:num_permutations){
    # initialize motif counter to store number of variants where each motif is
    # gained, lost, or maintained
    motif_count <- motif_count_init
    for (i in 1:length(sem_scores@variants)){
      variant_id <- sem_scores@variants$id[i]
      motif_count <- permuteVariant(sem_scores, motif_count, variant_id)
    }
    
    # count if number of variants where motif is altered exceeds actual
    for (sem in rownames(permutation_counts)) {
      permutation_counts <- countAlteredMotifs(permutation_counts, motif_count, sem)
    }
  }
  
  permutation_counts <- calculatePvalues(permutation_counts = permutation_counts, 
                                         num_permutations = num_permutations)
  result <- permutation_counts[, c("actual_changed", "adjpval_changed",
                                   "actual_gained", "adjpval_gained",
                                   "actual_lost", "adjpval_lost",
                                   "actual_maintained", "adjpval_maintained")] |>
    as.data.frame()
  
  colnames(result) <- c("changed", "adjpval_changed",
                        "gained", "adjpval_gained",
                        "lost", "adjpval_lost",
                        "maintained", "adjpval_maintained")
  result <- result[order(result$adjpval_changed, decreasing = FALSE), ]
  
  return(result)
}

