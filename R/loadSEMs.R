#' Add metadata to SNPEffectMatrix
#' 
#' @param sem a SNPEffectMatrix object
#' @param metadata a data.frame or named vector with metadata for the given sem.
#' 
#' @return SNPEffectMatrix object
addSemMetadata <- function(sem, metadata){
  sem@tf_name <- metadata$tf_name
  sem@ensembl_id <- metadata$ensembl_id
  sem@uniprot_id <- metadata$uniprot_id
  sem@cell_type <- metadata$cell_type
  return(sem)
}

#' Load .sem and baseline data from files or from github if file path is not
#' provided. Github default data: github.com/Boyle-Lab/SEMpl/raw/master/SEMs/
#'
#' @param sem_dir path to directory containing .sem files.
#' Defaults to pulling data from github.
#' @param baseline_file path to baselines file.
#' Defaults to pulling data from github.
#' @param metadata_file a csv file containing metadata for each sem with columns
#' sem_id, tf_name, ensembl_id, uniprot_id, cell_type where sem_id must match
#' the basenames of the sem files.
#'
#' @return named list of SNPEffectMatrix objects
#' 
#' @export
loadSEMs <- \(sem_dir=NULL, baseline_file=NULL, metadata_file=NULL) {
  url_base <- "https://github.com/Boyle-Lab/SEMpl/raw/master/SEMs/"
  
  if (!is.null(metadata_file)) {
    metadata <- utils::read.delim(metadata_file, sep = ",")
  } else {
    metadata <- NULL
  }

  # if baseline_file not provided, read from github
  if (is.null(baseline_file)) {
    baselines_url <- paste0(url_base, "BASELINE/SEMs_baseline_norm.txt")
    baselines <- utils::read.delim(url(baselines_url),
                                   sep = "\t", header = FALSE)
  } else {
    baselines <- utils::read.delim(baseline_file, header = FALSE)
  }

  # if sem_dir not provided, read from github
  if (is.null(sem_dir)) {
    sem_urls <- lapply(baselines[, 1],
                       function(name) {paste0(url_base, name, ".sem")}) |>
      unlist()
  } else {
    sem_files <- list.files(sem_dir, pattern = ".sem", full.names = TRUE)
  }

  sem_list <- list()
  for (i in 1:nrow(baselines)) {
    if (is.null(sem_dir)) {
      sem_matrix <- utils::read.delim(url(sem_urls[i]),
                                      sep = "\t")[-1, -1]
    } else {
      sem_matrix <- utils::read.delim(sem_files[i],
                                      sep = "\t")[-1, -1]
    }

    sem_id <- as.character(baselines[i, 1])

    sem_list[[sem_id]] <- SNPEffectMatrix(sem_matrix,
                               baseline = as.numeric(baselines[i, 2]),
                               sem_id = sem_id)
    
    if (!is.null(metadata)) {
      sem_list[[sem_id]] <- addSemMetadata(sem_list[[sem_id]], 
                     metadata[metadata$sem_id == basename(sem_files[i]), ])
    }
  }
  
  return(sem_list)
}
