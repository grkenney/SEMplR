# Add metadata to SNPEffectMatrix
# 
# @param x a SNPEffectMatrix object
# @param sem_metadata a data.frame or named vector with metadata for the 
# given sem.
# 
# @return SNPEffectMatrix object
addSemMetadata <- function(x, sem_metadata){
  x@tf <- sem_metadata$tf_name
  x@ensembl <- sem_metadata$ensembl_id
  x@uniprot <- sem_metadata$uniprot_id
  x@cellType <- sem_metadata$cell_type
  return(x)
}

#' Load .sem and baseline data from files or from github if file path is not
#' provided. Github default data: github.com/Boyle-Lab/SEMpl/raw/master/SEMs/
#'
#' @param semDir path to directory containing .sem files.
#' Defaults to pulling data from github.
#' @param baselineFile path to baselines file.
#' Defaults to pulling data from github.
#' @param metadataFile a csv file containing metadata for each sem with columns
#' sem_id, tf_name, ensembl_id, uniprot_id, cell_type where sem_id must match
#' the basenames of the sem files.
#'
#' @return named list of SNPEffectMatrix objects
#' 
#' @export
loadSEMs <- \(semDir=NULL, baselineFile=NULL, metadataFile=NULL) {
  url_base <- "https://github.com/Boyle-Lab/SEMpl/raw/master/SEMs/"
  
  if (!is.null(metadataFile)) {
    meta <- utils::read.delim(metadataFile, sep = ",")
  } else {
    meta <- NULL
  }

  # if baselineFile not provided, read from github
  if (is.null(baselineFile)) {
    baselines_url <- paste0(url_base, "BASELINE/SEMs_baseline_norm.txt")
    baselines <- utils::read.delim(url(baselines_url),
                                   sep = "\t", header = FALSE)
  } else {
    baselines <- utils::read.delim(baselineFile, header = FALSE)
  }

  # if semDir not provided, read from github
  if (is.null(semDir)) {
    sem_files <- lapply(baselines[, 1],
                       function(name) {paste0(url_base, name, ".sem")}) |>
      unlist()
  } else {
    sem_files <- list.files(semDir, pattern = ".sem", full.names = TRUE)
  }

  sem_list <- list()
  for (i in 1:length(sem_files)) {
    if (is.null(semDir)) {
      sem_matrix <- utils::read.delim(url(sem_files[i]),
                                      sep = "\t")[, -1]
    } else {
      sem_matrix <- utils::read.delim(sem_files[i],
                                      sep = "\t")[, -1]
    }

    sem_id <- tools::file_path_sans_ext(basename(sem_files[i]))
    baseline <- baselines[baselines[, 1] == sem_id, 2]

    sem_list[[sem_id]] <- SNPEffectMatrix(sem_matrix,
                               baseline = baseline,
                               semId = sem_id)
    
    if (!is.null(meta)) {
      sem_list[[sem_id]] <- addSemMetadata(x = sem_list[[sem_id]], 
                                           sem_metadata = meta[meta$sem_id == basename(sem_files[i]), ])
    }
  }
  
  return(sem_list)
}
