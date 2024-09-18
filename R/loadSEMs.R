#' Load .sem and baseline data from files or from github if file path is not
#' provided. Github default data: github.com/Boyle-Lab/SEMpl/raw/master/SEMs/
#'
#' @param sem_dir path to directory containing .sem files.
#' Defaults to pulling data from github.
#' @param baseline_file path to baselines file.
#' Defaults to pulling data from github.
#'
#' @export
loadSEMs <- \(sem_dir=NULL, baseline_file=NULL) {
  url_base <- "https://github.com/Boyle-Lab/SEMpl/raw/master/SEMs/"

  # if baseline file not provided, load data from github url,
  # otherwise, load from filepath provided
  if (is.null(baseline_file)) {
    baselines_url <- paste0(url_base, "BASELINE/SEMs_baseline_norm.txt")
    baselines <- utils::read.delim(url(baselines_url),
                                   sep = "\t", header = FALSE)

  } else if (!file.exists(baseline_file)) {
    stop("Path to baseline file does not exist: ", baseline_file)

  } else {
    baselines <- utils::read.delim(baseline_file, header = FALSE)
  }

  # if sem directory not provided, load data from github
  # otherwise, load data from .sem files in directory path provided
  if (is.null(sem_dir)) {
    sem_urls <- lapply(baselines[, 1],
                       function(name) {paste0(url_base, name, ".sem")}) |>
      unlist()
  } else {
    sem_files <- list.files(sem_dir, pattern = ".sem", full.names = TRUE)
  }

  sem_list <- list()
  for (i in 1:nrow(baselines)) {
    if (is.null(baseline_file)) {
      sem_matrix <- utils::read.delim(url(sem_urls[i]),
                                      sep = "\t")[-1, -1]
    } else {
      sem_matrix <- utils::read.delim(sem_files[i],
                                      sep = "\t")[-1, -1]
    }

    tfname <- as.character(baselines[i, 1])

    sem_list[[tfname]] <- SNPEffectMatrix(sem_matrix,
                               tf_name = tfname,
                               baseline = as.numeric(baselines[i, 2]),
                               pwm_filename = "")
  }

  return(sem_list)
}
