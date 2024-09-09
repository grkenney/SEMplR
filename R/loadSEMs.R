#' Define a function to calculate risk/non-risk binding propensity
#'
#' @param semFiles character vector of paths to each `.sem` file
#'
#' @export
loadSEMs <- \(semFiles) {

  ## Read in SEM files
  sems <- lapply(semFiles, data.table::fread, header = TRUE)

  ## Add names to sems
  names(sems) <-
    lapply(sems, \(x) colnames(x)[[1]]) |>
    unlist() |>
    paste0("_", 1:length(sems))

  ## Remove name column and convert to matrix
  sems <- lapply(sems, \(x) as.matrix(x[,-1]))

  return(sems)
}
