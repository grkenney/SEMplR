#' Load matrix and baseline data from .sem files
#'
#' @param semFile A path to a .sem file. Expects header of file to be in format
#' `#BASELINE:{bl}` where `bl` is the numeric baseline value. If matrix does
#' not include baseline header, must be specified in `bl`.
#' @param semId A unique id for the sem as a `character`. Defaults to
#' semFile file name without the extension.
#' @param bl `numeric` baseline value for the SEM. Overrides baseline specified
#' in semFile header.
#'
#' @return A SNPEffectMatrix object
#' 
#' @export
#' 
#' @examples
#' # write a tmp file to hold a SEM
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#' tf <- tempfile()
#' write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)
#' 
#' loadSEM(tf, semId = "sem_id", bl = 1)
#' 
loadSEM <- \(semFile, semId=NULL, bl=NULL) {
  s <- data.table::fread(file = semFile, sep = "\t")
  
  if (is.null(semId)) {
    semId <- gsub(pattern = ".sem", 
                  replacement = "", 
                  x = basename(semFile))
  }
  
  # if bl param is null and baseline is in header, use header baseline
  # if bl is null and baseline is not in header, stop
  # else, use bl param as baseline by default
  file_header <- readLines(semFile)[1]
  if (is.null(bl) & grepl("#BASELINE:", file_header)) {
    bl <- readLines(semFile)[1] |>
      strsplit(":") |>
      lapply("[[", 2) |>
      unlist()
  } else if (is.null(bl)) {
    stop("No baseline given. Baseline must be specified in semFile header or",
         " in bl parameter.")
  }
  
  return(SNPEffectMatrix(s, bl, semId))
}
