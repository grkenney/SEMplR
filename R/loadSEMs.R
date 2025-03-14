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


#' Load .sem files and meta data into a SNPEffectMatrixCollection
#'
#' @param semFiles A list of paths to .sem files. Expects header of .sem files
#' to be in format `#BASELINE:{bl}` where `bl` is the numeric baseline value.
#' If matrix does not include baseline header, must be specified in `bl`.
#' @param semMetaData A `data.table` with meta data on each SEM
#' @param semMetaKey The name of a column in semData that matches the semIds in 
#' the sems list as a `character`. If column entries have a .sem suffix, a new 
#' column named SEM_KEY will be created without the .sem suffixes.
#' @param semIds Unique id for the sem as a `character` vector in same order
#' as sems. Defaults to semFile file name without the extension.
#' @param bls `numeric` vector or baseline values for the SEMs.
#' Overrides baseline specified in semFile header.
#'
#' @return A SNPEffectMatrix object
#' 
#' @export
loadSEMCollection <- \(semFiles, semMetaData=NULL, semMetaKey="", 
                       semIds=NULL, bls=NULL) {
  s <- lapply(1:length(semFiles), 
              \(i) loadSEM(semFile = semFiles[i],
                           semId = semIds[i],
                           bl = bls[i]))
  
  sc <- SNPEffectMatrixCollection(sems = s, 
                                  semData = semMetaData, 
                                  semKey = semMetaKey)
  return(sc)
}
