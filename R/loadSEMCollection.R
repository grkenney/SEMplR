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
#' 
#' @examples
#' # write a tmp file to hold a SEM
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#' tf <- tempfile()
#' write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)
#' 
#' # build a meta data table
#' md <- data.table::data.table(transcription_factor = c("tf1"), 
#'                              cell_type = c("HepG2"), 
#'                              sem_id = c("sem_id"))
#'                              
#' loadSEMCollection(tf, semMetaData = md, semIds = "sem_id", 
#'                   semMetaKey = "sem_id", bls = 1)
#' 
loadSEMCollection <- \(semFiles, semMetaData=NULL, semMetaKey="", 
                       semIds=NULL, bls=NULL) {
  s <- lapply(1:length(semFiles), 
              \(i) loadSEM(semFile = semFiles[i],
                           semId = semIds[i],
                           bl = bls[i]))
  
  if (!is.null(semMetaData) & is.null(semMetaKey)) {
    stop("If providing semMetaData, must specify a column name linking the meta
         data to the semIds in semMetaKey")
  }
  
  sc <- SNPEffectMatrixCollection(sems = s, 
                                  semData = semMetaData, 
                                  semKey = semMetaKey)
  return(sc)
}
