# Validate SEM format
.validateSEM <- function(sem, semFile) {
    ..expected_cols <- NA
    expected_cols <- c("A", "C", "G", "T")

    # check that all nucleotides are represented in columns of SEM
    if (!all(expected_cols %in% colnames(sem))) {
        rlang::abort(c(
            paste0(
                semFile, " does not contain all expected nucleotide columns"
            ),
            "i" = "SEM files must have 4 columns with names: 'A', 'C', 'G', 'T'"
        ))
    }

    # check that each nucleotide only has one column
    if (sum(colnames(sem) %in% expected_cols) > 4) {
        rlang::abort(c(
            paste0(
                semFile, " has more than one column for one or more nucleotides"
            ),
            "i" = paste0(
                "Each nucleotide ('A', 'C', 'G', 'T')",
                "must have only one column"
            )
        ))
    }

    return(sem[, ..expected_cols])
}


# Extract the baseline from file header if not explicitly stated
.extractBaseline <- function(semFile) {
    # if bl param is null and baseline is in header, use header baseline
    # if bl is null and baseline is not in header, stop
    # else, use bl param as baseline by default
    file_header <- readLines(semFile)[1]
    if (grepl("#BASELINE:", file_header)) {
        bl <- readLines(semFile)[1] |>
            strsplit(":") |>
            lapply("[[", 2) |>
            unlist()
    } else {
        rlang::abort(paste0(
            "No baseline given. Baseline must be specified in semFile header ",
            "or in bl parameter."
        ))
    }
    return(bl)
}

# Load matrix and baseline data from .sem files
# semFile: path to a .sem file. Expects header of file to be in format
# `#BASELINE:{bl}` where `bl` is the numeric baseline value. If matrix does
# not include baseline header, must be specified in `bl`.
# semId: A unique id for the sem as a `character`. Defaults to
# semFile: file name without the extension.
# bl: baseline value for the SEM. Overrides baseline specified
# in semFile header.
# delim: delimiter of the SEM file
.loadSEM <- function(semFile, semId = NULL, bl = NULL, delim = "\t") {
    s <- data.table::fread(file = semFile, sep = delim)

    s <- .validateSEM(s, semFile)

    # if no semId given, use the basename of the file
    if (is.null(semId)) {
        semId <- tools::file_path_sans_ext(basename(semFile))
    }

    # if baseline is not provided, extract from file
    if (is.null(bl)) {
        bl <- .extractBaseline(semFile)
    }

    return(SNPEffectMatrix(sem = s, baseline = bl, semId = semId))
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
#'
#' @examples
#' # write a tmp file to hold a SEM
#' m <- matrix(rnorm(16), nrow = 4)
#' colnames(m) <- c("A", "C", "G", "T")
#' tf <- tempfile()
#' write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)
#'
#' # build a meta data table
#' md <- data.table::data.table(
#'     transcription_factor = c("tf1"),
#'     cell_type = c("HepG2"),
#'     sem_id = c("sem_id")
#' )
#'
#' loadSEMCollection(tf,
#'     semMetaData = md, semIds = "sem_id",
#'     semMetaKey = "sem_id", bls = 1
#' )
#'
loadSEMCollection <- function(
  semFiles, semMetaData = NULL, semMetaKey = "",
  semIds = NULL, bls = NULL
) {
    s <- lapply(
        seq_along(semFiles),
        function(i) {
            .loadSEM(
                semFile = semFiles[i],
                semId = semIds[i],
                bl = bls[i]
            )
        }
    )

    if (!is.null(semMetaData) & semMetaKey == "") {
        rlang::abort(paste0(
            "If providing semMetaData, must specify a column name ",
            " linking the meta data to the semIds in semMetaKey"
        ))
    }

    sc <- SNPEffectMatrixCollection(
        sems = s,
        semData = semMetaData,
        semKey = semMetaKey
    )
    return(sc)
}
