# validate user supplied txdb
.validateTxDb <- \(txdb) {
  if (!is(txdb, "TxDb")) {
    rlang::abort(c("Invalid txdb.",
                   "x" = paste0("txdb is class '", is(txdb)[1],
                                "'. An object of 'TxDb' class is expected.")))
  }
}


# validate user supplied txdb
.validateOrgDb <- \(orgdb) {
  if (!is(orgdb, "OrgDb")) {
    rlang::abort(c("Invalid orgdb",
                   "x" = paste0("orgdb is class '", is(orgdb)[1],
                                "'. An object of 'OrgDb' class is expected.")))
  }
}


# given a TxDb object, extract the genome build
.extractTxDbBuild <- \(txdb) {
  build <- txdb$packageName |>
    sub(pattern = "^TxDb\\.[^.]+\\.UCSC\\.", replacement = "") |>
    sub(pattern = "\\..*$", replacement = "")
  return(build)
}


#' Build an organism-specific id-mapping object
#'
#' @description
#' \code{buildMappingObject()} installs (if needed) and loads the appropriate
#' \pkg{OrgDb} and \pkg{TxDb} Bioconductor packages for a given species and
#' genome build, then returns an \code{src_organism} object along with metadata
#' about which packages and parameters were used.
#'
#' @param organism Character scalar specifying the species, e.g.
#' \dQuote{Homo sapiens}.
#' @param genomeBuild Character scalar giving the genome build
#' (e.g. \dQuote{hg38}),
#'   or \dQuote{auto} (default) to pick the latest supported build.
#' @param txdb Character scalar naming a specific \pkg{TxDb} package, or
#'   \dQuote{auto} (default) to select the most recent \pkg{TxDb} for the build.
#'
#' @return A named list with components:
#' \describe{
#'   \item{so_obj}{An \code{src_organism} data source
#'   \item{orgdb}{The loaded \pkg{OrgDb} package namespace object.}
#'   \item{organism}{The species string used.}
#'   \item{genomeBuild}{The resolved genome build string.}
#'   \item{txdb}{The resolved \pkg{TxDb} package name.}
#' }
#'
#' @examples
#' mapping <- buildMappingObject(
#'   organism    = "Homo sapiens",
#'   genomeBuild = "auto",
#'   txdb        = "auto"
#' )
#'
#' @export
buildMappingObject <- \( txdb,
                         orgdb ) {

  # validate input dbs
  .validateTxDb(txdb)
  .validateOrgDb(orgdb)

  build <- .extractTxDbBuild(txdb)

  # Return both the data‐source objects and the parameters that created them
  mapping_dbs <- list(
    genomeBuild = build,
    orgdb       = orgdb,
    txdb        = txdb
  )
  return(mapping_dbs)
}
