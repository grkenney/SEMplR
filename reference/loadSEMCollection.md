# Load .sem files and meta data into a SNPEffectMatrixCollection

Load .sem files and meta data into a SNPEffectMatrixCollection

## Usage

``` r
loadSEMCollection(
  semFiles,
  semMetaData = NULL,
  semMetaKey = "",
  semIds = NULL,
  bls = NULL
)
```

## Arguments

- semFiles:

  A list of paths to .sem files. Expects header of .sem files to be in
  format `#BASELINE:{bl}` where `bl` is the numeric baseline value. If
  matrix does not include baseline header, must be specified in `bl`.

- semMetaData:

  A `data.table` with meta data on each SEM

- semMetaKey:

  The name of a column in semData that matches the semIds in the sems
  list as a `character`. If column entries have a .sem suffix, a new
  column named SEM_KEY will be created without the .sem suffixes.

- semIds:

  Unique id for the sem as a `character` vector in same order as sems.
  Defaults to semFile file name without the extension.

- bls:

  `numeric` vector or baseline values for the SEMs. Overrides baseline
  specified in semFile header.

## Value

A SNPEffectMatrix object

## Examples

``` r
# write a tmp file to hold a SEM
m <- matrix(rnorm(16), nrow = 4)
colnames(m) <- c("A", "C", "G", "T")
tf <- tempfile()
write.table(m, tf, quote = FALSE, sep = "\t", row.names = FALSE)

# build a meta data table
md <- data.table::data.table(
    transcription_factor = c("tf1"),
    cell_type = c("HepG2"),
    sem_id = c("sem_id")
)

loadSEMCollection(tf,
    semMetaData = md, semIds = "sem_id",
    semMetaKey = "sem_id", bls = 1
)
#> An object of class SNPEffectMatrixCollection
#> SEMs(1): sem_id
#> semData(3): transcription_factor, cell_type, sem_id
```
