# SNPEffectMatrixCollection object and constructor

Constructs a SNPEffectMatrixCollection class object.

## Usage

``` r
SNPEffectMatrixCollection(sems, semData = NULL, semKey = "")
```

## Arguments

- sems:

  A list of SNPEffectMatrix objects

- semData:

  A `data.table` containing meta data for each SEM

- semKey:

  A column name in semData corresponding to the semIds in each
  SNPEffectMatrix object in the sems parameter

## Value

a SNPEffectMatrixCollection object

## Examples

``` r
#
# build two matries for a motif of length 4
m1 <- matrix(rnorm(16), nrow = 4)
m2 <- matrix(rnorm(16), nrow = 4)
colnames(m1) <- c("A", "C", "G", "T")
colnames(m2) <- c("A", "C", "G", "T")

# build two SNPEffectMatrix objects
s1 <- SNPEffectMatrix(m1, baseline = 1, semId = "sem_id_1")
s2 <- SNPEffectMatrix(m2, baseline = 1, semId = "sem_id_2")

# build a meta data table
md <- data.table::data.table(
    transcription_factor = c("tf1", "tf2"),
    cell_type = c("HepG2", "K562"),
    sem_id = c("sem_id_1", "sem_id_2")
)

# build a SNPEffectMatrixCollection object
SNPEffectMatrixCollection(list(s1, s2), semData = md, semKey = "sem_id")
#> Removing .sem suffixes from semKey. Formatted key now stored in column 'SEM_KEY'.
#> An object of class SNPEffectMatrixCollection
#> SEMs(2): sem_id_1, sem_id_2
#> semData(4): transcription_factor, cell_type, sem_id, SEM_KEY
```
