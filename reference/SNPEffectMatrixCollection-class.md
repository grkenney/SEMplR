# A collection of SNPEffectMatrix objects and their corresponding meta data

A collection of SNPEffectMatrix objects and their corresponding meta
data

## Slots

- `sems`:

  A list of SNPEffectMatrix objects

- `semData`:

  A `data.table` with meta data on each SEM

- `semKey`:

  The name of a column in semData that matches the semIds in the sems
  list as a `character`. If column entries have a .sem suffix, a new
  column named SEM_KEY will be created without the .sem suffixes.

## Examples

``` r
# build a SEM example matrix
m <- matrix(rnorm(16), nrow = 4)
colnames(m) <- c("A", "C", "G", "T")

# build a SNPEffectMatrix object
sm <- SNPEffectMatrix(sem = m, baseline = 1, semId = "sem_name")

# create a meta data table
md <- data.table::data.table(tf = "tf_name", sem = "sem_name")

# build a collection with 1 SEM
SNPEffectMatrixCollection(sems = list(sm), semData = md, semKey = "sem")
#> An object of class SNPEffectMatrixCollection
#> SEMs(1): sem_name
#> semData(2): tf, sem
```
