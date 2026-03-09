# SNPEffectMatrix object and constructor

Constructs a SNPEffectMatrix class object.

## Usage

``` r
SNPEffectMatrix(sem, baseline, semId)
```

## Arguments

- sem:

  A `data.table` object of format base position (rows) by nucleic acid
  (columns). Column names must inclue A, C, T, and G; other columns will
  be ignored.

- baseline:

  A `numeric` scrambled baseline, representing the binding score of
  randomly scrambled kmers of the same length.

- semId:

  A `character` unique identifier for the SEM.

## Value

a SNPEffectMatrix object

## Examples

``` r
# build a matrix for a motif of length 4
m <- matrix(rnorm(16), nrow = 4)
colnames(m) <- c("A", "C", "G", "T")

# build a SNPEffectMatrix object
SNPEffectMatrix(m, baseline = 1, semId = "sem_id")
#> An object of class SNPEffectMatrix
#> semId:  sem_id
#> baseline:  1
#> sem:
#>              A          C           G          T
#>          <num>      <num>       <num>      <num>
#> 1: -0.52201251  0.4681544  1.88850493 -0.8267890
#> 2: -0.05260191  0.3629513 -0.09744510 -1.5123997
#> 3:  0.54299634 -1.3045435 -0.93584735  0.9353632
#> 4: -0.91407483  0.7377763 -0.01595031  0.1764886
```
