# Access sems matrix from a SNPEffectMatrixCollection object

Access sems matrix from a SNPEffectMatrixCollection object

## Usage

``` r
getSEMs(x, semId = NULL)

# S4 method for class 'SNPEffectMatrixCollection'
getSEMs(x, semId = NULL)
```

## Arguments

- x:

  SNPEffectMatrixCollection object

- semId:

  optional `character` corresponding to an SEM in the
  SNPEffectMatrixCollection object. See semIds with `names(getSEMs(x))`.
  Defaults to returning all SEMs.

## Value

A position (rows) by nucleic acid (columns) `data.table` is returned

## Examples

``` r
## Create example SEM
df <- data.frame(
    A = c(1, 2, 3),
    C = c(1, 2, 3),
    G = c(1, 2, 3),
    T = c(1, 2, 3)
)
s <- SNPEffectMatrix(df, 1.205, "motif_id")
sc <- SNPEffectMatrixCollection(s)

## Access count matrix
getSEMs(sc)
#> $motif_id
#> An object of class SNPEffectMatrix
#> semId:  motif_id
#> baseline:  1.205
#> sem:
#>        A     C     G     T
#>    <num> <num> <num> <num>
#> 1:     1     1     1     1
#> 2:     2     2     2     2
#> 3:     3     3     3     3
#> 
```
