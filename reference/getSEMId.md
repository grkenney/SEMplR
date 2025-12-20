# Access semId from a SNPEffectMatrix object

Access semId from a SNPEffectMatrix object

## Usage

``` r
getSEMId(x)

# S4 method for class 'SNPEffectMatrix'
getSEMId(x)
```

## Arguments

- x:

  SNPEffectMatrix object

## Value

The `character` id is returned

## Examples

``` r
# Isolate a single SNPEffectMatrix object from the default
# SNPEffectMatrixCollection
sm <- getSEMs(SEMC)[[1]]

# Access the SEM id
getSEMId(sm)
#> [1] "TFAP2B"
```
