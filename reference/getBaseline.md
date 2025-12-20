# Access baseline from a SNPEffectMatrix object

Access baseline from a SNPEffectMatrix object

## Usage

``` r
getBaseline(x)

# S4 method for class 'SNPEffectMatrix'
getBaseline(x)
```

## Arguments

- x:

  SNPEffectMatrix object

## Value

The `numeric` baseline value is returned

## Examples

``` r
# Isolate a single SNPEffectMatrix object from the default
# SNPEffectMatrixCollection
sm <- getSEMs(SEMC)[[1]]

# Access the baseline
getBaseline(sm)
#> [1] -1.161048
```
