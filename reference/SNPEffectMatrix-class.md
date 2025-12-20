# SNP Effect Matrix (SEM)

A SNP Effect Matrix (SEM) estimates the binding affinity of every
possible mutation in a particular transcription factor (TF) binding
motif. Read more about SEMs here:
https://doi.org/10.1093/bioinformatics/btz612 This class contains three
slots: the matrix, the baseline value, and a unique id

## Slots

- `sem`:

  The SEM itself as a data table Rows represent sequence position
  (variable length), columns represent effects due to each nucleotide
  base A, C, G, T (fixed length: 4)

- `baseline`:

  A scrambled baseline, representing the binding score of randomly
  scrambled kmers of the same length. This is the binding cutoff for a
  TF

- `semId`:

  basename of the sem file

## Examples

``` r
# build a SEM example matrix
m <- matrix(rnorm(16), nrow = 4)
colnames(m) <- c("A", "C", "G", "T")

SNPEffectMatrix(sem = m, baseline = 1, semId = "sem_id")
#> An object of class SNPEffectMatrix
#> semId:  sem_id
#> baseline:  1
#> sem:
#>               A          C          G           T
#>           <num>      <num>      <num>       <num>
#> 1: -2.437263611 -1.8218177 -0.5536994  0.51242695
#> 2: -0.005571287 -0.2473253  0.6289820 -1.86301149
#> 3:  0.621552721 -0.2441996  2.0650249 -0.52201251
#> 4:  1.148411606 -0.2827054 -1.6309894 -0.05260191
```
