# Accessor scores slot in a SEMplScores object

Accessor scores slot in a SEMplScores object

## Usage

``` r
scores(x)

# S4 method for class 'SEMplScores'
scores(x)
```

## Arguments

- x:

  a SEMplScores object

## Examples

``` r
library(VariantAnnotation)

# load default SEMs
data(SEMC)

# create a VRanges object
vr <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C"
)

# calculate binding propensity
s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

scores(s)
#>                   varId          semId           refSeq           altSeq
#>                  <char>         <char>           <char>           <char>
#>   1: chr12:94136009:G>C            AHR      TTTGAGGCATC      TTCAGGCATCT
#>   2: chr12:94136009:G>C AHR:ARNT:HIF1A        GGCTTTGAG        GGCTTTCAG
#>   3: chr12:94136009:G>C         ARID3A           TTTGAG           TTCAGG
#>   4: chr12:94136009:G>C           ARNT        TTTGAGGCA        TTTCAGGCA
#>   5: chr12:94136009:G>C          ARNTL      GAGGCATCTGC      TTCAGGCATCT
#>  ---                                                                    
#> 219: chr12:94136009:G>C            ZFX       TGAGGCATCT       TCAGGCATCT
#> 220: chr12:94136009:G>C          ZNF18     GCTTTGAGGCAT     GGCTTTCAGGCA
#> 221: chr12:94136009:G>C         ZNF217         AGGCTTTG         AGGCTTTC
#> 222: chr12:94136009:G>C         ZNF281  GGAGAAGGCTTTGAG  AAGGAGAAGGCTTTC
#> 223: chr12:94136009:G>C         ZSCAN4 GCTTTGAGGCATCTGC GCTTTCAGGCATCTGC
#>        refScore   altScore    refNorm    altNorm refVarIndex altVarIndex
#>           <num>      <num>      <num>      <num>       <int>       <int>
#>   1:  -1.495324  -1.393158 -0.4349553 -0.3934902          17          18
#>   2:  -1.304220  -1.308236 -0.3888573 -0.3905561          14          14
#>   3:  -1.400445  -1.428525 -0.4994304 -0.5090790          17          18
#>   4:  -6.892799  -4.947151 -0.9693833 -0.8820613          17          17
#>   5:  -5.810051  -4.861625 -0.9460537 -0.8958962          20          18
#>  ---                                                                    
#> 219:  -1.459472  -1.478959 -0.5682106 -0.5740037          19          19
#> 220:  -5.410220  -7.107575 -0.9264060 -0.9773072          15          14
#> 221:  -1.661981  -1.523076 -0.5654352 -0.5215141          13          13
#> 222:  -4.347612  -4.740351 -0.9355197 -0.9508865           8           6
#> 223: -15.439087 -13.002975 -0.9998307 -0.9990837          15          15
```
