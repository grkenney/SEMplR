# Access ranges slot in a SEMScores object

Access ranges slot in a SEMScores object

## Usage

``` r
getRanges(x)

# S4 method for class 'SEMScores'
getRanges(x)
```

## Arguments

- x:

  a SEMScores object

## Value

A GRanges or VRanges object

## Examples

``` r
library(VariantAnnotation)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘matrixStats’
#> The following objects are masked from ‘package:Biobase’:
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> The following object is masked from ‘package:Biobase’:
#> 
#>     rowMedians
#> Loading required package: SummarizedExperiment
#> Loading required package: Rsamtools
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
#> 
#> Attaching package: ‘VariantAnnotation’
#> The following object is masked from ‘package:base’:
#> 
#>     tabulate

# load default SEMs

# create a VRanges object
vr <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C"
)

# calculate binding propensity
s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

getRanges(s)
#> VRanges object with 1 range and 3 metadata columns:
#>       seqnames    ranges strand         ref              alt     totalDepth
#>          <Rle> <IRanges>  <Rle> <character> <characterOrRle> <integerOrRle>
#>   [1]    chr12  94136009      *           G                C           <NA>
#>             refDepth       altDepth   sampleNames softFilterMatrix |
#>       <integerOrRle> <integerOrRle> <factorOrRle>         <matrix> |
#>   [1]           <NA>           <NA>          <NA>                  |
#>                      ref_seq                alt_seq                 id
#>                  <character>            <character>        <character>
#>   [1] CCGTCAAGGAGAAGGCTTTG.. CCGTCAAGGAGAAGGCTTTC.. chr12:94136009:G>C
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>   hardFilters: NULL
```
