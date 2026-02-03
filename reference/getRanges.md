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
