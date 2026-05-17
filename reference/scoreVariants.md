# Calculate risk/non-risk binding propensity for all SEM motifs and variants provided

Calculate risk/non-risk binding propensity for all SEM motifs and
variants provided

## Usage

``` r
scoreVariants(
  x,
  sem,
  genome,
  refCol = NULL,
  altCol = NULL,
  varId = NULL,
  rc = TRUE
)
```

## Arguments

- x:

  `VRanges` object

- sem:

  a list of `SNPEffectMatrix` objects

- genome:

  A `BSgenome` object for the genome build to use. ie.
  [`BSgenome.Hsapiens.UCSC.hg19::Hsapiens`](https://rdrr.io/pkg/BSgenome.Hsapiens.UCSC.hg19/man/package.html)

- refCol:

  If providing a GRanges, the meta data column name with the reference
  (ref) allele. Ignored if providing a VRanges object.

- altCol:

  If providing a GRanges, the meta data column name with the alternative
  (alt) allele. Ignored if providing a VRanges object.

- varId:

  A column name in the meta data of x to use as a unique id.

- rc:

  plot the reverse complement SEMs

## Value

a SEMScores object with slots for the ranges of the provided variants
with an added `sequence` column with the sequence scored, the SEM
metadata, and the resulting scoring table. These slots are accessible
with the
[`getRanges()`](https://grkenney.github.io/SEMPLR/reference/getRanges.md),
[`semData()`](https://grkenney.github.io/SEMPLR/reference/semData.md),
and [`scores()`](https://grkenney.github.io/SEMPLR/reference/scores.md)
accessor functions respectively.

The scoring table will contain the following columns:

- varId: a unique identifier for the sequence scored

- SEM: the identifier for the SEM

- rc: the orientation of the sequence (fwd or rev)

- refSeq/altSeq: the motif sequence with the highest SEM scores within
  the sequence scored

- refScore/altScore: the unnormalized SEM score

- refNorm/altNorm: the SEM score normalized to it's corresponding
  baseline

- refVarIndex/altVarIndex: the index of the motif with the highest SEM
  score within the sequence scored

## Examples

``` r
library(VariantAnnotation)

# load default SEMs

# create a VRanges object
x <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C"
)

# calculate binding propensity
scoreVariants(x, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#> An object of class SEMScores
#> ranges(1): 
#> semData(12): transcription_factor, ensembl_id ... dnase_ENCODE_accession, PWM_source
#> scores(446):
#>                   varId            SEM     rc           refSeq           altSeq
#>                  <char>         <char> <char>           <char>           <char>
#>   1: chr12:94136009:G>C            AHR    fwd      TTTGAGGCATC      TTCAGGCATCT
#>   2: chr12:94136009:G>C            AHR    rev      TTTGAGGCATC      TTTCAGGCATC
#>   3: chr12:94136009:G>C AHR:ARNT:HIF1A    fwd        GGCTTTGAG        GGCTTTCAG
#>   4: chr12:94136009:G>C AHR:ARNT:HIF1A    rev        GGCTTTGAG        TTTCAGGCA
#>   5: chr12:94136009:G>C         ARID3A    fwd           TTTGAG           TTCAGG
#>  ---                                                                           
#> 442: chr12:94136009:G>C         ZNF217    rev         CTTTGAGG         CTTTCAGG
#> 443: chr12:94136009:G>C         ZNF281    fwd  GGAGAAGGCTTTGAG  AAGGAGAAGGCTTTC
#> 444: chr12:94136009:G>C         ZNF281    rev  GGCTTTGAGGCATCT  GGCTTTCAGGCATCT
#> 445: chr12:94136009:G>C         ZSCAN4    fwd GCTTTGAGGCATCTGC GCTTTCAGGCATCTGC
#> 446: chr12:94136009:G>C         ZSCAN4    rev GCTTTGAGGCATCTGC GCTTTCAGGCATCTGC
#>        refScore   altScore    refNorm    altNorm refVarIndex altVarIndex
#>           <num>      <num>      <num>      <num>       <int>       <int>
#>   1:  -1.495324  -1.393158 -0.4349553 -0.3934902          17          18
#>   2:  -1.354220  -1.158371 -0.3768978 -0.2862997          17          17
#>   3:  -1.304220  -1.308236 -0.3888573 -0.3905561          14          14
#>   4:  -1.266234  -1.086663 -0.3725524 -0.2893862          14          17
#>   5:  -1.400445  -1.428525 -0.4994304 -0.5090790          17          18
#>  ---                                                                    
#> 442:  -1.202679  -1.517704 -0.4025273 -0.5197292          16          16
#> 443:  -4.347612  -4.740351 -0.9355197 -0.9508865           8           6
#> 444:  -5.713722  -4.886938 -0.9749858 -0.9556316          14          14
#> 445: -15.439087 -13.002975 -0.9998307 -0.9990837          15          15
#> 446: -14.469707 -12.918336 -0.9996685 -0.9990284          15          15
```
