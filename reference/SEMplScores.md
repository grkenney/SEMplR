# SEMplScores object and constructor

Constructs a SEMplScores class object.

## Usage

``` r
SEMplScores(ranges = NULL, semData = NULL, scores = NULL)
```

## Arguments

- ranges:

  A `GRanges` or `VRanges` object to hold one or more variants

- semData:

  A named list of SNPEffectMatrix objects

- scores:

  (optional) A `data.table` object for motif information and binding
  scores

## Value

a SEMplScores object

## Examples

``` r
# load default SEMs
data(SEMC)

# create a VRanges object
vr <- VariantAnnotation::VRanges(
    seqnames = c("chr12", "chr19"),
    ranges = c(94136009, 10640062),
    ref = c("G", "T"), alt = c("C", "A")
)

SEMplScores(ranges = vr, semData = semData(SEMC))
#> An object of class SEMplScores
#> ranges(2): 
#> semData(12): transcription_factor, ensembl_id ... dnase_ENCODE_accession, PWM_source
#> scores(0):
```
