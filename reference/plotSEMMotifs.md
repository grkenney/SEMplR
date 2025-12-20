# Plot non-alt versus alt binding propensity for a single variant

Plot non-alt versus alt binding propensity for a single variant

## Usage

``` r
plotSEMMotifs(
  s,
  variant,
  label = "transcription_factor",
  labsize = 4,
  cols = c("#F8766D", "dodgerblue2"),
  ptsize = 1
)
```

## Arguments

- s:

  a SEMplScores object with scores populated

- variant:

  variant id to plot

- label:

  column in sem_metadata slot of semplObj to use for point labels

- labsize:

  numeric size of the point labels

- cols:

  vector of length 2 with colors to use for plotting gained and lost
  motifs respectively

- ptsize:

  numeric size of the ggplot points

## Value

a `ggplot` scatter plot of non-alt binding propensity versus alt binding
propensity

## Examples

``` r
library(VariantAnnotation)
data(SEMC)

# create a VRanges object
vr <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C"
)

# calculate binding propensity
s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

plotSEMMotifs(s, "chr12:94136009:G>C", label = "transcription_factor")

```
