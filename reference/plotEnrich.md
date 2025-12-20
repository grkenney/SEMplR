# Plot the results of `enrichSEMs`

Generates a circular dendrogram, clustering SNP Effect Matrices on
similarity and a heatmap representing the -log10 transformed adjusted
p-value of a SEMplR enrichment.

## Usage

``` r
plotEnrich(
  e,
  sem,
  label = "transcription_factor",
  method = "WPCC",
  threshold = 0.05,
  textCols = c("darkgrey", "black"),
  textCex = 0.7,
  heatmapCols = c("white", "red"),
  pvalRange = c(0, 20)
)
```

## Arguments

- e:

  The resulting data.table from `enrichSEMs`

- sem:

  A `SNPEffectMatrixCollection` object

- label:

  Column in semData(sem) to use for tree labels

- method:

  Method to use for SEM comparison. See ?universalmotif::compare_motifs
  for options.

- threshold:

  The adjusted p-value threshold for coloring SEMs

- textCols:

  A vector of two colors to label non-significant and significant SEMs
  respectively.

- textCex:

  Text size of SEM labels.

- heatmapCols:

  A vector of two colors to use for the heatmap, ordered low to high
  `-log10(padj)`.

- pvalRange:

  A vector of 2 numerics to use as the scale range for the heatmap of
  `-log10(padj)`.

## Value

a `ggtree` object

## Examples

``` r
# load SEMs
data(SEMC)

# note that this is a small example for demonstration purposes
# in actual enrichment analyses sets of 100+ ranges are recommended

# create a GRanges object
gr <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = 94136009
)

# calculate binding propensity
sb <- scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

e <- enrichSEMs(sb, SEMC)
#> Building background set (this may take several minutes) ...
plotEnrich(e, SEMC)

```
