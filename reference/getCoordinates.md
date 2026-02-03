# Retrieve Promoter Regions and Sample Background Elements

Given a mapped set of foreground IDs (with Entrez and mappedID), this
function:

1.  Fetches genomic coordinates for each feature (gene or transcript),
    according to the chosen `TSS.method`.

2.  Extracts promoter windows around those coordinates.

3.  Optionally reduces overlaps and selects one promoter per gene.

4.  Samples background promoter elements matching the foreground

## Usage

``` r
getCoordinates(
  mapped,
  txdb,
  transcript = FALSE,
  n_ratio = 1,
  promoterWindow = c(upstream = 300, downstream = 50),
  standardChroms = TRUE,
  reduceOverlaps = TRUE,
  overlapMinGap = 0,
  onePromoterPerGene = FALSE
)
```

## Arguments

- mapped:

  A list as returned by
  [`mapIDs()`](https://grkenney.github.io/SEMPLR/reference/mapIDs.md),
  containing at least `fg_ids` and `bg_ids`.

- txdb:

  A TxDb object. (e.g. `TxDb.Hsapiens.UCSC.hg38.knownGene`).

- transcript:

  Logical; `TRUE` for transcript-level coordinates, `FALSE` for
  gene-level.

- n_ratio:

  Numeric; ratio of ranges to retain in the background set. The number
  of background ranges will be equal to `n_ratio` multiplied by the
  number of foreground ranges.

- promoterWindow:

  Numeric named vector of lengths: `c(upstream, downstream)` (default
  `c(300,50)`).

- standardChroms:

  Logical; restrict to standard chromosomes.

- reduceOverlaps:

  Logical; merge any overlapping promoter windows.

- overlapMinGap:

  Numeric; minimum gap when reducing overlaps.

- onePromoterPerGene:

  Logical; if `TRUE`, choose one promoter per gene.

## Value

A named `list` from the final `.defineBackgroundElements()`:

- backgroundElements:

  `GRanges` of sampled background promoters

- foregroundElements:

  `GRanges` of foreground promoters

- backgroundUniverse:

  `GRanges` of the pruned universe

- matchObject:

  `MatchedGRanges` if `bgMethod="matched"`, else `NULL`

## Examples

``` r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

my_genes <- c("ENSG00000139618", "ENSG00000157764")
ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_genes,
    id_type = "ENSEMBL"
)
#> Mapping foreground ids to ENTREZIDs...
#> 'select()' returned 1:1 mapping between keys and columns
#> Successfully mapped 100% of the provided foreground ids.
#> Building background id set...
#> 'select()' returned 1:many mapping between keys and columns
#> Checking for inflation...
#> 'select()' returned 1:1 mapping between keys and columns
filtered <- poolFilter(ids, geneType = "protein-coding")
coords <- getCoordinates(mapped = filtered, txdb = txdb)
#> Combining overlapping promoter ranges within genes. (This step may take 1-2 minutes)...
#> Defining background elements...
```
