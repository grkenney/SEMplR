# Run the Full Differential‐Expression Motif Enrichment Pipeline

`enrichmentSets()` is a one‐stop wrapper that takes a vector of
user‐provided IDs from a differential expression analysis and returns
promoter regions plus matched background sets for downstream motif
enrichment. It performs:

1.  ID mapping

2.  Optional biotype filtering

3.  Promoter coordinate extraction

4.  Background sampling within that same call

## Usage

``` r
enrichmentSets(
  txdb,
  orgdb,
  id_type,
  foreground_ids,
  background_ids = NULL,
  transcript = FALSE,
  threshold = 0.9,
  stripVersions = TRUE,
  inflateThresh = 1,
  geneType = NULL,
  overlapMinGap = 0,
  onePromoterPerGene = FALSE,
  n_ratio = 1,
  promoterWindow = c(upstream = 300, downstream = 50),
  standardChroms = TRUE,
  reduceOverlaps = TRUE
)
```

## Arguments

- txdb:

  A TxDb object. (e.g. `TxDb.Hsapiens.UCSC.hg38.knownGene`).

- orgdb:

  An OrgDb object. (e.g. `org.Hs.eg.db`).

- id_type:

  Type of identifier supplied in foreground and background IDs. See the
  `keytypes(orgdb)` for available id type options for each genome.

- foreground_ids:

  Character vector of gene or transcript IDs (e.g. Ensembl, RefSeq, gene
  symbols) to analyze.

- background_ids:

  Character vector of gene or transcript IDs to use as background set.

- transcript:

  Logical; `TRUE` to treat inputs as transcript‐level, `FALSE` for
  gene‐level.

- threshold:

  Numeric in range 0 to 1. Min fraction of IDs that must map to pick a
  keytype (default 0.9).

- stripVersions:

  Logical; strip version suffixes (e.g. ".1") from Ensembl/RefSeq IDs.

- inflateThresh:

  Numeric in range 0 to 1; max allowed transcript:gene inflation before
  auto‐collapsing (default 1).

- geneType:

  Optional character; biotype filter (e.g. `"protein-coding"`).

- overlapMinGap:

  Numeric; minimum gap when reducing overlapping promoters.

- onePromoterPerGene:

  Logical; if `TRUE`, pick only one promoter per gene.

- n_ratio:

  Numeric; ratio of ranges to retain in the background set. The number
  of background ranges will be equal to `n_ratio` multiplied by the
  number of foreground ranges.

- promoterWindow:

  Numeric named vector `c(upstream, downstream)`; promoter flank widths
  in bp (default `c(300,50)`).

- standardChroms:

  Logical; restrict to standard chromosomes.

- reduceOverlaps:

  Logical; merge overlapping promoter windows.

## Value

A `list` containing:

- backgroundElements:

  `GRanges` of sampled background promoters

- foregroundElements:

  `GRanges` of foreground promoters

- backgroundUniverse:

  `GRanges` of the pruned promoter pool

## Examples

``` r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#> Loading required package: GenomicFeatures
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: GenomicRanges
#> Loading required package: AnnotationDbi
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
library(org.Hs.eg.db)
#> 

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db

# Note that this is a minimal set for demonstration only
genes <- c("ENSG00000139618", "ENSG00000157764")
background <- c("ENSG00000111640", "ENSG00000075624")

# Minimal run with defaults:
results <- enrichmentSets(genes,
    background,
    txdb = txdb,
    orgdb = orgdb,
    id_type = "ENSEMBL"
)
#> Mapping foreground ids to ENTREZIDs...
#> 'select()' returned 1:1 mapping between keys and columns
#> Successfully mapped 100% of the provided foreground ids.
#> Mapping background ids to ENTREZIDs...
#> 'select()' returned 1:1 mapping between keys and columns
#> Successfully mapped 100% of the provided foreground ids.
#> Checking for inflation...
#> 'select()' returned 1:1 mapping between keys and columns
#> Combining overlapping promoter ranges within genes. (This step may take 1-2 minutes)...
#> Defining background elements...
```
