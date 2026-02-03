# Filter Foreground and Background ID Sets by Gene Type

`poolFilter()` takes the mapped foreground and background ID data frames
(as produced by
[`mapIDs()`](https://grkenney.github.io/SEMPLR/reference/mapIDs.md))
and, if requested, filters both sets to only include genes (or their
transcripts) of a specified biotype (e.g. “protein-coding”).

## Usage

``` r
poolFilter(mapped, geneType = NULL)
```

## Arguments

- mapped:

  A list returned by
  [`mapIDs()`](https://grkenney.github.io/SEMPLR/reference/mapIDs.md),
  containing at least:

  - `fg_ids`: data.frame with columns `entrez` and `mappedID`
    (foreground).

  - `bg_ids`: data.frame with columns `entrez` and `mappedID`
    (background).

  - `so_obj`: a `src_organism` object for transcript lookups.

  - `orgdb`: the loaded OrgDb package object.

  - `transcript`: logical, whether IDs are transcripts.

- geneType:

  Optional character scalar; if not NULL, only genes of this biotype
  (`GENETYPE` in the OrgDb) will be retained. Valid values vary by
  organism (e.g. “protein-coding”, “lncRNA”, etc.).

## Value

The original `mapped` list, but with `fg_ids` and `bg_ids` replaced by
filtered versions (only rows matching `geneType`, if provided).

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
```
