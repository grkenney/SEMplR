# Map User IDs to Entrez and Determine Best Keytype

Given a vector of user-supplied gene/transcript IDs, finds the
AnnotationDbi keytype (e.g. ens, refseq, symbol, etc.) that maps the
highest fraction of inputs, and returns both foreground and full
background sets as Entrez IDs. Optionally collapses transcript-style
inputs to genes when requested or when reverse-mapping inflation exceeds
a threshold.

## Usage

``` r
mapIDs(
  orgdb,
  id_type,
  foreground_ids,
  background_ids = NULL,
  threshold = 0.9,
  transcript = FALSE,
  stripVersions = TRUE,
  inflateThresh = 1
)
```

## Arguments

- orgdb:

  An OrgDb object. (e.g. `org.Hs.eg.db`).

- id_type:

  Type of identifier supplied in foreground and background IDs.

- foreground_ids:

  Character vector of gene or transcript IDs (e.g. Ensembl, RefSeq, gene
  symbols) to analyze.

- background_ids:

  Character vector of gene or transcript IDs to use as background set.

- threshold:

  Fraction in range 0 to 1; minimum mapping rate to accept a keytype
  without falling back (default 0.9).

- transcript:

  Logical; if `TRUE`, analyze as transcript-level IDs (default `FALSE`).

- stripVersions:

  Logical; strip version suffixes (e.g. ".1") from Ensembl/RefSeq IDs.

- inflateThresh:

  Fraction in range 0 to 1; if reverse-mapping shows excessive,
  inflation automatically collapse transcripts to genes (default 1 ie.
  100%).

## Value

A named list combining the original `mapping` components with:

- `fg_ids`:

  data.frame(entrez, mappedID) for your foreground set

- `bg_ids`:

  data.frame(entrez, mappedID) for the full background

- `userIDtype`:

  the chosen keytype (e.g. "ensembl")

- `transcript`:

  logical, whether transcript-level mapping was used

## Examples

``` r
library(org.Hs.eg.db)
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

# Transcript Ids
my_transcripts <- c("ENST00000245479", "ENST00000633194")
ids <- mapIDs(
    orgdb = orgdb,
    foreground_ids = my_transcripts,
    id_type = "ENSEMBLTRANS",
    transcript = TRUE
)
#> Mapping foreground ids to ENTREZIDs...
#> 'select()' returned 1:1 mapping between keys and columns
#> Successfully mapped 100% of the provided foreground ids.
#> Building background id set...
#> 'select()' returned 1:many mapping between keys and columns
```
