# Calculates binding enrichment for motif(s)

Perform a binomial test to determine if SNP Effect Matrices are bound
more often than expected.

## Usage

``` r
enrichSEMs(x, sem, background = NULL, seqs = NULL, nFlank = 0, genome = NULL)
```

## Arguments

- x:

  The scoring table produced by `scoreBinding`

- sem:

  A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object

- background:

  A `GRanges` object or a list of DNA sequences to use as a background
  set for the binomial test. The length of each sequence must match the
  length of sequences in `x`. By default, will scramble the provided
  sequences.

- seqs:

  The sequences scored in `scoreBinding`

- nFlank:

  Number of flanking nucleotides added to the sequences. Defaults to the
  length of the longest motif. If no flanks were added (ie, sequences
  were scored rather than a `GRanges`), use `nFlank = 0`.

- genome:

  A `BSGenome` object. Requried if background isn't specified.

## Value

a `list` of matrices

## Examples

``` r
# load SEMs

# note that this is a small example for demonstration purposes
# in actual enrichment analyses sets of 100+ ranges are recommended

# create a GRanges object
gr <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = 94136009
)

# calculate binding propensity
sb <- scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

enrichSEMs(sb, SEMC)
#> Building background set (this may take several minutes) ...
#>         SEM pvalue  padj n_bound n_bound_bg
#>      <char>  <num> <num>   <int>      <int>
#>   1: TFAP2B      1     1       0          0
#>   2:   ARNT      1     1       0          0
#>   3:   ATF1      1     1       0          0
#>   4:   ATF2      1     1       0          0
#>   5:   ATF3      1     1       0          0
#>  ---                                       
#> 219: ZBTB7A      1     1       0          0
#> 220:    ZFX      1     1       0          0
#> 221: ZNF281      1     1       0          0
#> 222:  ZNF18      1     1       0          0
#> 223: ZSCAN4      1     1       0          0
```
