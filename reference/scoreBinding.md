# Calculate binding propensity for all SEM motifs and genomic positions provided

Calculate binding propensity for all SEM motifs and genomic positions
provided

## Usage

``` r
scoreBinding(x, sem, genome, nFlank = NULL, seqId = NULL, rc = TRUE)
```

## Arguments

- x:

  `GRanges` object or a vector of DNA sequences

- sem:

  A `SNPEffectMatrix` or `SNPEffectMatrixCollection` object

- genome:

  A `BSgenome` object for the genome build to use. ie.
  [`BSgenome.Hsapiens.UCSC.hg19::Hsapiens`](https://rdrr.io/pkg/BSgenome.Hsapiens.UCSC.hg19/man/package.html).
  Required if providing a GRanges object. Ignored if providing a vector
  of sequences.

- nFlank:

  Number of flanking nucleotides to add to provided range. By default
  will add flank equal to the length of the longest motif. Ignored if
  providing a vector of sequences.

- seqId:

  Column in `GRanges` object to use for unique id. By default, ids will
  be generated from the `seqnames` and `ranges.` Ignored if not
  providing a `GRanges` object.

- rc:

  plot the reverse complement SEMs

## Value

If a `GRanges` object is provided, return a `SEMScores` object. If a
list of sequences is provided, just return the scoring table

## Examples

``` r
# load SEMs

# create a GRanges object
gr <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = 94136009
)

# calculate binding propensity
scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#> An object of class SEMScores
#> ranges(1): chr12:94136009
#> semData(12): transcription_factor, ensembl_id ... dnase_ENCODE_accession, PWM_source
#> scores(446):
#>               seqId    SEM     rc      score  scoreNorm index              seq
#>              <char> <char> <char>      <num>      <num> <int>           <char>
#>   1: chr12:94136009 TFAP2B    fwd  -1.689754 -0.3068238    15       GCTTTGAGGC
#>   2: chr12:94136009   ARNT    fwd  -6.892799 -0.9693833    17        TTTGAGGCA
#>   3: chr12:94136009   ATF1    fwd  -7.079925 -0.9420095    16      CTTTGAGGCAT
#>   4: chr12:94136009   ATF2    fwd  -4.890126 -0.9098440    16      CTTTGAGGCAT
#>   5: chr12:94136009   ATF3    fwd  -8.605675 -0.9885365    14      GGCTTTGAGGC
#>  ---                                                                          
#> 442: chr12:94136009 ZBTB7A    rev  -1.967170 -0.6612178    18        TTGAGGCAT
#> 443: chr12:94136009    ZFX    rev  -1.162039 -0.4693499    19       TGAGGCATCT
#> 444: chr12:94136009 ZNF281    rev  -5.713722 -0.9749858    14  GGCTTTGAGGCATCT
#> 445: chr12:94136009  ZNF18    rev  -6.739405 -0.9707101    14     GGCTTTGAGGCA
#> 446: chr12:94136009 ZSCAN4    rev -14.469707 -0.9996685    15 GCTTTGAGGCATCTGC
```
