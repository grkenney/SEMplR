# Calculate binding propensity for all SEM motifs and genomic positions provided

Calculate binding propensity for all SEM motifs and genomic positions
provided

## Usage

``` r
scoreBinding(x, sem, genome, nFlank = NULL, seqId = NULL)
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

## Value

If a `GRanges` object is provided, return a `SEMplScores` object. If a
list of sequences is provided, just return the scoring table

## Examples

``` r
# load SEMs
data(SEMC)

# create a GRanges object
gr <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = 94136009
)

# calculate binding propensity
scoreBinding(gr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
#> An object of class SEMplScores
#> ranges(1): chr12:94136009
#> semData(12): transcription_factor, ensembl_id ... dnase_ENCODE_accession, PWM_source
#> scores(223):
#>               seqId    SEM      score  scoreNorm index              seq
#>              <char> <char>      <num>      <num> <int>           <char>
#>   1: chr12:94136009 TFAP2B  -1.689754 -0.3068238    15       GCTTTGAGGC
#>   2: chr12:94136009   ARNT  -6.892799 -0.9693833    17        TTTGAGGCA
#>   3: chr12:94136009   ATF1  -7.079925 -0.9420095    16      CTTTGAGGCAT
#>   4: chr12:94136009   ATF2  -4.890126 -0.9098440    16      CTTTGAGGCAT
#>   5: chr12:94136009   ATF3  -8.605675 -0.9885365    14      GGCTTTGAGGC
#>  ---                                                                   
#> 219: chr12:94136009 ZBTB7A  -1.859506 -0.6349682    12        AAGGCTTTG
#> 220: chr12:94136009    ZFX  -1.459472 -0.5682106    19       TGAGGCATCT
#> 221: chr12:94136009 ZNF281  -4.347612 -0.9355197     8  GGAGAAGGCTTTGAG
#> 222: chr12:94136009  ZNF18  -5.410220 -0.9264060    15     GCTTTGAGGCAT
#> 223: chr12:94136009 ZSCAN4 -15.439087 -0.9998307    15 GCTTTGAGGCATCTGC
```
