# Get sequence for genomic ranges and variants

Query sequences associated with a given range in a reference genome and
optionally include nucleotides up and downstream of the range of
interest sequences. Store the results in the metadata of the range
object.

## Usage

``` r
getRangeSeqs(x, genome, up = 0, down = 0, refCol = NULL, altCol = NULL)
```

## Arguments

- x:

  A `VRanges` or `GRanges` object with one or more variants. seqnames
  and ranges fields are required.

- genome:

  A `BSgenome` object for the genome build to use.

- up:

  Numeric, number of bases to return upstream of variant

- down:

  Numeric, number of bases to return downstream of variant

- refCol:

  Column name in meta data storing the ref allele

- altCol:

  Column name in meta data storing the alt allele

## Value

a `x` with added meta data columns. By default, only a "sequence" column
is added to the meta data. If variant information is supplied, either in
a `VRanges` object or through the "allele" parameter, "ref_seq" and
"alt_seq" columns are added to the meta data.

## Examples

``` r
# set reference genome
b <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

# Query sequnece for a GRanges object
gr <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = IRanges::IRanges(94136009, 94136019)
)
getRangeSeqs(gr, genome = b)
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames            ranges strand |    sequence
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]    chr12 94136009-94136019      * | GAGGCATCTGC
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# Query variant sequences for a VRanges object with 10bp flanks
vr <- VariantAnnotation::VRanges(
    seqnames = "chr12",
    ranges = IRanges::IRanges(94136009),
    ref = "G", alt = "C"
)
getRangeSeqs(vr, genome = b, up = 10, down = 10)
#> VRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand         ref              alt     totalDepth
#>          <Rle> <IRanges>  <Rle> <character> <characterOrRle> <integerOrRle>
#>   [1]    chr12  94136009      *           G                C           <NA>
#>             refDepth       altDepth   sampleNames softFilterMatrix |
#>       <integerOrRle> <integerOrRle> <factorOrRle>         <matrix> |
#>   [1]           <NA>           <NA>          <NA>                  |
#>                     ref_seq               alt_seq
#>                 <character>           <character>
#>   [1] AGAAGGCTTTGAGGCATCTGC AGAAGGCTTTCAGGCATCTGC
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>   hardFilters: NULL
```
