# Default SNP Effect Matrix Data Collection

A collection of pre-computed SNP Effect Matrix objects to be used for
motif scoring

## Usage

``` r
SEMC
```

## Format

### `SEMC`

A SNPEffectMatrixCollection object containing 223 SEMs as
SNPEffectMatrix objects and a data frame with 223 rows and 13 columns
containing meta data:

- transcription_factor:

  Transcription factor name

- ensembl_id:

  Ensembl id

- ebi_complex_ac:

- uniprot_ac:

  Uniprot accession id

- PWM_id:

  Position weighted matrix id

- SEM:

  SNP Effect Matrix file

- SEM_baseline:

  SNP Effect Matrix baseline

- cell_type:

  Cell Type

- neg_log10_pval:

  -log10(p value) from SEMpl calculation

- chip_ENCODE_accession:

  ENCODE accession for ChIP data used in SEMpl

- dnase_ENCODE_accession:

  ENCODE accession for DNase data used in SEMpl

- PWM_source:

  Position weighted matrix source

## Source

<https://data.igvf.org/multireport/?type=ModelSet&software_versions.software.title=SEMpl>
