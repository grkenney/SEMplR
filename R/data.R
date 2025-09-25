#' Default SNP Effect Matrix Data Collection
#'
#' A collection of pre-computed SNP Effect Matrix objects to be used for motif
#' scoring
#'
#' @format ## `SEMC`
#' A SNPEffectMatrixCollection object containing 223 SEMs as SNPEffectMatrix
#' objects and a data frame with 223 rows and 13 columns containing meta data:
#' \describe{
#'   \item{transcription_factor}{Transcription factor name}
#'   \item{ensembl_id}{Ensembl id}
#'   \item{ebi_complex_ac}{}
#'   \item{uniprot_ac}{Uniprot accession id}
#'   \item{PWM_id}{Position weighted matrix id}
#'   \item{SEM}{SNP Effect Matrix file}
#'   \item{SEM_baseline}{SNP Effect Matrix baseline}
#'   \item{cell_type}{Cell Type}
#'   \item{neg_log10_pval}{-log10(p value) from SEMpl calculation}
#'   \item{chip_ENCODE_accession}{ENCODE accession for ChIP data used in SEMpl}
#'   \item{dnase_ENCODE_accession}{ENCODE accession for DNase data used
#'   in SEMpl}
#'   \item{PWM_source}{Position weighted matrix source}
#'   \item{SEM_KEY}{SEM id, used as data table key}
#'   ...
#' }
#' @source <>
"SEMC"
