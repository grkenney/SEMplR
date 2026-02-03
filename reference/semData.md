# Accessor semData slot in a SEMScores object

Accessor semData slot in a SEMScores object

Access semData from a SNPEffectMatrixCollection object

## Usage

``` r
semData(x)

# S4 method for class 'SEMScores'
semData(x)

# S4 method for class 'SNPEffectMatrixCollection'
semData(x)
```

## Arguments

- x:

  SNPEffectMatrixCollection object

## Value

A `data.table` is returned

## Examples

``` r
library(VariantAnnotation)

# load default SEMs

# create a VRanges object
vr <- VRanges(
    seqnames = "chr12",
    ranges = 94136009,
    ref = "G", alt = "C"
)

# calculate binding propensity
s <- scoreVariants(vr, SEMC, BSgenome.Hsapiens.UCSC.hg19::Hsapiens)

semData(s)
#> Key: <transcription_factor>
#>      transcription_factor                                      ensembl_id
#>                    <char>                                          <char>
#>   1:                  AHR                                 ENSG00000106546
#>   2:       AHR:ARNT:HIF1A ENSG00000106546;ENSG00000143437;ENSG00000100644
#>   3:               ARID3A                                 ENSG00000116017
#>   4:                 ARNT                                 ENSG00000143437
#>   5:                ARNTL                                 ENSG00000133794
#>  ---                                                                     
#> 219:                  ZFX                                 ENSG00000005889
#> 220:                ZNF18                                 ENSG00000154957
#> 221:               ZNF217                                 ENSG00000171940
#> 222:               ZNF281                                 ENSG00000162702
#> 223:               ZSCAN4                                 ENSG00000180532
#>      ebi_complex_ac           uniprot_ac                PWM_id
#>              <char>               <char>                <char>
#>   1:                              P35869                M00778
#>   2:                P35869;P27540;Q16665                M00976
#>   3:                              Q99856              MA0151.1
#>   4:                              P27540  ARNT_HUMAN.H11MO.0.B
#>   5:                              O00327 BMAL1_HUMAN.H11MO.0.A
#>  ---                                                          
#> 219:                              P17010   ZFX_HUMAN.H11MO.0.A
#> 220:                              P17022 ZNF18_HUMAN.H11MO.0.C
#> 221:                              O75362                M01306
#> 222:                              Q9Y2X9 ZN281_HUMAN.H11MO.0.A
#> 223:                              Q8NAM6  Zscan4_pwm_secondary
#>                         SEM SEM_baseline cell_type neg_log10_pval
#>                      <char>        <num>    <char>          <num>
#>   1:             M00778.sem   -0.6717610     HepG2      18.350950
#>   2:             M00976.sem   -0.5938010     HepG2      14.290190
#>   3:           MA0151.1.sem   -0.4020880     HepG2       2.980258
#>   4: ARNT_HUMAN.GM12878.sem   -1.8632604   GM12878      12.081160
#>   5:  BMAL1_HUMAN.HepG2.sem   -1.5977194     HepG2      20.220840
#>  ---                                                             
#> 219:    ZFX_HUMAN.HepG2.sem   -0.2478717     HepG2       3.443155
#> 220: ZNF18_HUMAN.HEK293.sem   -1.6459529    HEK293       3.793244
#> 221:             M01306.sem   -0.4596240   GM12878       2.843803
#> 222:  ZN281_HUMAN.HepG2.sem   -0.3926150     HepG2       3.428189
#> 223:   ZSCAN4_secondary.sem   -2.9110180    HEK293      11.973690
#>      chip_ENCODE_accession dnase_ENCODE_accession  PWM_source
#>                     <char>                 <char>      <char>
#>   1:           ENCFF242PUG            ENCFF001UVU    TRANSFAC
#>   2:           ENCFF242PUG            ENCFF001UVU    TRANSFAC
#>   3:           ENCFF460WFF            ENCFF001UVU      JASPAR
#>   4:           ENCFF538DNI            ENCFF759OLD HOCOMOCOv11
#>   5:           ENCFF363NHH            ENCFF897NME HOCOMOCOv11
#>  ---                                                         
#> 219:           ENCFF329DCV            ENCFF897NME HOCOMOCOv11
#> 220:           ENCFF888BHA            ENCFF285OXK HOCOMOCOv11
#> 221:           ENCFF807UJN            ENCFF097LEF    TRANSFAC
#> 222:           ENCFF948PYK            ENCFF897NME HOCOMOCOv11
#> 223:           ENCFF485WHZ            ENCFF910QHN    UNIPROBE

## Create example SEM
df <- data.frame(
    A = c(1, 2, 3),
    C = c(1, 2, 3),
    G = c(1, 2, 3),
    T = c(1, 2, 3)
)
sem_data <- data.frame(SEM = "motif_id", TF = "tf_id")
s <- SNPEffectMatrix(df, 1.205, "motif_id")
sc <- SNPEffectMatrixCollection(s, semData = sem_data, semKey = "SEM")

## Access count matrix
semData(sc)
#> Key: <SEM>
#>         SEM     TF
#>      <char> <char>
#> 1: motif_id  tf_id
```
