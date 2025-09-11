# SEMplR <a href="https://grkenney.github.io/SEMplR"><img src="man/figures/SEMplR-B.png" align="right" height="200" alt="SEMplR website" style="float:right; height:200px;" /></a>


SEMplR (SNP Effect Matrix Pipeline in R) is an R package that predicts transcription factor (TF) binding. SEMplR can be used to predict binding affinity of TFs at genomic loci or predict the affect of genetic variation on TF binding.

This package extends the functionality of the [SEMpl](https://github.com/Boyle-Lab/SEMpl) (SNP Effect Matrix) command line tool developed by the Boyle Lab at the University of Michigan. If you use SEMplR in your work, please also cite SEMpl:

> Sierra S Nishizaki, Natalie Ng, Shengcheng Dong, Robert S Porter, Cody Morterud, Colten Williams, Courtney Asman, Jessica A Switzenberg, Alan P Boyle, Predicting the effects of SNPs on transcription factor binding affinity, Bioinformatics, Volume 36, Issue 2, 15 January 2020, Pages 364–372, https://doi.org/10.1093/bioinformatics/btz612


## Installation

```
devtools::install_github("grkenney/SEMplR")
```

## Basic Usage

Below are some examples of basic usage. Please see the [vignettes](https://grkenney.github.io/SEMplR/) for more detailed workflow examples.


### Predicting transcription factor binding

SEMplR accepts GRanges objects or lists of sequences to score. Given the two loci below, the `scoreBinding` produces a table with 446 rows (an entry for each loci and SEM combination).

```
gr <- GenomicRanges::GRanges(seqnames = c("chr12", "chr19"),
                             ranges = c(94136009, 10640062))
result <- scoreBinding(gr, 
                       semList = sc, 
                       bs_genome_obj = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
```

### Predicting effect of genetic variation on transcription factor binding

SEMplR accepts both VRanges and GRanges objects, specifying a reference an alternative allele. Every variant is scored against every SEM and a scoring is done for each allele independently.

The resulting object contains three slots containing the variants scored, SEM meta data, and the scoring table. These can be accessed with the `variants()`, `semData()`, and `scores()` functions respectively.

```
vr <- VRanges(seqnames = c("chr12", "chr19"),
              ranges = c(94136009, 10640062), 
              ref = c("G", "T"), alt = c("C", "A"))
sempl_obj <- scoreVariants(vr = vr,
                           semList = sc,
                           bs_genome_obj=BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
```



