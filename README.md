The microRMS R package
================
Lars Snipen

Installation
============

Start R and run

``` r
devtools::github_install("larssnip/microRMS")
```

Several of the R functions calls upon the software VSEARCH that must be installed and available on the system, see [VSEARCH on GitHub](https://github.com/torognes/vsearch).

Tutorial - composition
======================

This is a short step-by-step tutorial on a small toy example to illustrate how the various functions would be used in a typical RMS study where we want to estimate the abundance of various taxa in a microbial community.

The RMS database
----------------

First, we need a collection of sequenced genomes that cover the community we want to study. Just as for shotgun metagenome data, this is a 'closed reference' type of taxonomic assignment. The use fo RMS to discover new taxa is out of the scope of this tutorial.

From every genome, we collect all RMS fragments, cluster these, and create the *fragment cluster copy number matrix* which is the central data structure of this method.

### Collecting fragments

First, we collect the RMS fragments from each genome, and store these as fasta-formatted files. Genomes are stored as fasta files, and each genome results in a new fasta file with RMS fragments.

``` r
library(tidyverse)
library(microseq)
# genome.files contain the (full path) names of the existing genome fasta files
# fragment.files contain the (full path) names of the fragment files to be created here
for(i in 1:length(genome.files)){
  getRMSfragments(genome.files[i]) %>% 
    writeFasta(out.file = fragment.files[i])
}
```

This job is done once for each genome, and stored. Thus, if you later want to add more genomes, you only do this for the new genomes.

### Making the database

This means we collect all fragments from the genomes. First, we write all genome-wise fragment files into one temporary fasta file:

``` r
for(i in 1:length(fragment.files)){
  ok <- file.append(file.path(tmp.fdir, "all.frg"),
                    genome.files[i], overwrite = T)
}
```

Set `tmp.dir` to some folder of temporary results that you typically delete in the end.

The creation of the database object requires the use of the `vsearch` software for clustering. The function `RMSdbase()` will not work unless `vsearch` is a valid command in your system/computer.

``` r
rmsdb.obj <- RMSdbase(file.path(tmp.fdir, "all.frg"),
                      identity = 0.99,
                      min.length = 30,
                      max.length = 500)
```

Here we use only default settings for the options `identity`, `min.length`and `max.length`.

The resulting object `rmsdb.obj` is a `list` with two tables and a mtarix, the latter is the fragment cluster copy number matrix.

Processing reads
----------------

Estimation
----------
