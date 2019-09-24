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

Tutorial - microbial community composition
==========================================

This is a short step-by-step tutorial on a small toy example to illustrate a typical RMS study where we want to estimate the abundance of various taxa in a microbial community.

Download the archive [RMStutorial.zip](), and unzip it to some folder. Start R/RStudio and make the `RMStutorial` your working directory for this R session. Create a new R-script, copy the code chunks below into it, and save it in the RMStutorial folder. Step through this, and inspect the results.

Note that after some code chunks there is some output, usually starting with `##`. You do not copy this into your R-script.

Several of the function in this R package calls upon the software `vsearch`, see <https://github.com/torognes/vsearch>.

Creating the RMS object
=======================

First, we need a collection of sequenced genomes that cover the community we want to study. Just as for shotgun metagenome data, this is a 'closed reference' type of taxonomic assignment. The use of RMS to discover new taxa is out of the scope of this tutorial.

For this tutorial, we have collected 10 genomes in the `gnm/` folder. From every genome, we collect all RMS fragments, cluster these, and create the fragment cluster *copy number matrix* which is the central data structure of this method.

All information is stored in an *rms object*, which is simply a `list` with several tables and matrices. It is convenient to have all data structures assembled into one object like this, and since this is a simple list you have full access to all components.

Collecting RMS fragments
------------------------

First, we collect the RMS fragments from each genome, and store these as fasta-formatted files in a separate folder. You also find the file `gnm/genome_table.txt` among the genomes, containing a small table with one row for each genome. It is required you have such a metadata table with information about the genomes. This table *must contain* the columns `genome_id` and `genome_file`, but may contain any other genome metadata columns in addition to these.

The `genome_id` should be a text without spaces, and unique to each genome. It is added to the header-lines of the fragment fasta files, to indicate the genome of origin for all fragments. In this example the species name (note underscores, not spaces) has been used. The filenames in the `genome_file` column must all be unique, and we could have used them, or a prefix of them, as well.

The `genome_file` column specify the genome fasta files, but without the path to where they are located. *Keep it this way*. The reason is we re-use these filenames for the fragment fasta files, but in another folder. Thus, we supply the path as a separate argument when needed.

The first code chunk reads genomes, collect their RMS fragments, and store these as fasta files in the `frg/` folder:

``` r
library(tidyverse)
library(microrms)
gnm.dir <- "gnm"   # genome fasta files are found here
frg.dir <- "frg"   # fragment fasta files ends up here
genome.tbl <- suppressMessages(read_delim("gnm/genome_table.txt", delim = "\t"))
for(i in 1:nrow(genome.tbl)){
  readFasta(file.path(gnm.dir, genome.tbl$genome_file[i])) %>% 
    getRMSfragments(genome.id = genome.tbl$genome_id[i]) %>% 
    writeFasta(out.file = file.path(frg.dir, genome.tbl$genome_file[i]))
}
```

    ## Genome Bacteroides_uniformis : found 1511 RMS-fragments
    ## Genome Bacteroides_vulgatus : found 2033 RMS-fragments
    ## Genome Parabacteroides_distasonis : found 1621 RMS-fragments
    ## Genome Bacteroides_dorei : found 2265 RMS-fragments
    ## Genome Bacteroides_thetaiotaomicron : found 1959 RMS-fragments
    ## Genome [Eubacterium]_rectale : found 794 RMS-fragments
    ## Genome Roseburia_inulinivorans : found 1283 RMS-fragments
    ## Genome Bacteroides_xylanisolvens : found 2366 RMS-fragments
    ## Genome Bacteroides_ovatus : found 2781 RMS-fragments
    ## Genome Blautia_obeum : found 1893 RMS-fragments

We observe that all these genome shave plenty of RMS fragments, given the default restriction enzymes used here (see `?getRMSfragments`).

This job is done once for each genome, and the fragment files are stored. Thus, if you later want to add more genomes, you only run this for the new genomes. Notice that the fragment files have names identical to the genome files, and must be in a separate folder. Thus, the `genome_file` column in `genome.tbl` is used to name both genome and fragment files.

The RMS object
--------------

Select the genomes you want to be able to recognize later. This means you may `slice()` or `filter()` your `genome.tbl` to only contain the rows (=genomes) you are interested in. Here we use all genomes.

The creation of the RMS object requires the use of the `vsearch` software for clustering. The function `RMSobject()` will not work unless `vsearch` is a valid command in your system/computer.

``` r
rms.obj <- RMSobject(genome.tbl, frg.dir)
```

    ## VSEARCH clustering...
    ## ...produced 17648 clusters
    ## ...the cluster table...done
    ## ...the copy number matrix
    ## genome 1 / 10 
    genome 2 / 10 
    genome 3 / 10 
    genome 4 / 10 
    genome 5 / 10 
    genome 6 / 10 
    genome 7 / 10 
    genome 8 / 10 
    genome 9 / 10 
    genome 10 / 10 
    ## ...the genome table...done

We notice the fragments clustered into more than 17000 clusters, using the default `identity`, see `?RMSobject` for details.

The resulting object `rms.obj` is a `list` with two tables and a matrix. The `Genome.tbl` is a copy of the `genome.tbl`, but has got two new columns, and should be inspected right away. The column `N_clusters` lists the number of fragment clusters in each genome, and `N_unique` how many of these are unique to each genome:

``` r
print(rms.obj$Genome.tbl)
```

    ## # A tibble: 10 x 6
    ##    genome_id     genome_file   tax_id organism_name     N_clusters N_unique
    ##    <chr>         <chr>          <dbl> <chr>                  <int>    <int>
    ##  1 Bacteroides_~ GCA_00347154~    820 Bacteroides unif~       2250     1945
    ##  2 Bacteroides_~ GCA_00347513~    821 Bacteroides vulg~       1943     1895
    ##  3 Parabacteroi~ GCA_00346362~    823 Parabacteroides ~       2016     1720
    ##  4 Bacteroides_~ GCA_00346646~ 357276 Bacteroides dore~       2346     2070
    ##  5 Bacteroides_~ GCA_00346936~    818 Bacteroides thet~       1870     1847
    ##  6 [Eubacterium~ GCA_00347477~  39491 [Eubacterium] re~       1600     1576
    ##  7 Roseburia_in~ GCA_00347003~ 360807 Roseburia inulin~        778      754
    ##  8 Bacteroides_~ GCA_00346444~ 371601 Bacteroides xyla~       2754     2480
    ##  9 Bacteroides_~ GCA_00347064~  28116 Bacteroides ovat~       1504     1479
    ## 10 Blautia_obeum GCA_00346464~  40520 Blautia obeum AF~       1269     1249

In this case there are plenty of unique clusters for all genomes. If two or more genomes are very closely related, this number will shrink towards zero, making the recognition of each impossible.

Dendrogram based on Jaccard distances
-------------------------------------

You get a more detailed picture of the relation between genomes by computing their Jaccard distance from the shared/unique clusters, and then plot it as a dendrogram:

``` r
library(ggdendro)
d <- dist(t(rms.obj$Cpn.mat), method = "binary")
tree <- hclust(d, method = "average")
ggd <- ggdendrogram(dendro_data(tree),
                    rotate = T,
                    theme_dendro = F) +
  labs(x = "", y = "Jaccard distance")
print(ggd)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

All branches look deep and nice here, as expected. If branches get too shallow (Jaccard distance close to zero) in this tree, the genomes in that clade will be diffuclt/impossible to separate since they share too many clusters. The solution to this is to prune the genome collection, keeping only one representative of such a tight clade in the database, as a representative of the entire clade.

Processing reads
================

The read processing means essentially taking the fastq files from sequencing as input and producing a fasta file as output, for each sample.

The data from sequencing are the paired fastq files in the folder `fastq/`. Again, it is required to have a table like the one in `fastq/sample_table.txt` with metadata about each sample. There should always be a column `sample_id` with a unique text for each sample. Also, the `sample.tbl` must have two columns `R1_file` and `R2_file` specifying the corresponding fastq filenames.

In addition we also create the required column `reads_file` below, containing the name of the resulting fasta files with reads for each sample. Below we create this from the `sample_id`, which means this text must be useful as a filename prefix (e.g. no `/` or spaces inside). The table may contain this column already, with no need to create it. Again, paths *should not* be part of any filenames, we supply them as separate inputs.

There is no R function for doing the read processing, since this may be done in many different ways. Here is an R script with some suggested code for doing this processing using the `vsearch` software. The (long) output from this is hidden in this document:

``` r
#################
### The settings
fq.dir <- "fastq"   # path to folder with (input) fastq files
fa.dir <- "fasta"   # path to folder with (output) fasta files
tmp.dir <- "tmp"
PCR.forward.primer <- "GACTGCGTACCAATTC"
PCR.reverse.primer <- "GATGAGTCCTGAGTAA"
min.read.length <- 30
maxee <- 0.02

#####################
### The sample table
suppressMessages(read_delim("fastq/sample_table.txt", delim = "\t")) %>% 
  mutate(reads_file = str_c(sample_id, ".fasta")) -> sample.tbl

############################################################
### Looping over all samples
### 1) Filtering by maxee, discarding read-pairs
### 2) Merging read-pairs
### 3) Trimming primers from merged reads
### 4) Trimming primers from un-merged reads
### 5) Writing all reads to fasta-file
### 6) De-replicating and saving one fasta-file per sample
###########################################################
Nf <- str_length(PCR.forward.primer)
Nr <- str_length(PCR.reverse.primer)
for(i in 1:nrow(sample.tbl)){
  cat("\n\n##### VSEARCH quality filtering sample", sample.tbl$sample_id[i], "...\n")
  cmd <- paste("vsearch",
               "--fastq_filter", file.path(fq.dir, sample.tbl$R1_file[i]),
               "--reverse",      file.path(fq.dir, sample.tbl$R2_file[i]),
               "--fastq_maxee_rate", maxee,
               "--fastqout", file.path(tmp.dir, "filtered_R1.fq"),
               "--fastqout_rev", file.path(tmp.dir, "filtered_R2.fq"))
  system(cmd)
  
  cat("\n\n##### VSEARCH mergings read-pairs...\n")
  cmd <- paste("vsearch",
               "--fastq_mergepairs", file.path(tmp.dir, "filtered_R1.fq"),
               "--reverse",          file.path(tmp.dir, "filtered_R2.fq"),
               "--fastq_allowmergestagger",
               "--fastq_minmergelen", min.read.length,
               "--fastaout", file.path(tmp.dir, "merged.fa"),
               "--fastqout_notmerged_fwd", file.path(tmp.dir, "notmerged_R1.fq"),
               "--fastqout_notmerged_rev", file.path(tmp.dir, "notmerged_R2.fq"))
  system(cmd)
  
  cat("\n\n##### VSEARCH trimming primers from merged reads...\n")
  cmd <- paste("vsearch",
               "--fastx_filter", file.path(tmp.dir, "merged.fa"),
               "--fastq_stripleft",  Nf,
               "--fastq_stripright", Nr,
               "--fastq_minlen", min.read.length,
               "--relabel", "'size=2;pair'",
               "--fastaout", file.path(tmp.dir, "merged_filt.fa"))
  system(cmd)
  
  cat("\n\n##### VSEARCH trimming primers from un-merged reads...\n")
  cmd <- paste("vsearch",
               "--fastq_filter", file.path(tmp.dir, "notmerged_R1.fq"),
               "--fastq_stripleft", Nf,
               "--fastq_minlen", min.read.length,
               "--relabel", str_c("'size=1;notmerged_R1_'"),
               "--fastaout", file.path(tmp.dir, "notmerged_R1_filt.fa"))
  system(cmd)
  cmd <- paste("vsearch",
               "--fastq_filter", file.path(tmp.dir, "notmerged_R2.fq"),
               "--fastq_stripleft", Nr,
               "--fastq_minlen", min.read.length,
               "--relabel", str_c("'size=1;notmerged_R2_'"),
               "--fastaout", file.path(tmp.dir, "notmerged_R2_filt.fa"))
  system(cmd)
  
  cat("\n\n##### VSEARCH adding all reads to one fasta-file...\n")
  ok <- file.append(file.path(tmp.dir, "merged_filt.fa"),
                    file.path(tmp.dir, "notmerged_R1_filt.fa"))
  cmd <- paste("vsearch",
               "--fastx_revcomp", file.path(tmp.dir, "notmerged_R2_filt.fa"),
               "--fastaout", file.path(tmp.dir, "notmerged_R2_filt_rc.fa"))
  system(cmd)
  ok <- file.append(file.path(tmp.dir, "merged_filt.fa"),
                    file.path(tmp.dir, "notmerged_R2_filt_rc.fa"))
 
  cat("\n\n##### VSEARCH de-replicating sample", sample.tbl$sample_id[i], "...\n")
  cmd <- paste("vsearch",
               "--derep_fulllength", file.path(tmp.dir, "merged_filt.fa"),
               "--minuniquesize", 1,
               "--minseqlength", min.read.length,
               "--sizein --sizeout",
               "--relabel", str_c(sample.tbl$sample_id[i], ":uread_"),
               "--output", file.path(fa.dir, sample.tbl$reads_file[i]))
  system(cmd)
}
```

Note that you should have created the folder `fasta` and `tmp` before running this script. The first is where the resulting fasta files appear. The second is just temporary files. They may be nice to have for debugging, in case something goes wrong, but should be deleted in the end when everything runs smoothly.

Before we are done with this step, we add the `sample.tbl` to our `rms.obj`. Note that `sample.tbl` must have at least the two columns `sample_id` and `reads_file` for the downstream analysis (the `R1_file` and `R2_file` may still be present but are no longer needed):

``` r
rms.obj <- addSampleTable(rms.obj, sample.tbl)
```

In this way we have the information about our samples in the same object as we have all other information.

Mapping reads to clusters
=========================

The next step is to map reads from each sample to the fragment clusters, and obtain a *readcount matrix*. This matrix has one column for each sample, and one row for each fragment cluster, not unlike an OTU or ASV matrix for 16S amplicon data.

The readcount matrix
--------------------

The `readMapper()` function needs the `rms.obj` with information about fragment clusters (`rms.obj$Cluster.tbl`) and samples (`rms.obj$Sample.tbl`), and uses `vsearch` to search with reads against the fragment cluster centroids. The path to the fasta files with processed reads is also required, since the `rms.obj$Sample.tbl` has the filenames, but not their path.

``` r
rms.obj <- readMapper(rms.obj, fa.dir)
```

The matrix `Readcount.mat` is added as a new component to the returned `rms.obj`.

Length normalization
--------------------

We must expect the readcounts from RMS amplicons to have some length bias, due to the PCR amplification. We may plot and see if this is indeed the case:

``` r
rms.obj$Cluster.tbl %>% 
  select(Length) %>% 
  bind_cols(as_tibble(rms.obj$Readcount.mat)) %>% 
  gather(key = "Sample", value = "Readcounts", -Length) %>%  
  ggplot() +
  geom_point(aes(x = Length, y = Readcounts), alpha = 0.3) +
  scale_y_log10() +
  facet_wrap(~Sample) -> plt
print(plt)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

Since we have limited fragments to lengths of 30 to 500 only, the bias is not severe. Still, we may try to normalize:

``` r
rms.obj.norm <- normLength(rms.obj)
```

We make the same plot with the normalized data:

``` r
rms.obj.norm$Cluster.tbl %>% 
  select(Length) %>% 
  bind_cols(as_tibble(rms.obj.norm$Readcount.mat)) %>% 
  gather(key = "Sample", value = "Readcounts", -Length) -> tbl
plt %+% tbl %>% print
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

There is some effect on the most extreme lengths, as usual. Note the log-transformed y-axis here.

If we decide to stick to the normalized data, we do not need both `rms.obj` and `rms.obj.norm`:

``` r
rms.obj <- rms.obj.norm
rm(rms.obj.norm)
```

Estimating abundances
=====================

Finally, we estimate the abundances of the various genomes in our samples.

The Constrained Least Square estimation
---------------------------------------

Abundance estimation is done by the `rmscols()` function. The supplied `rms.obj` must contain the cluster copy number matrix `rms.obj$Cpn.mat` and the readcounts `rms.obj$Readcount.mat`. The idea is to look for a linear combination of genome abundances that, given cluster copy numbers, best explains the observed readcounts in a sample. The output has been hidden from the following code:

``` r
relative.abundances <- rmscols(rms.obj)
```

The `relative.abundances` is a matrix with one row for each genome, and one column for each sample. The numbers in a column are the relative contributions of each genome to this sample.

We may plot the results as stacked bar charts:

``` r
relative.abundances %>% 
  as_tibble(rownames = "Genome") %>% 
  gather(key = "Sample", value = "Estimated", -Genome) -> long.tbl
p1 <- ggplot(long.tbl) +
  geom_col(aes(x = Sample, y = Estimated, fill = Genome), color = "black")
print(p1)
```

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

Comparing to gold standard
--------------------------

The sequencing fastq files in `fastq/` are simulated data. The `art` software (<https://www.ncbi.nlm.nih.gov/pubmed/22199392>) was used to simulate Illumina HiSeq reads. Biases typical for RMS data (i.e. due to fragment length and fragment-specific amplification bias) were added as described in [Snipen et al, 2019]().

The file `gold_standard.txt` contains the actual relative abundances of all genomes in all samples. Let us compare the estimated abundances from above to this:

``` r
suppressMessages(read_delim("gold_standard.txt", delim = "\t")) %>%
  rename(Genome = genome_id) %>% 
  gather(key = "Sample", value = "Gold.standard", -Genome) %>% 
  full_join(long.tbl, by = c("Genome", "Sample")) %>% 
  gather(key = "Type", value = "Abundance", -Sample, -Genome) %>% 
  ggplot() +
  geom_col(aes(x = Type, y = Abundance, fill = Genome), color = "black") +
  facet_wrap(~Sample)
```

![](README_files/figure-markdown_github/unnamed-chunk-15-1.png)

New data
========

This forms a code template, and by replacing the genomes in `gnm/` by your genomes, and the sequencing data in `fastq/` by your own data, it should be possible to run an analysis.

Beware that for large collections of genomes (thousands), computations will be much slower than in this tutorial. You may also need a computer with a lot of memory, even if the copy number matrix has been implemented as a sparse matrix here. This saves memory, but results in slower computations.

We hope to improve on all aspects of this methods in the near future.
