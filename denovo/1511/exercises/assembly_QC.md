---
layout: default
title:  'Exercise: Pre-assembly Analyses'
---

##Exercise: Pre-assembly Analyses

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: henrik.lantz@bils.se

**Acknowledgements**

Credits to to Björn Nystedt, Doug Scofield, Nat Street, Francesco Vezzi, Amaryllis Vidali, Andrea Zuccolo and others! 

###Setting up your project environment

1.  Login to Uppmax and get your own node to run the lab

  First follow the instructions [here](../../common/login_instructions) to log in to UPPMAX.

2.  Go to your glob directory (in your home) and copy the files for the lab to your local folder

  ```
  cd ~/glob/
  rsync -r -v --progress  /proj/g2015027/assemblyQC . # ~2 min to copy cd assemblyQC/
  ```

3.  Load the bioinfo module: This will give you access to other bioinformatics modules used below.

  ```
  module load bioinfo-tools
  ```

You are now ready to start the analyses!

### FastQC

First we can get a number of statistics from the reads by running FastQC. This is dead simple and fast.

  ```
# Enter the course directory and create a new folder
cd ~/assemblyQC
mkdir FastQC
# Load FastQC
module load bioinfo-tools
module load FastQC
# Run FastQC separately on both read files
fastqc ../trim/Rbac_1.fq
fastqc ../trim/Rbac_2.fq
  ```

Transfer the html-reports to you own computer using scp. Ask the teachers if you are unsure how to.

#### Questions

- Do you see anything that looks problematic in the reports? If so, how can you improve the read files?

### Adapter trimming and QC filtering 

Much of the current large-scale assemblers are built on deBruijn kmer graphs. Hence, errors and adapters remained in your seqs can cause more problems than in e.g read mapping and SNP-calling, since base quality values are not taken into account in the assembly process.

Here, we will adopt a strategy to first search a subset of reads for adapters with Cutadapt (using a large collection of possible Illumina adapters). We will then use Trimmomatic to trim adapters and perform quality trimming based on the read quality scores. The reason for not running Cutadpat on the entire dataset is that it is quite slow if many adaptors are used in the search.

```
cd trim/
# Note that the Rhodobacter genome data is placed here (slightly edited to adher to older fastq format).
/proj/g2015027/assemblyQC/scripts/runCutadapt.sh
```

Take a look at the output file 'Rbac_cutadapt.report.short'.

We will now convert the Cutadapt output to a fasta file with adapters, and then run Trimmomatic using this file for trimming. We will at the same time also use Trimmomatic to trim the reads based on their quality values.

```
module add BioPerl
perl /proj/g2015027/assemblyQC/scripts/cutadaptReport2conf.pl Rbac_cutadapt.report.short > adapter_seqs.fa 
/proj/g2015027/assemblyQC/scripts/runFastQTrim.sh # ~6 min to run
```

Take a look at the html summary file produced by the script above: 'fastqc_summary/Rbac.html.

#### Questions

- How many adaptors were found in the read subset by Cutadapt?
- What fraction of the 25,000 reads in the subset had adapters?
- How many reads survived the quality trimming?
- Did the quality trimming have any major effect on the overall quality? Compare with the FastQC-report you did before trimming

###  Kmers

We will count the abundance (coverage) of each kmer in the complete set of reads using [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/). For a small dataset like this you can run it quickly (~3 min, <24 GB RAM). For large datsets, this can demand quite a lot of RAM.

Here we will run with the kmer size 31 bp:

```
# Goto the kmer folder
cd ~/glob/assemblyQC/kmer/
#
# Link the data files (both raw and filtered data; unpack the filtered)
ln -s ../trim/Rbac_1.fq .
ln -s ../trim/Rbac_2.fq .
gunzip -c ../trim/trimmomatic_dir/Rbac.trimmomatic_1.fq.gz > Rbac.trimmomatic_1.fq
gunzip -c ../trim/trimmomatic_dir/Rbac.trimmomatic_2.fq.gz > Rbac.trimmomatic_2.fq
#
# Load jellyfish and run the kmer analysis
module add jellyfish/1.1.11
# Raw data 
# I precomputed the "jellyfish count" run for you; takes too long for the lab
## jellyfish count -m 31 -c 4 -s 2G -t 8 --both-strands -o Rbac_raw_jelly --timing=rawTiming --stats=rawStats Rbac_1.fq Rbac_2.fq 
jellyfish histo -o Rbac_raw_jelly.hist Rbac_raw_jelly_0
# Trimmed data 
# I precomputed the "jellyfish count" run for you; takes too long for the lab
## jellyfish count -m 31 -c 4 -s 2G -t 8 --both-strands -o Rbac_trim_jelly --timing=trimTiming --stats=trimStats Rbac.trimmomatic_1.fq Rbac.trimmomatic_2.fq 
jellyfish histo -o Rbac_trim_jelly.hist Rbac_trim_jelly_0
```
To inspect the histograms you can plot them
Rscript --vanilla /proj/g2015027/assemblyQC/scripts/plotJelly.R
Inspect the histogram ('Rbac_jelly_hist.pdf'; copy to your own computer first, ask the teachers for help if needed)

#### Questions

- Does it look like you expected?
- What is a good cutoff for noise?
- What is the Cpeak value (the single-copy kmer coverage)
- What is the upper bound of the single-copy peak? 

 We will now estimate the genome size by using the formula

`Genome size = (total nb of kmers - noise kmers ) / Cpeak`

I've made a small script to make this calculation, along with a rough repeat estimation (the total number of kmers is found in the 'Stats' file produced by Jellyfis, while the rest of the parameters are taken from the histogram plot). You can run this for both the raw and trimmed data:

```
perl /proj/g2015027/assemblyQC/scripts/kmerStats.pl Rbac_raw_jelly.hist  <total nb of kmers> <noise cutoff> <Cpeak> <single-copy upper bound>
perl /proj/g2015027/assemblyQC/scripts/kmerStats.pl Rbac_trim_jelly.hist  <total nb of kmers> <noise cutoff> <Cpeak> <single-copy upper bound>
```

#### Questions
- What is the estimated genome size? Is there a difference between estimations from the raw and the trimmed data?
- What is the expected repeat content of the genome, based on 31-mers? 

You can check the true genome size of Rhodobacter sphaeroides here: http://gage.cbcb.umd.edu/results/index.html 

### Compare kmer-spectrum and GC-content

```
# compare read 1 vs read 2 or lib A vs lib B
# Density plot
cd ~/assemblyQC
mkdir kat
cd kat
kat comp -p -t 8 -C -D -o katout ../trim/Rbac_1.fq ../trim/Rbac_2.fq
# Spectra plot (must run density computation first)
kat plot spectra-mx -n -o katout_s.png katout-main.mx


# Compare GC content
kat gcp -t 8 -C -o katout ../trim/Rbac_1.fq ../trim/Rbac_2.fq

```


### De novo repeat library

One interesting aspect of a novel genome is the level of repeats. The repeat fraction differ enormously, sometimes even between relatively closely related species. Since repeats break up the assembly, and also tend to appear very fragmented themselves, it is wise to try to get an assembly-independent repeat library and repeat-quantification.

Here, we will use [Repeat Explorer](http://repeatexplorer.umbr.cas.cz/static/html/help/manual.html) to build a repeat library from a small subset of the reads (e.g. 0.1X genome coverage), and then use that library to estimate the repeat content of the genome directly from the reads.

Since the genome we are looking at contains very few repeats, we will use novel marine organism as a test case here instead. The organism has an estimated genome size of 1.4 Gbp (from kmer analyses). We will run Repeat Explorer on a subset of reads representing 0.01X coverage of the genome. Repeat Explorer is slow (even for this small subset) so I've done the run for you. The results are here:

```
cd ~/glob/assemblyQC/repeats/
```

You find two output folders, 'seqClust' and 'summary'. The most important output file is the reapeat library, ie a set of representative sequenced for high-copy repeat families

```
seqClust/assembly/contigs.info.minRD5_sort-GR
```

We can begin by running some stats on the repeat library.

```
module add BioPerl
/proj/g2015027/assemblyQC/scripts/contigStats.pl seqClust/assembly/contigs.info.minRD5_sort-GR 600000 > repeatlib.stats
```

#### Questions

- How many copies in the genome would you expect a repeat family to (at least) have to be present in the repeat lib?
- How many repeat seqs are there in the library?
- How large are they? (max, N50) 

I've also run a Repeatmasker job on a benchmark set of 500,000 independent reads from the same species. The results are in `repeatmaskerReads500k/`.

#### Questions
- How large fraction of all the reads were masked by the repeat library?
- What does it tell you about the repeat content of this organism?
- Why are the repeat hits not classified by type in Repeatmasker? 
