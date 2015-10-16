---
layout: default
title:  'Kallisto'
---

#Kallisto and Sleuth (Beta format)

## Transcript-level quantification with Kallisto

*Kallisto* is an "alignment free" RNA-seq quantification method that runs very fast with a small memory footprint, so that it can be run on most laptops. It is a command-line program that can be downloaded as binary executables for Linux or Mac, or in source code format. For a first insight in the program, read [here](https://liorpachter.wordpress.com/2015/05/10/near-optimal-rna-seq-quantification-with-kallisto/) and for the preprint article, see [here](http://arxiv.org/abs/1505.02710).

Kallisto is geared towards quantification on the transcript (isoform) level, rather than the gene level (although the latter can also be done by post-processing Kallisto output.) However, read assignment to transcript isoforms cannot (in general) be done unambiguously, so there is an intrinsic "quantification noise" or variability in this process. Kallisto can thus be run either in a single step (which is very fast) or in "bootstrap" mode (which takes longer, but can be done on several processors in parallel) in order to get uncertainty estimates for the expression levels - a kind of error bars for the quantification process. Running with bootstraps is mandatory if you want to perform differential expression analysis of isoforms with Sleuth (see below). 

Kallisto is primarily meant for quantification of an existing set of FASTA sequences, that is, it does not perform transcript assembly and it cannot quantify the expression of novel transcripts that are not in the transcript index that you provide to it. With that said, you can of course use contigs from an assembly that you have produced in some other program in your Kallisto index. It would also be possible to use the software for e g metagenomics or metatranscriptomics quantification.

## Differential expression with Sleuth

*Sleuth* is a companion package for Kallisto which is used for differential expression analysis of transcript quantifications from Kallisto. While you could use other differential expression packages such as limma or DESeq2 to analyze your Kallisto output, Sleuth also takes into consideration the inherent variability in transcript quantification as explained above. Sleuth also allows the modeling of covariates such as batch, individual, tissue type etc. in the same way as DESeq2/edgeR/limma, which is useful for experimental designs with multiple varying factors. 

Unlike Kallisto, Sleuth is an R package. At the end of a Sleuth analysis, it is possible to view a dynamical graphical presentation of the results where you can explore the differentially expressed transcripts in an intuitive way.

## Example commands for running Kallisto

This is just showing you a suggestion for steps to run Kallisto. More documentation will probably be added later when we know for sure this is a good method for RNA-seq downstream analysis. (We are currently benchmarking both Kallisto and Sleuth!)

These steps are all done at the command-line.

We start by downloading Kallisto and sratools. Note that the latter is not needed for Kallisto to run - we want it in order to be able to download files from SRA in this examples.

		wget https://github.com/pachterlab/kallisto/releases/download/v0.42.3/kallisto_mac-v0.42.3.tar.gz

You can of course download a different version if you want to - and if you are using Linux, you need to replace "mac" with "linux" in the command above.

		tar zvxf kallisto_mac-v0.42.3.tar.gz 
		wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-mac64.tar.gz

Again, if you are not using Mac, you need to change the file name above to something appropriate from http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2.

		tar zxvf sratoolkit.2.5.2-mac64.tar.gz

Now we will download and merge human cDNA and ncDNA files from ENSEMBL in order to build a Kallisto transcript index. Note that we can concatenate the two gzipped files without unpacking them first. We use both the protein-coding transcripts and the non-coding ones to be able to capture more of the transcriptome.

	wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
	wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
	cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.rna.fa.gz

Now we can build the transcriptome index. Let's also time it to get a sense of how long it takes.
	
	time kallisto/kallisto index -i hsGRCh38_kallisto Homo_sapiens.GRCh38.rna.fa.gz

The next steps will take quite a while (maybe several hours). What they entail is (1) downloading FASTQ files from SRA and (2) quantifying the FASTQ files against our Kallisto index with bootstrapping for later use in Sleuth. You might want to either substitute your own FASTQ files here, or to run Kallisto without bootstraps (skip the -b 100 option). The samples we use here are from an experiment on prostate tumors, with three pairs of normal/tumor samples for a total of six samples (there are many more in the original study)

		for f in SRR057629 SRR057630 SRR057632 SRR057649 SRR057650 SRR057651
		do sratoolkit.2.5.2-mac64/bin/fastq-dump --split-files $f
		time kallisto/kallisto quant -i hsGRCh38_kallisto -t 4 -b 100 ${f}_1.fastq ${f}_2.fastq -o ${f}
		rm ${f}_?.fastq
		done	

In this example, we put "-t 4" so we can use up to four processors in the bootstrapping. You may want to modify this value according to the machine you are working on. If you wanted to run Kallisto without bootstraps and just get expression values on a pair of FASTQ files, you would run

		kallisto/kallisto quant -i hsGRCh38_kallisto <FILE1>.fastq <FILE2>.fastq -o <OUTPUT_DIR_NAME>

## Example commands for running Sleuth

Here we give an example workflow for a DE analysis in Sleuth. This is based on the prostate cancer scenario for which we downloaded files above, so if you are using other data, you will of course need to adjust some commands accordingly.

This part is done entirely in R, so start your R environment and begin by installing the dependencies. This only needs to be done the first time, of course.

		source("http://bioconductor.org/biocLite.R")
		biocLite("rhdf5")
		install.packages("devtools") 
		devtools::install_github("pachterlab/sleuth")
 
Now load the package and use a function that we borrowed from the Sleuth documentation for connecting ENSEMBL transcript names to common gene names. This will turn out to be useful at the end, when we look at the dynamic visualization of the results.

		library("sleuth")

		tx2gene <- function(){
		mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
		t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                	"external_gene_name"), mart = mart)
		t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     	ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
		return(t2g)
		}

		t2g <- tx2gene()

The actual Sleuth analysis starts by defining the path to the directory where your Kallisto output folders are. You will need to replace "/PATH/TO/YOUR/FOLDER" with the actual path on your machine. 

		base_dir <- "PATH/TO/YOUR/FOLDER" 
		samples <- c("SRR057629","SRR057630","SRR057632","SRR057649","SRR057650","SRR057651")
		kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

Now fill in metadata about the samples. In this case pretend we just know it, rather than fetching metadata from SRA, which is also possible in principle.

		s2c <- data.frame(sample=samples,individual=as.factor(rep(c(2,3,6),2)), condition=c(rep("tumor",3),rep("normal",3)))
		

The next command will read the Kallisto output files, connect them with metadata, and set up a linear model for analyzing the expression data.
 
		so <- sleuth_prep(kal_dirs, s2c, ~individual+condition, target_mapping = t2g)

Next we fit the linear model and test for one of the model coefficients, in this case "condition" (tumor vs normal).

		so <- sleuth_fit(so)
		so <- sleuth_test(so, which_beta = 'conditiontumor') 

Now we should be able to visualize the results:

		sleuth_live(so)

