---
layout: default
title:  'Mapping reads'
---


#Mapping reads using a reference and converting them to the BAM format


###Data available for exercise

All fastqfiles if you do not have any your self that you will need for this exercise can be found in
 
	/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/referenceBased/data

on UPPMAX and through this [URL](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/referenceBased/data).

If you want to map more files for practice you can continue with the files found in
	
	/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/RAB11FIP5_fastqFiles

on UPPMAX and through this [URL](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/RAB11FIP5_fastqFiles).
 
A pre-build hg19 human genome for HiSAT is found here on uppmax
 
	/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/reference/hg19_hisat2

on UPPMAX and through this [URL](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_hisat2).

A pre-build hg19 human genome for STAR is found here on uppmax
 
	/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75

on UPPMAX and through this [URL](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75).
 

 

##Mapping short reads to a reference using HiSAT2

If you are using our data you will map the reads to the hg19 reference genome using a popular RNA-seq aligner **HiSAT2**. There are many features that can be tweaked using HiSAT2. For more information on all flags that can be used go [here](https://ccb.jhu.edu/software/hisat2/manual.shtml).

Read below for the flags we use for this exercise. Remember to change names accordingly so that you can run your program properly and know which files you have used.

To load the HiSAT2 package on Uppmax, execute::

     
     module use /proj/b2013006/sw/modules
     module load hisat2/2.0.0-beta

Now you can map the reads from one of the samples (or several; it's up to you which ones(s)) to map using a command such as the one below.

	mkdir outDir
    
    hisat2 -p N --dta-cufflinks -x path/to/reference/fileName -1 path/to/fastqFile/sample_1.fastq -2 path/to/fastqFile/sample_2.fastq -S /path/to/outDir/fileName.sam &> /path/to/outDir/fileName.sam.info 
    
    
flags used are 

*  ``-p N`` is the number of threads that will be used by the program.  
*  ``--dta-cufflinks`` will generate a output that is optimal for downstream analysis with cufflinks     
* ``-x /path/to/HiSAT2genome/fileName`` is the path to a pre-rendered reference library that HiSAT2 uses to map reads to the genome.``fileName`` is the part of the files in the folder without the endings ``.bt2``  
*  `` -1 /path/to/read1/sample_1.fastq `` is where you should add your forward fastq files that you will map to the reference.  
*  `` -2 /path/to/read1/sample_2.fastq `` is where you should add your reverse fastq files that you will map to the reference.  
*  ``-s /path/to/output/fileName.sam`` is the fileName of the samfile that will tell were the reads mapped to the reference.     
*  ``> /path/to/output/fileName.sam.info`` will give the statistics of how many of the reads that mapped that is generated from HiSAT2.  

This should run fairly quickly and put a file called ``fileName.sam`` in 
the directory that you specified with ``-s``. 


##Mapping short reads to a reference using STAR

Here we will map the reads to the hg19 reference genome using a popular RNA-seq 
aligner, **STAR**. There are many many features that can be tweaked using STAR. For more information concerning different features that can be used see the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

Read below for the flags we use for this exercise. Remember to change names accordingly 
so that you can run your program properly and know which files you have used.

To load the STAR package on Uppmax, execute

     module load bioinfo-tools
     module load star
     
  
Now you can map the reads from one of the samples (or several; it's up to you 
which ones(s)) to map using a command such as the one below.
  
	mkdir outDir
    
	STAR  --genomeDir /proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75  --readFilesIn path/to/fastqFile/sample_1.fastq path/to/fastqFile/sample_2.fastq --runThreadN N --outSAMstrandField intronMotif --outFileNamePrefix outDir/
	
flags used are 

* ``--genomeDir /path/to/STARgenome`` is the path to a pre-rendered reference library that STAR uses to map reads to the genome. 

*  ``--readFilesIn /path/to/read1 [/path/to/read2 ]`` is where you should add your fastq files that you will map to the reference.

*  ``--runThreadN N`` is the number of threads that will be used by the program.

*  ``--outSAMstrandField intronMotif`` adds SAM information so that unstranded reads function with cufflinks 

*  ``--outFileNamePrefix outDir`` tells you were the output data ends up. 


  
This should run fairly quickly and put a file called ``Aligned.out.sam`` in 
the directory that you specified with ``--outFileNamePrefix``. 



##Converting SAM files to BAM files



Next step is to convert the SAM  file to a BAM file and rename it to a more representable name. The most common tool used for this is [Samtools](http://www.htslib.org/doc/samtools.html). For more information regarding Samtools follow the link.  

To load the Samtools package on Uppmax, execute

     module load bioinfo-tools
     module load samtools


A good naming praxis is to name the file to correspond what you mapped. As an example if you mapped sample 12 using HiSAT2 you should rename the mapped SAM file to a file with the name ``sample12_RAB11FIP5.HiSAT2.bam``.  The renaming and BAMfile conversion can be done in one step. Then to view them you also have to sort the hits and index them. You can also get a report on your mapped reads using samtools **flagstat**. Since the BAM file contains all the information that you have in the SAM file remember to remove the sam file and the unsorted bam file once you are finished. 

	samtools view -bSh -o /path/to/outDir/properName.bam /path/to/fileName.sam
	#Go to outdir
	cd /path/to/outDir/
	#Convert to properName, sort, index and get stats
	
	samtools sort properName.bam  properName.sorted
	samtools index properName.sorted.bam
	samtools flagstat properName.sorted.bam > properName.sorted.flagstat
	
	#remove all the temporary files
	rm /path/to/fileName.sam
	rm /path/to/outDir/properName.bam


The sorted, indexed bam file can then be viewed in IGV. 

