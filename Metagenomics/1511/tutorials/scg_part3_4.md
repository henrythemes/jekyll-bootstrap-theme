---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

## 3.4. Assessing assembly quality using Quast


After assembling the reads into contigs, you will use this tool called 'Quast' to calculate the basic metrics such as the length of largest contig, N50, etc.  
To run 'Quast', first we will create a link to the different contig files in the assembly folder. 
Use the following commands to help yourself:  

```sh
cd assemblies
ln -s `pwd`/IDBA${trim}/contig.fa G5_${sample}${trim}_IDBA_contig.fasta
ln -s `pwd`/Spades${trim}/contigs.fasta G5_${sample}${trim}_Spades_contig.fasta
ln -s `pwd`/Ray${trim}/Contigs.fasta G5_${sample}${trim}_Ray_contig.fasta
```

then run quast on all assemblies at once:

```sh
quast.py -T 8 -o assembly_metrics *.fasta
```

Results from 'Quast' will be found in the *assembly_metrics* folder. Take a look at the files produced by the program. 
The summary file containing the stats is *report.txt*. You should see the assembly metrics such as N50, G+C%, largest contig, total length (i.e., total length of all contigs added).  

Report the results in the spreadsheet:

```sh
less assembly_metrics/report.txt
```


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part3_3"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part3_5"> Next page</a></span>
</div>

