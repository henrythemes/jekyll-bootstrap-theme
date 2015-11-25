---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

## 3.7. Identifying ribosomal RNAs

First, you will run this tool called **rnammer** to predict the regions within assembled contigs that contain ribosomal RNA sequences (rRNA), such as 16S, 23S, and 5S rRNA sequences.  
Go into the folder containing the assembled contigs and type this command (which you should already be in):  

```sh
#Fix some module dependencies (downgrade gcc, openmpi, hmmer)
module load rnammer/1.2
#End Fix
cd ..
mkdir 16s_blast
cp assemblies/*.fasta ./16s_blast
cd 16s_blast
rnammer -S bac -m lsu,ssu,tsu \
-gff G5_${sample}${trim}_IDBA_contig.rnammer.gff \
-f G5_${sample}${trim}_IDBA_contig.rnammer.fasta \
< G5_${sample}${trim}_IDBA_contig.fasta
```

Take a look at the file (*'contigs.rnammer.gff'*) produced by 'rnammer'.  
Can you identify the positive matches predicted by 'rnammer'?  
Can you interpret the result output?  
Next, you will run 'blastn' against the 'Silva' database. 
Depending on whether or not you have identified 16S or 23S sequences, you will need to run 'blastn' on a different database.  

If you have identified a 16S rRNA sequence (check the gff file), type:

```sh
DB=/proj/g2015028/nobackup/single_cell_exercises/databases/SSURef_NR99_115_tax_silva_trunc.dna.fasta
blastn -db $DB -evalue 1e-6 -num_threads 8 -query G5_${sample}${trim}_IDBA_contig.rnammer.fasta \
-out G5_${sample}${trim}_IDBA_contig.rnammer.16S_silva.blastn
```

If you have identified 23S rRNA sequence, then type:

```sh
DB=/proj/g2015028/nobackup/single_cell_exercises/databases/LSURef_115_tax_silva_trunc.dna.fasta
blastn -db $DB -evalue 1e-6 -num_threads 8 -query G5_${sample}${trim}_IDBA_contig.rnammer.fasta \
-out G5_${sample}${trim}_IDBA_contig.rnammer.23S_silva.blastn
```

After running Blastn, can you identify what organism G5 belongs to?

Now you can go back to the [questions](scg_part3)!

<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part3_5"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part4"> Next page</a></span>
</div>
