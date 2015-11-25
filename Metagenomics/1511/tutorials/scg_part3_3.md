---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

We will run three different assemblers on the data. You will choose either HiSeq or MiSeq dataset and trim it or use raw reads depending on how you organize your group. Make sure your variables (sample and trimming) are set properly. You will use one dataset so you can have one assembly folder and create three subfolders for the Ray, IDBA and Spades results. You do not need to create the folders, the assemblers will do that.

```sh
mkdir -p assemblies/
```

## 3.3a. Assemble your data Using SPAdes

Spades is a prokaryotic genome assembler that was specifically designed to be able to handle uneven coverage in single cell datasets. It works with Ion Torrent, PacBio and Illumina paired-end, mate-pairs and single reads. Spades is based on a de Bruijn graph and involves removing graph structures that result from erroneus reads. It uses multiple kmer to optimize assembly of different regions. You can read up about Spades at the [website](http://bioinf.spbau.ru/spades). Spades offer excellent user support.

* if you work with the **raw** reads  

>```sh
>time spades.py --sc -t 8 -m 24 \
>-1 G5_${sample}_R1_001.fastq \
>-2 G5_${sample}_R2_001.fastq \
>-o assemblies/Spades
>```

* if you work with the **trimmed** reads  

>```sh
>time spades.py --sc -t 8 -m 24  \
>-1 trimmed/G5_${sample}${trim}_1P.fastq \
>-2 trimmed/G5_${sample}${trim}_2P.fastq \
>-s trimmed/G5_${sample}${trim}_U.fastq \
>-o assemblies/Spades${trim}
>```


*Notice: Those previous commands will launch SPAdes assembly but also check how long the assembly takes. After the assembly has completed, check the time it took for the program to run. You should look at the 'real' time. 
Record the time in the spreadsheet table.
In general, typing the command ```time``` before other commands will help you check how long the computation took.*


## 3.3b. Assemble your data Using IDBA-UD

IDBA offers a collection of assemblers of which IDBA-UD was designed to handle data with uneven coverage depth. It is also a de Bruijn graph based assembler that iterates over mutliple kmers. You can find the website [here](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/) and although manual is rather limited, there is a google group offering support. 

IDBA-UD does not take the quality files and in the first step you have to convert the fastq file into fasta format compatible with the assembler.

* if you work with the **raw** reads

>```sh
>fq2fa --merge G5_${sample}_R1_001.fastq G5_${sample}_R2_001.fastq G5_${sample}_idba_input.fa
>time idba_ud -r G5_${sample}_idba_input.fa -o assemblies/IDBA --maxk 124
>```

* if you work with the **trimmed** reads

>```sh
>fq2fa --merge trimmed/G5_${sample}${trim}_1P.fastq trimmed/G5_${sample}${trim}_2P.fastq trimmed/G5_${sample}${trim}_idba_input.fa
>time idba_ud -r trimmed/G5_${sample}${trim}_idba_input.fa -o assemblies/IDBA${trim} --maxk 124 
>```

## 3.3c. Assemble your data Using Ray

Ray is an assembler that can be highly parallelized and can therefore be a good option in terms of running times, especially if you have access to large number of CPUs on a computer cluster. Memory can be another limiting factor when running assemblies for large datasets. We did not do a comparison here to check memory usage for these prepared datasets but that is another thing to keep in mind for your future datasets. Ray also has a version for metagenomes but not for the single cell genomes characterized by the extreme coverage differences, so in this case it surves the purpose of a regular assembler used for comparison with the two specialized assemblers.

* if you work with the **raw** reads

>```sh
>time mpirun -n 8 Ray -k 31 -p G5_${sample}_R1_001.fastq G5_${sample}_R2_001.fastq -o assemblies/Ray &> assemblies/ray.log 
>```

* if you work with the **trimmed** reads

>```sh
>time mpirun -n 8 Ray -k 31 -p trimmed/G5_${sample}${trim}_1P.fastq trimmed/G5_${sample}${trim}_2P.fastq -o assemblies/Ray${trim} &> assemblies/ray.log 
>```


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part3_2"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part3_4"> Next page</a></span>
</div>

