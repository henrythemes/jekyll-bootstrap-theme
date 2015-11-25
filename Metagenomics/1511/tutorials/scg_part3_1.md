---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

## 3.1. Organzing working folder

The following set of commands are to be typed in your compute node (for example mXX - look up using ```jobinfo -u username``` command). 
Make sure you are typing them in the compute node and not log in node. Go back to Part 1 to check how to log in to your compute node.  
Before starting the exercises, you should make a folder named *single_cell_exercises* in your home directory where the exercises will be run.
Then create 2 folders *dataset1* and *dataset2* where the raw data will be linked

```sh
mkdir ~/single_cell_exercises
cd ~/single_cell_exercises/
mkdir dataset1 dataset2
```

Next, make symbolic links of sequences in those folders:

```sh
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset1/* dataset1/
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset2/* dataset2/
```
**Please, do not modify those files**

Check that the data is present in those 2 datasets folders

```sh
ls dataset1
ls dataset2
```

You should now see 2 files per dataset, a forward fastq file ```_R1_001.fastq``` and its reverse ```_R2_001.fastq```  

Later in some commands we use the variables *sample* and *trim*, the following commands set those variables. 
We will also load the softwares we need to complete all this tutorial.  
**In case you loose your connection, you will need to redo this step again.**  

#### For *Hiseq* data without trimming:
```sh
sample=Hiseq
trim=''
cd ~/single_cell_exercises/dataset1
source /proj/g2015028/nobackup/single_cell_exercises/scripts/modules_load
```

#### For *Hiseq* data with trimming:
```sh
sample=Hiseq
trim=_Trimmomatic
cd ~/single_cell_exercises/dataset1
source /proj/g2015028/nobackup/single_cell_exercises/scripts/modules_load
```

#### For *Miseq* data without trimming:
```sh
sample=Miseq
trim=''
cd ~/single_cell_exercises/dataset2
source /proj/g2015028/nobackup/single_cell_exercises/scripts/modules_load
```

#### For *Miseq* data with merging:
```sh
sample=Miseq
trim=_Trimmomatic
cd ~/single_cell_exercises/dataset2
source /proj/g2015028/nobackup/single_cell_exercises/scripts/modules_load
```

Now we can check that all softwares have been correctly loaded, please type the following command in your terminal and check that all element of the list below is listed.  

```sh
module list
```

* bioinfo-tools
* trimmomatic/0.32
* spades/3.1.1
* IDBA/1.1.1-384
* quast/2.3
* prodigal/2.60
* blast/2.2.29+
* bwa/0.7.5a
* picard/1.92
* artemis/15.0.0
* MEGAN/4.70.4
* Ray/2.3.1-mpiio

**Please contact the teachers if one of those softwares do not appear in your list.**


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part3"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part3_2"> Next page</a></span>
</div>

