---
layout: default
title:  'Advanced Linux Answers'
---

# Advlinux, cheat sheet

#### Assignment – Ownership and permissions

Q:  
Change permission from 

```bash
----------
```

to

```bash
-rw-rw—wx
```

One answer:

```bash
chmod ug+rw filename  
chmod o+wx filename  
```

#### Assignment - Symbolic links - files
Q Did editing the information in the link change the information in the original?  
A: Yes.

Q: what happens to a symbolic link when we move whatever its pointing to?  
A: it breaks.

Q: is the link re-activated?  
A: yes.

#### Assignments - Grep
command to grep all lines with chr1 and send to other file.

ex:

```bash
grep "chr1" sample_1.sam > chr1.txt
```

#### Assignment - WC

```bash
grep "chr1" sample_1.sam | wc -l
samtools view outbam.bam | grep "CATCATCAT" | wc -l
```

#### Extra material 1.
Task1:

```bash
grep "^@" sample_1.sam > at.txt
```

Task2:

```bash
grep "[0-9]\{3\}$" sample_1.sam
```

#### Extra material 2.

```bash
sed 's/chr1/chr2/' sample_1.sam > sample_2.sam
```

#### Extra material 3

```bash
for f in *.sam; do mv $f ${f/.sam}.bam; done
```
