---
layout: default
title:  'Part 5: Exploring your single cell genome'
---

# Part 5: Exploring your single cell genome

Now that you are fairly confident about the quality of your assembled single cell genome, 
you are going to have a look at the actual content of the assembly using the program MEGAN. 
MEGAN is a tool to visualize taxonomic and functional annotations identified through BLAST searches. 
It parses BLAST results and uses special algorithms and NCBI taxonomy to assign lineages to sequences with BLAST hits. 
This is useful for various reasons. It can be quite overwhelming to view result outputs from BLAST searches of metagenomic data, 
which can be very big usually. MEGAN condenses the output from these searches in a manageable view, summarizing sequence data into taxonomic groups. 
This makes it easier to compare large metagenomic data sets and associated functional annotations.

## 5.1 Contamination analysis in MEGAN
---

The production of a SAG is very sensitive to contamination. Even the slightest amounts of contaminating DNA can be amplified several thousands of times. 
A good way to detect contaminants in your assembly, is to compare your contigs to all the sequences in NCBI’s ‘nt’ database.  
To do this, perform the following BLASTn search:

```sh
module load blast/2.2.29+
blastn -query contigs.fasta -db /proj/g2014180/nobackup/single_cell_exercises/databases/nt -evalue 1e-5 -num_threads 8 -out contigs_vs_nt.blastn #[This shouldn’t take more than 4 minutes]
```

Once the search is done, take a quick look of contigs_vs_nt.blastn with less. Scroll a bit through the file by hitting space. 
When you are satisfied, hit q to quit less. This file contains all the information that we require to estimate which sequences in your assembly are 
contaminations, but as you may have guessed by looking at the file, it is an extremely laborious task to inspect each and every contig by hand. 
What is missing is a clear overview that allows you to quickly identify potential contaminations. This is where MEGAN comes in. 
MEGAN inspects the blastn file for you, and categorizes each contig to a certain taxonomic group.  
Enter the following command in your terminal window:  

```sh
module load MEGAN/4.70.4
MEGAN &
```

When you first launch MEGAN, it will ask you to enter the license. 
Usually, you need to register and obtain a license to use MEGAN here ( http://www-ab2.informatik.uni-tuebingen.de/software/megan/register/ ) 
but in this exercise, you can just click on 'skip entering the license' and it should take you to the main interface.
This will open MEGAN’s graphical user interface (GUI). Then select *File -> Import From BLAST*…
In the Import tab, click the folder icon. Then find *'contigs_vs_nt.blastn'* and hit Open. Click two times Next Step. 
In Step 3, click the save icon and specify the name of the MEGAN file that you want to create, *'contigs_vs_nt.rma'*. Hit Save and then Apply. 
Now MEGAN will assess all BLASTn hits and assign each of them to a taxonomic category. 
After MEGAN is done with loading, you will see taxonomic categorization of all your contigs visualized with a tree-like figure. 
To get more space in between branches of this figure, use the mouse scroll wheel. Click 'root'. 
A red box will appear to indicate root is selected, and the texts Ass=3 and Sum=197 (the number could be different on different assemblies). 
'Ass=' refers to the amount of contigs that have been confidently assigned to a given taxonomic node, 
and 'Sum=' refers to the total amount of contigs that have been confidently assigned to this node, plus any node below that node.
Right click 'root' and hit 'Uncollapse subtree'. This command will cause MEGAN to display the complete subtree below the selected node, 
in this case the root. The sizes of the nodes correspond to the relative amounts of contigs assigned to them.


### Questions
---

**Q5.1a** Find the node 'Bacteria'. How much contigs are assigned to this node? If you don't see 'Bacteria' you can click on *Tree->Collapse* At *Taxonomic Level->Kingdom* and you should see the node.  
**Q5.1b** How many contigs are assigned to this node and all of its subnodes in total?
*Q5.1c** Which subnode of 'Bacteria' contains all of the remaining bacterial contigs not assigned to 'Bacteria'?  
There are two special nodes displayed: 'Not assigned' and 'No hits'. The 'No hits' node contains all the contigs that BLASTn could not find a significant hit for in the nt database. The 'Not assigned' node contains contigs for which BLASTn did find a hit for in nt, but MEGAN could not confidently assign it to a taxonomic node.  
**Q5.2** Why do you think BLASTn could not find a hit in nt for the majority of the contigs?  
MEGAN allows you to inspect a selected node in more detail. Right click the 'Thermodesulfobacterium' node and hit Inspect. In the window that pops up, uncollapse Thermodesulfobacterium. Now you see the identities of all the contigs that were assigned to this specific node. Here you can also find all the hits and BLASTn alignments.  
Close the window.
Explore all the taxonomic nodes to which MEGAN has assigned contigs.  
**Q5.3** Did you detect any contamination in your assembled dataset?  
If so, which taxa (or taxon) is (or are) the source of this contamination?  
**Q5.4** Of how much contigs and how much sequence in total (in kbp) does this contamination consist of?  
**Q5.5** What do you think could be the cause of this contamination?  
Close MEGAN and return to the terminal.

## 5.2 Functional analysis of predicted protein content in MEGAN
---

In addition to being a great tool to detect contaminants, you can also use MEGAN to get a rough estimate of which proteins are encoded in your SAG. 
For this, MEGAN utilizes the KEGG (Kyoto Encyclopedia of Genes and Genomes) Pathways database in combination with a BLASTp file in which 
you searched NCBI's 'nr' database for protein sequences similar to the predicted protein sequences in your SAG (see 3c). 
MEGAN inspects the BLASTp outfile and 'maps' all BLASTp-hits to all KEGG Pathways. 
This allows you to get a clear overview of which biological pathways and processes are present in your SAG. 

Note however that the KEGG Pathways database does not contain all proteins and 
several proteins you might be interested in could be missed by this analysis.  
Execute the BLASTp search:  

```sh
blastp -query contigs.prodigal.faa -db /proj/g2014180/nobackup/single_cell_exercises/databases/db -evalue 1e-5 -num_threads 8 -out contigs.prodigal_vs_db.blastp
```

The database that you are searching now is a modified version of NCBI's 'nr'. 
We have severely reduced the size to reduce the search time to 10-15 minutes, while still finding all relevant proteins. 
Open MEGAN and import the BLASTp file as before but this time you need to select the Content tab and check Analyse KEGG content, before hitting Apply. 
Name your MEGAN file *'contigs.prodigal_vs_db.rma'*.
Explore the taxonomic classification using functions described in *5.1*.

### Questions
---

**Q5.6** Can you find the taxa (or taxon) that you identified as contamination in 5.1?  
Why (not)?  
**Q5.7** Apparently the fraction of ‘No hits’ is smaller in the BLASTp search compared to the BLASTn search. Can you think of a reason that explains this?  
Look for the KEGG button in the toolbar and hit it. A new window will pop up. Like with the taxonomy window, the KEGG Pathways are organized in a tree-like structure, making it easy to browse through all the different pathways.
For each pathway there is a corresponding graphical overview. Double click on the name of a pathway to see the graphical overview. Basically, each overview consists of the illustrated pathway and a list of proteins that make up this pathway. In these overviews it is possible to see which particular proteins of a given pathway are present in your SAG and which ones are absent. If present, the box that contains the name of protein is marked with green. The amount of proteins present in your SAG that ‘map’ to a given protein determine the intensity of the green color.  
Explore the KEGG Pathways.  
**Q5.8** Which RNA Polymerase subunits can you find?  
**Q5.9** Does this SAG encode genes for flagellar assembly? If so, which ones?   
**Q5.10** How many genes encoding for ribosomal proteins are present in the SAG?  
**Q5.11** Look at metabolism pathways. Would you say this bacterium is able to grow aerobically?  
**Q5.12** What other potentially interesting proteins can you find?  
Explain for each one why do you think they are interesting.  
**Q5.13** Based on what you find in the KEGG window in MEGAN, what can you say about the lifestyle of the bacterium from which this SAG was retrieved?  
For example, would you say it is motile, does it grow aerobically or anaerobically, is it free-living or symbiotic?  
Try to say as much as possible.  
**Q5.14** Given the estimated completeness of the genome (see 3f), how confident are you about statements regarding the nature of this bacterium?  
Close MEGAN and return to the terminal.  