---
layout: default
title:  'Part 5: Exploring your single cell genome'
---

# Part 5: Exploring your single cell genome

## 5.1 Contamination analysis in MEGAN

The production of a SAG is very sensitive to contamination. Even the slightest amounts of contaminating DNA can be amplified several thousands of times. 
A good way to detect contaminants in your assembly, is to compare your contigs to all the sequences in NCBI’s ‘nt’ database.  
To do this, perform the following BLASTn search:

```sh
blastn -query contigs.fasta -db /proj/g2015028/nobackup/single_cell_exercises/databases/nt -evalue 1e-5 -num_threads 8 -out contigs_vs_nt.blastn #[This shouldn’t take more than 4 minutes]
```

Once the search is done, take a quick look of contigs_vs_nt.blastn with less. Scroll a bit through the file by hitting space. 
When you are satisfied, hit q to quit less. This file contains all the information that we require to estimate which sequences in your assembly are 
contaminations, but as you may have guessed by looking at the file, it is an extremely laborious task to inspect each and every contig by hand. 
What is missing is a clear overview that allows you to quickly identify potential contaminations. This is where MEGAN comes in. 
MEGAN inspects the blastn file for you, and categorizes each contig to a certain taxonomic group.  
Enter the following command in your terminal window:  

```sh
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

**Q5.1a** Find the node 'Bacteria'. How much contigs are assigned to this node?  
If you don't see 'Bacteria' you can click on *Tree->Collapse* At *Taxonomic Level->Kingdom* and you should see the node.  

**Q5.1b** How many contigs are assigned to this node and all of its subnodes in total?
*Q5.1c** Which subnode of 'Bacteria' contains all of the remaining bacterial contigs not assigned to 'Bacteria'?  
There are two special nodes displayed: 'Not assigned' and 'No hits'. 
The 'No hits' node contains all the contigs that BLASTn could not find a significant hit for in the nt database. 
The 'Not assigned' node contains contigs for which BLASTn did find a hit for in nt, but MEGAN could not confidently assign it to a taxonomic node.  

**Q5.2** Why do you think BLASTn could not find a hit in nt for the majority of the contigs?  
MEGAN allows you to inspect a selected node in more detail. Right click the 'Thermodesulfobacterium' node and hit Inspect. 
In the window that pops up, uncollapse Thermodesulfobacterium. Now you see the identities of all the contigs that were assigned to this specific node. 
Here you can also find all the hits and BLASTn alignments.  
Close the window.
Explore all the taxonomic nodes to which MEGAN has assigned contigs.  

**Q5.3** Did you detect any contamination in your assembled dataset?  
If so, which taxa (or taxon) is (or are) the source of this contamination?  

**Q5.4** Of how much contigs and how much sequence in total (in kbp) does this contamination consist of?  

**Q5.5** What do you think could be the cause of this contamination?  
Close MEGAN and return to the terminal.

<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part5"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part5_2"> Next page</a></span>
</div>

