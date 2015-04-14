---
layout: default
title:  'Exercise - Functional annotation'
---

# Functional annotation

Functional annotation is the process during which we try to put names to faces - what do the genes do that we have annotated and curated? Basically all existing approaches accomplish this by means of similarity. If a translation product has strong similarity to a protein that has previously been assigned a function, the rationale is that the function in this newly annotated transcript is probably the same. Of course, this thinking is a bit problematic (where do other functional annotations come from...?) and the method will break down the more distant a newly annotated genome is to existing reference data. A complementary strategy is to scan for more limited similarity - specifically to look for the motifs of functionally characterized protein domains. It doesn't directy tell you what the protein is doing exactly, but it can provide some first indication.

In this exercise we will use a tool that combines the search for full-sequence simliarity by means of 'Blast' against large public databases with more targeted characterization of functional elements through the InterproScan pipeline. Interproscan is a meta-search engine that can compare protein queries against numerous databases, including the protein family databases PFam and ProDom. The output from Blast2Go can then be used to add some information to our annotation.
## Runnig Blast2Go using a GUI

Blast2Go is available over at [http://www.blast2go.com/b2ghome](http://www.blast2go.com/b2ghome). Here you can diretly launch the graphical user interface (i.e. GUI) and run your analysis. To launch Blast2Go, click on the big orange field to the right ('Please click here') - and then execute the downloaded file. If you are working on a Mac with the latest version of OSX, you may be greeted by an error pointing to security settings. You then need to go to your security panel and launch the application from there. An exception will be created and the GUI should start. Note that the startup may take a few minute while Blast2Go downloads some needed files to your computer. You may also have to make a 64bit browser your default browser (e.g. Safari or Firefox) and install a Java update.

The amount of memory you specify will depend on both how much RAM your computer has available and the size of the data you wish to analyze. For most test data sets, 1GB of Ram is probably fine. For a full proteome, you may wish to increase this number substantially though. However, for propper 'production' data sets, you should consider setting up your own Blast2Go mirror and run it from the command line (instructions are available on their website).

A 'full' Blast2Go analysis can run for several days and consume several GB of Ram. This is because it needs to submit all your proteins (or transcripts) to NCBI for blasting - and it can only do that in very small chunks - as well as run them against InterProscan, which in turn needs to launch a large number of queries against other databases.

Since we do not wish to spend too much time on this, we will again limit our analysis to chromosome 4. It is also robably best to choose the analysis with ab-initio predictions enabled (unless you found the other build to be more convincing). Maker produces a protein fasta file together with the annotation and this file should be located in your maker directory.

Copy the file to your computer and load it in Blast2Go.
### Loading sequences

*File-&gt;Load Sequences*

*Specify 'Protein' as sequence type.*

This will add all the protein sequences from the Fasta file to your Blast2Go analysis.
### Perform Blast searches

Next we will run all proteins against the SwissProt database (other data bases are also available of course). There are two strategies to do this:
#### From within the GUI:

This procedure is generally not recommend, since it is very very very very slow - use the command-line blast option discussed below. However, for completeness:

Blast -&gt; Run BLAST step

This will open an interface where you can specify the parameters of the Blast run:

*Set an Email address* (doesn't have to be your actual address though).

*Set the Blast program to 'blastp'* (we wish to run proteins against a protein database)

*Set the BlastDB to 'Swissprot'*

Leave the rest at their defaults and click on the arrow in the top left corner of the dialog (start button). This will start the analysis. All the data will be stored to a file, so you now need to specify a location to where Blast2Go should write the output.

The blast search in Blast2Go takes about 20-25secs per protein request - depending on how many sequences you have submitted, you can make a fairly educted guess regarding the running time.

Once the Blast search is finished, you can move on to the next step.
#### From the command line on Uppmax:

A more sensible, but slightly more technical approach is to run Blast on Uppmax and import the output into Blast2Go. This has a few added benefitis:

- it scales better, since you can run Blast on many cores

- it works even if the NCBI blast server is down or very slow

To run Blast on your data, use the Ncbi Blast+ package against a Drosophila-specific database (included in the folder we have provided for you, under blastdb/blastplus) - of course, any other NCBI database would also work:

*blastp -db path/to/blastdb -num\_threads 8 -query your\_proteins.fa -outfmt 5 -out blastoutput.xml*

The output can then be added to your Blast2Go analysis by doing:

*File -&gt; Import -&gt; Import Blast results -&gt; One XML file*

Since we are not using an offical NCBI blast database, we need to make a small change to the input rules:

Unckeck 'Join hit ID and description'.
## GO mapping

Blast2Go uses an internal database that related entries in public databases to established GO terms. By running the built-in option, we can add this information to our data set:

*Mapping -&gt; Run GO-Mapping step*

Again, this analysis can take quite some time. If you have a larger annotation project to run, you may want to consider setting up your own local Blast2Go database. This will cut down the running time significantly. Some information on how to do this can be found here: [http://www.blast2go.com/b2gsupport/resources/35-localb2gdb](http://www.blast2go.com/b2gsupport/resources/35-localb2gdb)
## Compile a first functional annotation

Once the GO-term mapping is done, proceed to compiling a first annotation:

*Annotation -&gt; Run Annotation step*
### Add InterproScan data

Blast searches provide an indication about potential homology to known proteins. Interproscan, on the otherhand, combines a number of searches for conserved motifs and curated data sets of protein clusters etc.

InterproScan can be run in a number of ways - through a website, from the command line on a linux server and of course from within Blast2Go:

Annotation -&gt; Interproscan - &gt; Run Interproscan

This will open a list of databases against which to search - the common choices here are HmmPfam, SuperFamil y and perhaps SignalPHmm. Since we can limit the number of sequences to run, it is fine to use all these databases (but perhaps only on 20-30 sequences).
## What's next?

Blast2Go is a convenient way to attach information to newly annotated transcripts. You can use it to visualize pathways, check for protein domains and other known sequence features. Next, you could write scripts of your own to merge Blast2Go output into your annotation. Incidentially, Maker comes with a small utility script that can take InterProscan output and add it to a Maker annotation file.

For this, we first need to extract the Interproscan output from your Blast2Go project:

File -&gt; Export -&gt; Export Interproscan Results -&gt; Export as Txt

If you now copy this file back to Uppmax and into the same folder where the corresponding Maker gene annotation lives:

*ipr\_update\_gff maker.gff interproscan.txt &gt; maker.with\_interpro.gff*

Where a match is found, the new file will now include a feature called Dbxref in the transcript feature field.

And of course, because of Makers' compatibility with GMOD standards, an annotation augmented in this way can be loaded into e.g. WebApollo and will save annotators a lot of work when e.g. adding meta data to transcript models.