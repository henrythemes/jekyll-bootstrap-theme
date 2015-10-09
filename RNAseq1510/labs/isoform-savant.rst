==========
Using Savant
==========

It is up to you whether you'd like to install Savant on your own machine or run it on UPPMAX. If you do the latter, you need to login to Uppmax with X windows support, e g::

      ssh -X someuser@kalkyl.uppmax.uu.se

It would probably be a good idea to run the program from within an interactive session, in order not to disturb the login nodes too much::

      interactive -A g2013179 --qos=short -t 15:00

This will give you an interactive session for up to 15 minutes. 

Call a script for launching Savant::
       
      sh /proj/g2013179/private/nobackup/Savant/Savant.sh 

This should hopefully launch Savant in a window.

Start by loading the hg19 genome (included with the software) by selecting File > Select genome...


Importing the peptide track
===========================

Now you can import the track with identified peptides from the MS experiment. Select File > Load Track from File and navigate to ``/proj/g2013179/private/RNAseqWorkshop/reference`` (which could be annoyingly difficult) and select the file ``human_A431_global-TDA-FDR1pc_green-known_red-novel.txt`` When Savant asks if you would like to format the file, select yes and choose BED as the format in the menu. 

Now, click the newly created track, which will reveal a menu on the right hand side. Choose Appearance > Enable ItemRGB. Now, the peptides that correspond to known coding sequences will be colored green, and novel peptides will be colored red.


