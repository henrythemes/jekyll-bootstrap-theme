---
layout: default
title:  'Connecting to UPPMAX'
---

# Part 1: Connecting to UPPMAX and reserving nodes

Before moving on to the exercises, you will first need to connect to Milou on Uppmax through ssh. 
If you have a Mac, first check if you will be able to use X forwarding. 
It is included in older macs but was removed in the latest releases where you need to install it (see: [http://support.apple.com/en-us/HT201341](http://support.apple.com/en-us/HT201341)). 
Then start the terminal (black screen icon). If you work on a PC, download and start MobaXterm (you can download it here: [http://mobaxterm.mobatek.net/](http://mobaxterm.mobatek.net/)). 
You will be using half a node (8 cores) per person. This means you will end up sharing one node for two people as one node consists of 16 cores.

It is essential to log in with X forwarding (-X) option to be able to display graphical interfaces across remote computers.
Type the following command but replace *username* with your login name.  

```bash
ssh -X *username*@milou.uppmax.uu.se
```

Next, request a compute node. Note that the reservation code for Wednesday is g2014179_wed.  

```bash
salloc -A g2014180 -t 08:00:00 -p core -n 8 --no-shell --qos=interact &
```

Do not repeat this command. Just type it once and see which node you are assigned to.Type the following command to see the login node:  

```bash
jobinfo -u *username*
```
or this command:

```bash
squeue -u *username*
```

The nodelist column gives you the name of the node that has been reserved for you (starts with 'q').  
Log in to compute node (replace mXX with the actual compute node you are assigned to)  

```bash
ssh -X mXX
```
Make sure that you can launch graphical tools in your node by typing this command:  

```bash
xclock
```

If you see the clock then you should be able to launch GUI-based tools such as MEGAN or Artemis remotely. Close the clock.
Next, before you are able to use a specific bioinformatics tool it needs to be first loaded using the 'module load' command.  
For example, You can view the list of available modules by typing:  

```bash
module avail
```

This will list all the basic tools available but you won't see any bioinformatics tools available.  
Type q to return to the command prompt. To see the bioinformatics tools installed on Milou, type:  

```bash
module load bioinfo-tools
```

Then type:  

```bash
module avail
```

Now, you’ll see the bioinformatics tools installed on Milou that are categorized by the type of main tasks they perform.  