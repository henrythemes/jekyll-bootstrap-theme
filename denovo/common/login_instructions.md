---
layout: default
title:  'UPPMAX Login Instructions'
---

###1. Connecting to UPPMAX

The first step of this lab is to open a ssh connection to UPPMAX. If you have a Mac, start the terminal (black screen icon). If you work on a PC, download and start [MobaXterm](http://mobaxterm.mobatek.net).

Now type (change username to your own username):

```
$ ssh -X username@milou.uppmax.uu.se
```

and give your password when prompted. As you type the password, nothing will show on screen. No stars, no dots. It is supposed to be that way. Just type the password and press enter, it will be fine.

You should now get a welcoming message from Uppmax to show that you have successfully logged in. 

###2. Getting a node of your own

Usually you would do most of the work in this lab directly on one of the login nodes at uppmax, but we have arranged for you to have half of one node (=8 cores) each to avoid disturbances. To get this reservation you need to use the salloc command like this:

```
$ salloc -A g2015027 -t 08:00:00 -p core -n 8 --no-shell --reservation=check_below &
```

where you should substitute “check_below” with one of these alternatives depending on the day.

**Monday:** g2015027_20151116  
**Tuesday:** g2015027_20151117  

Now check which node you got (replace username with your uppmax user name) like this:

```
$ squeue -u username
```

The nodelist column gives you the name of the node that has been reserved for you (starts with "q"). Connect to that node using:

```
$ ssh -X nodename
```

Note: there is a uppmax specific tool called jobinfo that supplies the same kind of information as squeue that you can use as well (`$ jobinfo -u username`).

You are now logged in to your reserved node, and there is no need for you to use the SLURM queuing system. You can now continue with the specific exercise instructions.

By Martin Dahlö, revised by Henrik Lantz, re-revised by Martin Norling. 
