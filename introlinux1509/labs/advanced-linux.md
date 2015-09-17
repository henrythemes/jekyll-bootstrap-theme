---
layout: default
title:  'Advanced Linux'
---

# Advanced Linux
**NOTE:** in syntax examples, the dollar sign ($) is not to be printed.
The dollar sign is usually an indicator that the text following it should be typed in a terminal window.

[Answers:](advanced-linux-answers) If you did not finish the lab and you get stuck on any of the questions, here are the solutions.
No cheating though, try until you cry first.

## Ownership & Permissions
As Linux can be a multi-user environment it is important that files and directories can have different owners and permissions to keep people from editing or executing your files.

## Owners
The permissions are defined separately for **users**, **groups** and **others**.

The **user** is the username of the person who owns the file.
By default the user who creates a file will become its owner.

The **group** is a group of users that co-own the file.
They will all have the same permissions to the file.
This is useful in any project where a group of people are working together.

The **others** is quite simply everyone else's permissions.

## Permissions
There are four permissions that a file or directory can have.
Note the one character designations/flags, **r**,**w**,**x**  and  **-**.

In all cases, if the file or directory has the flag it means that it is enabled.

**Read: r**

File: Whether the file can be opened and read.

Directory: Whether the contents of the directory can be listed.

**Write: w**

File: Whether the file can be modified.
Note that for renaming or deleting a file you need additional directory permissions.

Directory: Whether the files in the directory can be renamed or deleted.

**Execute: x**

File: Whether the file can be executed as a program or shell script.

Directory: Whether the directory can be entered through "cd".

**No permissions: -**

## Interpreting the permissions of files and directories
Start your terminal, log onto UPPMAX (check with squeue which core you had and ssh onto it, if for some reason your core is gone :

```bash
$ salloc -A g2015037 -t 04:30:00 -p core --no-shell --reservation=g2015037_17 &

```

make an empty directory we can work in and make a file.

```bash
$ cd ~/glob/ngs-intro
$ mkdir advlinux
$ cd advlinux
$ touch  filename
$ ls -lh
total 0
 -rw-r--r-- 1 S_D staff 0B Sep 21 13:54 filename
 ```

(-lh means long and human readable, displaying more information about the files or directories in an humanly understandable format)

The first segment, -rw-r--r--, describes the ownerships and permissions of our newly created file.
The very first character, in this case "-", shows the files type.
It can be any of these:

**d** = directory
**-** = regular file
**l** = symbolic link
**s** = Unix domain socket
**p** = named pipe
**c** = character device file
**b** = block device file

As expected the file we have just created is a regular file.
Ignore the types other than directory, regular and symbolic link as they are outside the scope of this course.

The next nine characters, in our case rw-r--r--, can be divided into three groups consisting of three characters in order from left to right.
In our case rw- r-- and r--.
The first group designates the **users** permissions, the second the **groups** permissions and the third the **others** permissions.
As you may have guessed the within group permissions are ordered, the first always designates read permissions, the second write and the third executability.

This translates our files permissions to say this:

```bash
-rw-r--r--
```

"It is a regular file.
The user has read & write permission, but not execute.
The group has read permission but not write and execute.
Everyone else (other), have read permission but not write and execute."

As another example, lets create a directory.

```bash
$ mkdir directoryname
$ ls -lh
total 0
drwxr-xr-x  2 S_D  staff    68B Sep 21 14:41 directoryname
-rw-r--r--  1 S_D  staff     0B Sep 21 13:54 filename
```

As you can see the first character correctly identifies it as **d**, a directory, and all user groups have **x**, execute permissions, to enter the directory by default.

## Editing Ownership & Permissions
The command to set file permission is "chmod" which means "**CH**ange **MOD**e".
Only the owner can set file permissions.

1. First you decide which group you want to set permissions for.
User, **u**, group, **g**, other, **o**, or all three, **a**.
1. Next you either add, **+**, remove, **-**, or wipe out previous and add new, **=**, permissions.
1. Then you specify the kind of permission: **r**,**w**,**x**, or **-**.

Lets revisit our example file and directory to test this.

```bash
$ ls -lh
total 0
drwxr-xr-x  2 S_D  staff    68B Sep 21 14:41 directoryname
-rw-r--r--  1 S_D  staff     0B Sep 21 13:54 filename
$ chmod a=x filename
$ ls -lh
total 0
drwxr-xr-x  2 S_D  staff    68B Sep 21 14:41 directoryname
---x--x--x  1 S_D  staff     0B Sep 21 13:54 filename
```

As you can see this affected all three, **a**, it wiped the previous permissions, **=**, and added an executable permission, **x**, to all three groups.

Try some others both on the file and directory to get the hang of it.

```bash
$ chmod g+r filename
$ chmod u-x filename
$ chmod ug=rx filename
$ chmod a=- filename
$ chmod a+w directoryname
```

### Assignment
The assignments do not have to be handed in in any form.
Just complete them for your own sake, to make sure that you have understood the material.
If you need help or further explanations to complete any assignments, please do not hesitate to contact us, it is what we are here for =).

In no more than two commands, get the file permissions from

```bash
---------- 
```

to

```bash
-rw-rw--wx 
```

Notice also that we here gave everyone writing permission to the file, that means that ANYONE can write to the file.
Not very safe.

## Symbolic links - Files
Much like a windows user has a shortcut on his desktop to WorldOfWarcraft.exe, being able to create links to files or directories is good to know in Linux.

An important thing to remember about symbolic links is that they are not updated, so if you or someone else moves or removes the original file/directory the link will stop working.

Lets remove our old file and directory.
(Notice: rm -r * removes all folder and subfolders from where you are standing.
Make sure you are standing in our advlinux/ directory.)

```bash
$ rm -r *
```

Now that the directory is empty lets make a new folder and a new file in that folder.

```bash
$ mkdir stuff
$ touch stuff/linkfile
```

Lets put some information into the file, just put some text, anything, like "slartibartfast" or something.

```bash
$ nano stuff/linkfile
```
Now lets create a link to this file in our original folder.
**ln** stands for link and **-s** makes it symbolic.
The other options are not in the scope of this course.

```bash
$ ln -s stuff/linkfile
$ ls -l
total 8
lrwxr-xr-x  1 S_D  staff    14B Sep 21 15:38 linkfile -> stuff/linkfile
drwxr-xr-x  3 S_D  staff   102B Sep 21 15:36 stuff
```

Notice that we see the type of the file is l, for symbolic link, and that we have a pointer after the links name for where the link goes, -> stuff/linkfile.

### Assignment
I want you to change the information in the file using the link file, then check the information in the original file.
Did editing the information in the link change the information in the original?

Change the information using the original file, then checking the link will result in the same result.
The files are connected.

Now move or delete the original file.

What happens to the link?

What information is there now if you open the link?

Now create a new file in stuff/ with exactly the same name that your link file is pointing too with new information in it.

Is the link re-activated to the new file? (You can investigate this by seeing if the new information is seen when opening the link file).

## Symbolic links - Directories
Now lets create a link to a directory.

Lets clean our workspace.

```bash
$ rm -r *
```

Now lets create a directory three, arbitrarily, directories away.

```bash
$ mkdir -p one/two/three
```

Now lets enter the directory and create some files there.

```bash
$ cd one/two/three
$ touch a b c d e
$ ls -lh
total 0
-rw-r--r--  1 S_D  staff     0B Sep 21 16:11 a
-rw-r--r--  1 S_D  staff     0B Sep 21 16:11 b
-rw-r--r--  1 S_D  staff     0B Sep 21 16:11 c
-rw-r--r--  1 S_D  staff     0B Sep 21 16:11 d
-rw-r--r--  1 S_D  staff     0B Sep 21 16:11 e
```

Return to our starting folder and create a symbolic link to folder three

```bash
$ cd ../../..
$ ln -s one/two/three
$ ls -lh
total 8
drwxr-xr-x  3 S_D  staff   102B Sep 21 16:11 one
lrwxr-xr-x  1 S_D  staff    13B Sep 21 16:13 three -> one/two/three
```

Once again we see that it is correctly identified as a symbolic link, l, that its default name is the same as the directory it is pointing too, same as the files link had the same name as the file by default previously, and that we have the additional pointer after the links name showing us where it's going.

### Assignment
Perform "ls" and "cd" on the link.
Does it act just as if you were standing in directory two/ performing the very same actions on three/?

While having entered the link directory with "cd", go back one step using "cd ..", where do you end up?

Moving, deleting or renaming the directory would, just like in the case with the file, break the link.

## Grep - Searching for text
If you have a very large file, perhaps so large that opening it in program would be very hard on your computer.
It could be a file containing biological data, it could be a logfile of a transfer where we want check for any errors.
No matter the reason a handy tool to know the existence of is the "grep" command.

What grep does is search for a specific string in one or many files.
Case sensitive/insensitive or regular expressions work as well.

Lets start, as always, by cleaning our directory.

```bash
$ rm -r *
```

Then lets create a file with some text in it that we can grep around with.
I have supplied some great text below.

```bash
$ nano textfile 

Cats sleep anywhere, any table, any chair.
Top of piano, window-ledge, in the middle, on the edge.
Open draw, empty shoe, anybody's lap will do.
Fitted in a cardboard box, in the cupboard with your frocks.
Anywhere! They don't care! Cats sleep anywhere.
```

Now lets see how the grep command works.
The syntax is:

```bash
grep "string" filename/filepattern
```

Some examples for you to try and think about:

```bash
$ grep "Cat" textfile

$ grep "cat" textfile
```

As you can see the last one did not return any results.
Add a -i for case insensitive search.

```bash
$ grep -i "cat" textfile
```

Now lets copy the file and check both of them together by matching a pattern for the filenames.

```bash
$ cp textfile textcopy
$ grep "Cat" text*
```

The * will ensure that any file starting with "text" and then anything following will be searched.
This example would perhaps be more real if we had several text files with different texts and we were looking for a specific string from any of them.

### Assignment
Copy the file sample_1.sam to your folder using the command below

```bash
$ cp /proj/g2015037/labs/linux_additional-files/sample_1.sam .
```

Use grep to search in the file for a specific string of nucleotides, for example:

```bash
$ grep "TACCACCGAAATCTGTGCAGAGGAGAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAG" sample_1.sam
```

Try with shorter sequences.
When do you start getting lots of hits?

This file is only a fraction of a genome, you would have gotten many times more hits doing this to a complete many GB large sam file.

Use grep to find all lines with chr1 in them.
This output is too much to be meaningful.
Send it to a file (>) where you have now effectively stored all the chromosome 1 information.

## Piping
A useful tool in linux environment is the use of pipes.
What they essentially do is connect the first command to the second command and the second command to the third command etc for as many commands as you want or need.

This is often used in UPPMAX jobs and other analysis as there are three major benefits.
The first is that you do not have to stand in line to get a core or a node twice.
The second is that you do not generate intermediary data which will clog your storage, you go from start file to result.
The third is that it may actually be faster.

The pipe command has this syntax

```bash
command 1 | command 2
```

The "|" is the pipe symbol(on mac keyboard alt+7), signifying that whatever output usually comes out of command 1 should instead be directly sent to command 2 and output in the manner that command 2 outputs.

In a hypothetical situation you have a folder with hundreds of files and you know the file you are looking for is very large but you can't remember its name.
Lets do a ls -lh and pipe the output to be sorted by file size.
-n means we are sorting numerically and not alphabetically, -k 5 says "look at the fifth column of output", which happens to be the file size of ls command.

```bash
$ ls -lh | sort -k 5 -n
```

An example use would be to align a file and directly send the now aligned file to be converted into a different format that may be required for another part of the analysis.

The next step requires us to use a bioinformatics software called Samtools.
To be able to use this program we first have to load the module for it.
We will cover this in the UPPMAX lectures, so if you are a bit too fast for you own good, you will just have to type this command:

```bash
$ module load bioinfo-tools samtools
```

Here is an example where we convert the samfile to a bamfile (**-Sb** literally means the input is **S**am and to output in **b**am) and pipe it to immediately get sorted, not even creating the unsorted bamfile intermediary.
Notice that samtools is made to take the single "-" after samtools sort as the position of the piped data from samtools view.

```bash
$ samtools view -bS sample_1.sam | samtools sort - outbam
```

This should have generated a file called "outbam.bam" in your current folder.

We will have some more examples of pipes in the next section.

## Word Count
Word count, or wc, is a useful command for counting the number of occurrences of a word in a file.
This is easiest explained with an example.
Lets return to our sample_1.sam.

```bash
$ wc sample_1.sam
233288 3666760 58105794 sample_1.sam
```

This can be interpreted like this:

Number of lines = 233288

Number of words = 3666760

Number of characters = 58105794

To make this more meaningful lets use the pipes and grep command seen previously to see how many lines and how many times the string of nucleotides "CATCATCAT" exist in the file.

```bash
$ grep "CATCATCAT" sample_1.sam | wc
60  957 15074
```

To see only the line count you can add -l after wc and to count only characters **-m**.

### Assignment
Output only the amount of lines that have "chr1" in them from sample_1.sam.

Hard assignment: Count the lines that have "CATCATCAT" in them from outbam.bam

## Extra material 1

This is some harder assignments, so don't worry about it if you didn't have time to do it.

Lets look at grep and use some regular expressions http://www.cheatography.com/davechild/cheat-sheets/regular-expressions/

**Task1:** Using sample_1.sam find all lines that start with "@" and put them in a file called "at.txt".

**Task2:** Find all the lines that end with at least 3 numbers from at.txt.
(hint: remember, **sometimes** you have to escape {} with \\{\\})

## Extra material 2
Sed is a handy tool to replace strings in files http://www.grymoire.com/Unix/Sed.html

**Task:** You have realized that all the chromosomes have been misnamed as "chr3" when they should be "chr4".
Use sed to replace "chr3" with "chr4" in sample_1.sam and output it to sample_2.sam.

Some food for though: The solution to this replaces the first instance on each line of chr3, what if we have multiple instances? What if we had wanted to replace "chr1", this would effect chr10-19 as well! There are many things to consider :).

## Extra material 3
Bash loops are great for moving or renaming multiple files as well as many many other uses.

Create a couple of files as seen below

```bash
$ touch one.bam two.sam three.bam four.sam five.bam six.sam
```

**Task:** All the files are actually in bam format.
What a crazy mistake! Create a bash loop that changes all files ending in .sam to end with .bam instead.

Three hints:

1.The bash loop syntax is this:

```bash
$ for _variable_ in _pattern_; do _command with $variable_; done
```

2. To rename file1 to file2 you write this:

```bash
$ mv file1 file2
```

which effectively is the same thing as

```bash
$ cp file1 file 2
$ rm file1
```

3. Ponder how this can be used to your advantage:

```bash
$ i=filename
$ echo ${i/name}stuff
filestuff
```

Good luck and have fun! :)

If you get completely lost on some of the questions and there is nobody around to answer, here is a [cheat sheet](advanced-linux-answers) :) 
