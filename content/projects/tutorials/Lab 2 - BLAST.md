
---
# BLAST like a computer
---
In today's lab we are going to be working with BLAST - *Basic Local Alignment Search Tool*. All of you should have come across BLAST in previous courses (BIO206). That is because it is one of the most important bioinformatic tools. When most people are introduced to BLAST they go to the [NCBI website](http://blast.ncbi.nlm.nih.gov/Blast.cgi) where they enter the sequence of their gene or protein into a little box, click on few options and then press GO. That method is really useful when you have one sequence, maybe two -- but what do you do if you have 3000 sequences?

I want you to see the utilty of BLAST and to learn how to parse and interpret BLAST output. Unlike the user interface available at NCBI, we are going to be running BLAST on the BIO362 server from within this Notebook. Running BLAST on your local computer using your own BLAST databases is called **Standalone BLAST**. The challenge is that instead of creating an output on the screen, we  will create a file and tell Python to read it through. So the key Python skill you will learn is how to read through a file, line by line and make sense of a lot of biological data.

# Reading and Writing Files
---

You may get quite far into most Python books before you learn how to open and read a file. However, in Bioinformatics, reading in data from files is such an important skill that we are going to learn it now. We're lucky in biology that many of the types of data that we work with are stored in text files which are relatively simple to process using Python. Chief among these, of course, are DNA and protein sequence data, which can be stored in a variety of formats. But there are many other types of data – RNA/DNA sequencing reads, quality scores, genetic variation data, phylogenetic trees, geographical sample data, genetic mapping information – all of which we can access with Python.

Another reason for our interest in file input/output is the need for our Python programs to work as part of a pipeline or work flow involving other existing tools. When it comes to using Python in the real world, we often want Python to either accept data from, or provide data to, another program. Often the easiest way to do this is to have Python read, or write, files in a format that the other program already understands.



## Using `open()` to Read a File
---


What is a text file? In programming, text files don't have to be words on a pagee. Instead, they are characters and lines – something that you could open up and look at in a text editor, even if it is difficult to decipher.   


Normally when we want to open a file we just click on it with a mouse. In Python we use a function to open a file. Guess what it's called ... `open()`. This function requires at least one argument – the name of the file **as a string.**


```python
#RUN THIS BLOCK OF CODE
my_file = open("dna.txt",'r')
```

`open()` then returns a `file` object. A file object is a new type that we haven't encountered before, and it's a little more complicated than the strings, numbers and lists that we saw last week. With strings and numbers, it was easy to understand what they represented. A file object, in contrast, represents something a bit less tangible – I like to think of it as a file that is now open and awaiting your instructions. It hasn't been read or looked at in anyway - but it is waiting for you to do something.

Most of the time when we want to interact with a file object we use its special `methods`. Recall from last lab that a method is called like so:

    my_variable.method()
    
Similarly, file objects have their own methods for specific things you may want to do to a file. Read it, write to it etc.  

As a side note - this is exactly *why* Python uses methods - methods apply specifically to that one type of object. *ie* it doesn't makes sense to use the string method `replace()` on a number or a file object. This helps you keep things straight in your code, and helps guide the logic behind your manipulation of the objects you're working with.

The first thing we need to be able to do is to read the contents of the file. The file object has a `.read()` method which does this. It doesn't take any arguments, and it returns a string that is the entire contents of the file. We can store the file's contents in a  variable and treat them just like any other string – for example, we can print them:


```python
#RUN THIS BLOCK OF CODE

my_file = open("dna.txt")
file_contents = my_file.read()
print(file_contents)
```

    ATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCA
    


## Files, contents and filenames
---

When learning to work with files, it's very easy to get confused between a file object, a filename, and the contents of a file. Take a look at the following bit of code:


```python
my_file_name = "dna.txt"
my_file_object = open(my_file_name)
my_file_contents = my_file_object.read()
```

 What's going on here? 
 - On line 1, we store the string `dna.txt` in a simple string, a few letters  - nothing else  
 - On line 2 we use the variable `my_file_name` as the argument in the function `open()` and store the new file object in another variable `my_file_object`
 - On line 3, we use the `read()` method on `my_file`, and store the resulting data, in string format, in a new variable `my_file_contents`.
 
Confused? The important thing to understand about this code is that there are three separate variables that point to  three very different things:
 - `my_file_name` is a string, and it stores the name of a file on disk 
 - `my_file_object` is a file object, and it represents the open file.
 - `my_file_contents` is a string, and it stores the data that is in the file.


---
# Exercise 1
---
We have the protein sequence of the opsin1 gene associated with colour blindness stored in a file called `opsin.txt` that is stored in the Lab 2 folder. 

### Instructions:
 - Using the file `opsin.txt` count the number of amino acids in the protein. 
 - Store the count as an **integer** in the variable **`AA_count`**       


*Hint*: You may want to print the contents of the file to see what is in it



```python
opsin_file = open("opsin.txt")
opsin = opsin_file.read()
AA_count = len(opsin)
print(AA_count)
```

    348



```python
#CHECK YOUR ANSWER
from bio362 import lab2
lab2.count_amino_acids(AA_count)
```

    Correct!


---
## Looping through the lines of a file
---
Many Python objects have an obvious way you might want to loop through them in a `for` loop. Just like strings or lists, you can iterate through a file. The designers of Python knew that a lot of people would want to work line-by-line through a file so they built in this ability.


```python
#RUN THIS BLOCK OF CODE

file = open("disease_genes.txt")
for line in file:
    print(line)# do something with the line
```

    sickle cell	gene=HBB	protein=hemoglobin subunit beta
    
    colour blindness	gene=OPN1SW	protein=short-wave-sensitive opsin 1
    
    colour blindness	gene=OPN1LW	protein=long-wave-sensitive opsin 1
    
    colour blindness	gene=OPN1MW	protein=medium-wave-sensitive opsin 1
    
    albinism	gene=TYRP1	protein=5,6-dihydroxyindole-2-carboxylic acid oxidase precursor
    


Note that you are looping through the file object and not the contents. If you had read the contents and saved them in a string you would be iterating through the string, one character at a time.

This is a little confusing, but it is important to know that **reading through a file is a one-way street**. Once we have iterated over a file object, Python "remembers" that it is already at the end of the file, so when we try to iterate over it again, there are no lines remaining to be read. One way round this problem is to **close and re-open** the file each time we want to iterate over it:


```python
#RUN THIS BLOCK OF CODE

# print the length of each line
file = open("disease_genes.txt")
for line in file:
    print("The length is " + str(len(line)))
file.close()

# print the first character of each line
file = open("disease_genes.txt")
for line in file:
    print("The first character is " + line[0])
file.close()

```

    The length is 53
    The length is 66
    The length is 65
    The length is 67
    The length is 84
    The first character is s
    The first character is c
    The first character is c
    The first character is c
    The first character is a


Sometimes, if you want to move backwards through a file, or loop through it multiple times, a better approach is to read all of the lines of the file into a list. The `readlines()` method (for `fileObject`s) returns a list with each line from your file as an element.


```python
#RUN THIS BLOCK OF CODE

# first store a list of lines in the file
file = open("disease_genes.txt")
all_lines = file.readlines()
# print the lengths
for line in all_lines:
    print("The length is " + str(len(line)))
# print the first characters
for line in all_lines:
    print("The first character is " + line[0])
# Number of lines
print("The number of lines is " + str(len(all_lines)))
```

    The length is 53
    The length is 66
    The length is 65
    The length is 67
    The length is 84
    The first character is s
    The first character is c
    The first character is c
    The first character is c
    The first character is a
    The number of lines is 5


## `strip()` 
---
To tell a computer to start a new line, there is actually a special character at the end of each line called a *line break*. The most common type you will see is written like this: "`\n`"

Strings in Python have a method `strip()` that will remove not only the line breaks from the end of your sequence but all of the spaces, tabs or line breaks at the beginning or end of your string. It's often useful to use `strip()` if you are interested in the contents of the line but not the line break.

Similarly, there is `lstrip()` and `rstrip()` that remove the spaces, tabs or line breaks from left and right of your string.

---
# Exercise 2
---
The most common way to store sequence data is in a format called "FASTA". The protein that you are going to BLAST in Exercise 3 is stored as a FASTA sequence. Briefly, "FASTA" files are defined as follows:
1. A sequence starts with a single header line that consists of ">" character and the sequence name  
2. At the end of header line there is a linebreak
3. Subsequent lines are the DNA, RNA or amino acid sequence. 
4. Multiple sequences can be stored in one file with each new sequence starting with a new header line

For example:

    >Robs DNA sequence
    AGTCGATGCTAGCTAGCTACGATGCTAGCTGCTAGTCGATGCGATAGCTAGCTA
    AGATCGTAGCGATCGTAGCTAGTCGTAGCATGCTAGTCGATCGTACGATGCTCT

In the a file called "`albinism_protein.fasta`" there is a single protein. 

In this exercise use what you learned about opening files and `fasta` files to extract the name of the sequence.

### Instructions:
Store the name as a **string** in the variable **`protein_name`. **   
*Hint*: Don't forget the ">" in the sequence header is **not** part of the name


```python
albinism_file = open("albinism_protein.fasta")
albinism = albinism_file.readline()
protein_name = albinism[1:].strip()
print(protein_name)
```

    lcl|NG_011705.1_prot_NP_000541.1_1 [gene=TYRP1] [protein=5,6-dihydroxyindole-2-carboxylic acid oxidase precursor] [protein_id=NP_000541.1]



```python
from bio362 import lab2
lab2.name_that_protein(protein_name)
```

    Correct!


---
## Splitting a `string` to make a `list`
---
In last weeks lab we made some lists manually. However, there are plenty of functions and methods in Python that produce lists as their output. One such method that is particularly interesting to biologists is the `split()` method, which works on strings. `split()` takes a single argument, called a delimiter, and splits the original string wherever it sees the delimiter, producing a list. (By default, `split()` uses a *space* as a delimiter if you don't provide an argument).
Here's an example:


```python
#RUN THIS BLOCK OF CODE

names = "melanogaster,simulans,yakuba,ananassae"
species = names.split(",")
print(str(species))
```

    ['melanogaster', 'simulans', 'yakuba', 'ananassae']


We can see from the output that the string has been split wherever there was a comma leaving us with a list of strings. Of course, once we've created our list `species`, we can iterate over it using a loop, just like any other list.

It is **important** to know that if you split up a string in this way and there are digits/numbers in the string, Python will not try to convert these into integers for you. The `split` method will return a **list of strings**.


```python
#RUN THIS BLOCK OF CODE

data_line = "Chimpanzee 20 Humans 23"
data_values = data_line.split()
print(data_values)

chimp_value = data_values[1]
human_value = data_values[3]
# adding strings
print("Adding Strings: " + human_value + chimp_value)
```

    ['Chimpanzee', '20', 'Humans', '23']
    Adding Strings: 2320


Note that Python `concatenated` the strings '20' and '23' rather than adding the numbers 20+23. Unlike some other programming languages (e.g. R), Python will not try to guess what you want. Python is very literal. Therefore, if you know which items are numbers you can force convert them from strings to numbers. To get around this we can force the string to be an integer using the function `int()`:


```python
#RUN THIS BLOCK OF CODE

data_line = "Chimpanzee 20 Humans 23"
data_values = data_line.split()
print(data_values)

chimp_value = int(data_values[1])
human_value = int(data_values[3])

print("Adding numbers:", human_value + chimp_value)
```

    ['Chimpanzee', '20', 'Humans', '23']
    Adding numbers: 43


Now that you can loop through a file, `split` each line into pieces and convert those pieces into numbers - you have all the tools necessary to analyze some real data.

---
# BLAST Time
---

In this lab, you are going to be running BLAST from the command-line. Many of you will have never opened the command-line, so I want to explain it briefly. The **command-line**, or shell, is a tool that lets you run your computer using typed text commands rather than your mouse. You can edit files, delete them or move them to a new folder. In a way, it's like your normal operating system. Just like clicking on a program, you can also run programs from the command-line. The difference is that rather than clicking on an icon you type the name of the program and any other pieces of information the program may want in order to run.

While it might just sound annoying at first to have to type everything out, the power is that you can automate your tasks and easily record what you did. For this reason, most bioinformaticians spend most of their work time in the command-line. In fact, one of the primary uses of Python is to communicate data between programs run in the command line and assemble multiple programs together into workflows or pipelines.

In this lab we are going to run standalone BLAST (as opposed to web-based BLAST) from the command-line. It is really important to realize the distinction between Python and the command-line - they are NOT the same. The **command line** is more like an operating system (Mac or Windows) and **Python** is just one app in that operating system.

One of the great features of Jupyter Notebooks is that you can instruct the command-line from within the notebook. The downside of this is that people get confused between what is Python and what is command-line.

In this lab we will run BLAST in the command line using special cells that start with %%bash to tell Jupyter Notebook to send this command to the computer. The commands will create text files which we will analyze with Python. This workflow is a really common way that a bioinformatician will spend their day!

## Running BLAST 
---
To run programs in the command line is pretty simple. You type the name of the program, then you provide with information about what data to take as input and what to provide as output. 

Like with Python functions, command line input and output options are called arguments. These arguments are given to the program using special `flags`. It's simple to explain through an example:

In the command line, if you want to run BLASTp - you type:

    blastp

Now this command on its own is not enough. It's like opening Word and expecting your report to be finished. You need to tell `blastp` what to sequences you want to BLAST:

    blastp -query my_sequence.fasta

The `-query` is the flag that tells BLAST the input. Flags starting with a `-` and followed by a letter or word are the standard in command-line programs. You can provide BLAST with a variety of options, similar to those available on the NCBI graphical user interface

    blastp -query my_sequence.fasta -db NT -out my_sequences.VS.NT.txt

This command provide an input file (query) `my_sequences.fasta`, and a database (-db) to BLAST against `NR` (all Proteins in GenBank) and makes an output file called `my_sequences.VS.NT.txt`. 

### Enough chit-chat, let's run BLAST.
The cell below this is a command line call of BLAST - **it's *NOT* a Python cell!**. You can press play as usual and it will run. The little asterisk `In [*]` on the left indicates that the computer is working on a cell. Wait until this is done before proceeding - this is the first time you've given the computer a workout! 

The BLAST is comparing one of the proteins from week one against all known proteins. The top will be recorded


```bash
%%bash
blastp -query albinism_protein.fasta -db NR -out albinism_vs_NR.txt -outfmt 6

```

Let's review the components of this command
- `blastp` : That's the name of the program for blasting proteins against proteins
- `-query ../data/albinism_protein.fasta` : that tells blast to use the sequence data stored in albinism_protein.fasta  as a BLAST input query. 
- `-db NR` : this tells BLAST to use a file called NR as a protein database to search. NR is a standard database and stands for "Non Redundant" protein database
- `-out albinism_vs_NR.txt` : this flag tells BLAST to store the output in a file called albinism_vs_NR.txt
- `-outfmt 6` : the flag `outfmt` tells BLAST what output format to use. Here 6 means it will be a table (see below)

## Explanation of BLAST Output - You will need this for the Exercises below!!!
When the `blastp` run above is done you will have created a new blast output. We selected a tabular output format (`-outfmt 6`). The standard columns of the table are as follows: 


|Column |Python Index|  Column Name | Column Description |
|---|---|---|---|
| 1 | 0 | Query | Name of the query sequence | 
| 2 | 1 | Subject   | Name of the subject "hit" sequence | 
| 3 | 2 | Identity| Percent of identical residues in alignment | 
| 4 | 3 | Length | Length of Alignment | 
| 5 | 4 | Mismatches | Number of mismatches in alignment | 
| 6 | 5 | Gaps | Number of gap openings | 
| 7 | 6 | Query start | First residue of query in alignment | 
| 8 | 7 | Query end |Last residue of query in alignment | 
| 9 | 8 | Subject start | Residue of subject that aligns to query start | 
| 10 | 9 | Subject end | Residue of subject that aligns to query end | 
| 11 | 10 | E-value | The Expectation (E) value of the alignment | 
| 12 | 11 | Bit-score | The bit score of the alignment | 



---
# Exercise 3
---
The BLASTp output created above includes many hundreds of sequences from NR that align to our albinism protein. Each of these subjects is given a unique identifier called a [GenInfo Identifier](http://www.ncbi.nlm.nih.gov/genbank/sequenceids/)  or `GI`. The GI is encoded in the name  of the subject sequence between the first and the second "|" or pipe characters. e.g.

    gi|675616239|ref|XP_008936982.1| 

where `675616239` is the GI number.

### Instructions:

 - Extract all the GI numbers for all the hits in the blast output and store them in a list called **`GIs`**

*Hint*: the list will be huge. Perhaps look at one output line and learn to extract the GI before proceeding to the whole file. 

*Hint*: don't forget that you can start with an empty list and then use the `append` method to add to it.

*Hint*: although it is referred to as a GI number - it's not really a number, in that it is meaningless to perform numeric opertations on it. So leave each GI as a string


```python
albinism_output = open("albinism_vs_NR.txt")
alb_output = albinism_output.readlines()
GIs = []
for line in alb_output:
    splits = line.split("|")
#    print(splits[2])
    GIs.append(splits[2])
print(GIs)
```

    ['4507757', '823672962', '823672960', '823672958', '426361297', '1034184682', '675678080', '119579116', '332222648', '297684426', '635073609', '37513', '402897419', '795359366', '47115303', '795305588', '966992657', '724920506', '795579654', '982301566', '795385651', '817277698', '403272756', '675637829', '562871930', '661902962', '955530696', '640813616', '532084618', '829791140', '742190486', '982933191', '984125567', '470605239', '1016641124', '982933193', '752421212', '633265765', '466040278', '961738497', '633265767', '15150345', '555957907', '395819086', '586526977', '31341878', '110350673', '591346802', '507662427', '68161268', '759184376', '803007103', '927098112', '55775720', '55775740', '241177898', '344271133', '961738495', '284810606', '961738493', '440909960', '937631044', '880946940', '670994291', '68534994', '560972011', '511892271', '744588909', '284810608', '971175580', '602710914', '355918818', '591346798', '591346800', '119579113', '823432506', '743752238', '303227955', '824556499', '958705495', '584033919', '68161270', '953884557', '594656891', '664750545', '68161272', '332980554', '593728800', '940738529', '1034184684', '471366020', '946761761', '556734181', '194097413', '550822370', '940738533', '674079060', '641708390', '594064762', '532017951', '940738531', '585149975', '354496542', '193795147', '1047835559', '634854809', '1047835545', '488545297', '970724949', '946674888', '537229845', '537229846', '593728802', '118500941', '118500937', '118500935', '505767256', '507549371', '118500939', '157820905', '321140144', '1048405384', '507933658', '589947657', '820968465', '617641158', '395516023', '852778610', '612040483', '918603565', '585642551', '586467895', '13654241', '512951870', '74201819', '512951868', '23395752', '21493015', '149412970', '558138017', '591356811', '530607134', '706133203', '564227740', '541976395', '557285109', '695164304', '700385467', '704175706', '529443897', '704235742', '697835956', '884921561', '704522424', '704306545', '701318270', '6714632', '696963949', '698436708', '733925638', '568926408', '678002076', '959035930', '1004000156', '1019366682', '675616239', '157932124', '675413564', '576938267', '699663339', '700353050', '525028570', '698396868', '699694528', '937500762', '527274621', '669277015', '686593220', '697476094', '723123542', '768396509', '944213342', '960988009', '998738663', '543357177', '771915238', '699613146', '542158926', '543260011', '697012224', '700407897', '701325070', '705673987', '690453750', '902921263', '224091225', '683923038', '701419285', '637260791', '694847526', '355567766', '351695892', '504137693', '12230741', '975094911', '663269391', '926525765', '927126384', '466002357', '565310616', '1002587447', '704575098', '546271930', '719789356', '431898623', '1044442773', '676718382', '679188251', '677426189', '677496345', '677365483', '677302920', '678177509', '677240380', '679216779', '675308898', '678992621', '677226995', '677296681', '678200576', '683461594', '697464825', '676414674', '677094546', '677475351', '677550051', '677467617', '683462160', '677379857', '432120097', '678974202', '679142169', '676581689', '678133183', '677536369', '521025665', '678117435', '677420036', '678129408', '6226724']



```python
from bio362 import lab2
lab2.GI_collector(GIs)
```

    Correct!


## Batch Entrez
---
Interestingly - the list of GI numbers you just made can be used to download thousands of sequences from NCBI in what is called a ["*Batch Entrez*"](http://www.ncbi.nlm.nih.gov/sites/batchentrez) download. You simply put these GIs in a list, load it into the website, choose what components you want to download and save it to a new file. 

## Comparing online to standalone BLAST
---
For the sake of interest (*ie. this is optional*) you can compare the BLAST we just did to what you would find if you BLASTed it online. Open a new browser tab, to the [NCBI BLASTp page](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). Print the protein we just BLASTed `"../data/lab2/albinism_protein.fasta"` and copy the data into the BLAST Query Sequence window. If you look near the bottom left of the page you'  see an option to expand the page `"+Algorithm parameters"`. In here you can see many of the options for customizing your BLAST. All of these options and more are available in standalone BLAST. In fact, you can tell standalone BLAST to generate an output in `html` format so that you can open it with a web browser and it will look similar to what you see through NCBI. 

There are many options for standalone BLAST output, and with a little effort you can get all the information you get online in a computer readable format. In most well written command-line programs you can see all the options by typing the name of the program and then `-help` or `-h` or `--help`. This is analogous to going to the help menu of an App in Windows or Mac.

In JupyterNotebook if you want to run the help command, make a new cell and run this:

    %%bash 
    # the strange command above tells jupyter to run this command in the command-line and NOT with Python.
    # You'll notice its presence in all the BLAST calls we do
    blastp -help


# Comparing genes
---
Let's be serious, it's really not that impressive to BLAST *one* sequence against a database. Any first year can do that online. In this example we will kick it up a notch. 

We are going to BLAST a set of 2716 Chimpanzee mRNAs against 57405 Human mRNAs to observe how they align - that's about 156 million possible comparisons! 

Run the following command:


```bash
%%bash 

blastn -query Pan_troglodytes.ref_mRNA.fasta -db Homo_sapien.ref_mRNA -out chimp_v_human.mRNA.txt -outfmt 6 -num_alignments 1 -max_hsps 1 -evalue 1e-10 

```

We've added a few new options this time:
- `blastn` : Note that blast ends with and n this time - this version of BLAST is for comparing nucleotide sequences
- `-num_alignments 1` : this tells BLAST to only report the 1 best alignment for each query, rather than hundreds per query.
- `-max_hsps 1` : This tells BLAST to only report the top *high scoring pair* HSP - or the best local alignment. With 1 alignment and 1 hsp we expect at most 1 line in the table for each chimpanzee mRNA
- `-evalue 1e-10` : This tells BLAST to only output alignments with E-values less than 1e-10 which is 1x10^-10


---
# Exercise 4
---
Using the output from the above BLAST "`chimp_v_human.mRNA.txt`" calculate the average percent of bases that differ between human and chimp mRNAs.

Note that a decimal number (e.g. 3.6) is called a `float` in Python. `float` as in floating point number. It's different from an `integer` because it is not a whole number. You can force a `string` or an `integer` to be a `float` the same way we learned for `str()` and `int()`.

    float('1.6')  
    
Go back to the cell with explanation of BLAST outputs to help you identify the correct column to extract.
### Instructions:
- Store the percent of sites that **differ** between chimpanzee and human mRNAs *without the percent symbol* in a variable called **`percent_different `**
- Store the variable as a float

*Hint*: Last week we saw that we can grow a string. It's even easier in Python to keep a running number. You must first define the number at the start and then you can augment it. E.g.:

    counter = 0
    for line in file:
        counter = counter + 1

The syntax `counter = counter + 1` can be replaced with the shortcut `counter += 1`



```python
mRNA_file = open("chimp_v_human.mRNA.txt")
mRNA = mRNA_file.readlines()
percent_difference = []
for line in mRNA:
    b = line.split()
    c = float(b[2])
    percent_difference.append(c)
percent_different = 100-(sum(percent_difference)/len(percent_difference))

print(percent_different)
```

    1.0431000370508627



```python
from bio362 import lab2
lab2.test_percent_different(percent_different)
```

    Correct!


---
# Exercise 5
---
In the evolution of protein coding sequences insertion and deletion mutations are generally deleterious because they interrupt the protein with frame shifting mutations. 

Quantify the ratio of indels to single nucleotide substitutions between human and chimpanzee mRNAs using the same BLAST output, `chimp_v_human.mRNA.txt`

### Instructions:
 - Store the ratio in a variable called **`indels_to_substitions`** 
 - Store the answer as a float    
  
 
*Hint*: Go back to the cell with explanation of BLAST outputs to help you identify the correct column to extract.


```python
mRNA_file = open("chimp_v_human.mRNA.txt")
mRNA = mRNA_file.readlines()
subs = []
indels = []
for line in mRNA:
    a = line.split()
    b = float(a[4])
    c = float(a[5])
    subs.append(b)
    indels.append(c)
indels_to_substitions = float(sum(indels)/sum(subs))
print(indels_to_substitions)
```

    0.07526385427243915



```python
from bio362 import lab2
lab2.test_indel_to_substition_ratio(indels_to_substitions)
```

    Correct!


---
# Wrap up.
---
You've learned some really key skills today:
1. In many of your courses and in research (ROPs/BIO481) you'll commonly be asked to read through files and analyze or alter them. If you're doing a lot of silly find and replace commands in Word or Excel you might want to ask - how could I do this in Python?
2. You've learned how to loop through a table, split the columns and do analysis on them. If you are confident enough to use Python for these tasks you'll save time in the long run. Just yesterday I needed to make a list of all the students in the class that have were enrolled  and who added or dropped, and it was much easier to read the table with Python than mess around in Excel.
3. You've had a very light introduction to what command-line programs look like. The vast majority of real computational research tools are command-line. In undergrad you're not exposed to them because instructors either don't know how to run them, or they realize it's too much work to teach you! I've made a tutorial on how to use the command-line that I will release in the near future.
4. You've seen a powerful way to align thousands of sequences against thousands of other sequences, namely BLAST. BLAST is excellent for finding local alignments between distantly related sequences and searching databases. Next week we will see how you align multiple sequences into a single large alignment.




---
# Practice problems
---


1. Use the file `disease_genes.txt` and extract a list of all the gene names (2nd column) but don't include the bit that says "gene="
2. Repeat the above task but extract the protein and put it into a list
3. Using the blast output `chimp_v_human.mRNA.txt` count the average bit score of BLAST hits alignments
4. Find the shortest and longest alignment in the BLAST output `chimp_v_human.mRNA.txt`
5. Does: query_end - query_start always equal the alignment length? Explain the answer biologically
6. Write 3 different ways of reading all the lines in `disease_genes.txt` into a list (each line is an element of the list) using `.read()`, `.readlines()` and a for loop over the lines. 
7. Use the `chimp_v_human.mRNA.txt` BLAST output to determine whether on average human or chimp homologs are longer

