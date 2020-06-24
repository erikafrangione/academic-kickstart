---
# Course title, summary, and position.
linktitle: Topic 1
summary: A showcase of my solutions to common bioinformatics problems
weight: 1

# Page metadata.
title: Overview
date: "2018-09-09T00:00:00Z"
lastmod: "2018-09-09T00:00:00Z"
draft: false  # Is this a draft? true/false
toc: true  # Show table of contents? true/false
type: docs  # Do not modify.

# Add menu entry to sidebar.
# - name: Declare this menu item as a parent with ID `name`.
# - weight: Position of link in menu.
menu:
  example:
    name: Topic 1
    weight: 1
---


# Solutions to Bioinformatics Python Problems
## Topic 1: Using DNA as a String

## 1. Find the length of a gene

#### Given: A DNA string *s*
#### Required: The length of *s*.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Use the len() function to find the length of the string
GeneLength = len(my_gene)
print(GeneLength)
```

    3302


## 2. Extracting Nucleotides from a gene

#### Given: A DNA string *s*
#### Required: The first 3 and last 3 nucleotides from *s*.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Extract the first 3 and last 3 nucleotides from the string using Python notation
first_3_nucleotides= my_gene[:3]
last_3_nucleotides= my_gene[-3:]

print(first_3_nucleotides)
print(last_3_nucleotides)
```

    ATG
    TAA


## 3. Concatenating bases in a gene

#### Given: A DNA string *s*
#### Required: Concatenate the first 10 bases to the last 10 bases of *s*.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Extract the first 10 and last 10 nucleotides from the string
first_10_nucleotides = my_gene[:10]
last_10_nucleotides = my_gene[-10:]

# Concatenate the extracted bases and store the result in a new variable
first_and_last_10 = first_10_nucleotides + last_10_nucleotides
print(first_and_last_10)
```

    ATGAGAAAAATTACTTTTAA


## 4. Calculating GC content in a gene

#### Given: A DNA string *s*
#### Required: Calculate the GC content of *s* and store it as a percentage.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Use the count() method to return the number of G's and C's in the gene
G_count = my_gene.count('G')
C_count = my_gene.count('C')
print(G_count)
print(C_count)

# Find the total length of the gene
Total_nuc = len(my_gene)

# Concatenate the G and C content, divided by the total length of the gene, and mulitplied by 100
GC_content = ((G_count + C_count) / Total_nuc) * 100
print(GC_content)
```

    810
    794
    48.576620230163535


## 5. Find the reverse complement of a gene

#### Given: A DNA string *s*
#### Required: Store the reverse complement of *s* as an uppercase string.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Convert the gene to lowercase
lower = my_gene.lower()
print(lower)

# Use the replace() method to convert each lowercase nucleotide to their uppercase complement
replace_letters = lower.replace("a", "T").replace("t", "A").replace("c", "G").replace("g","C")

# Reverse the DNA string to properly record the opposite strand
reverse_comp = replace_letters[::-1]
print(reverse_comp)
```

    atgagaaaaatgtcggaggaagagttttatctgttcaaaaatatctcttcagtggggccgtgggatgggcctcagtaccacattgcccctgtctgggccttctacctccaggcagctttcatgggcactgtcttccttatagggttcccactcaatgccatggtgctggtggccacactgcgctacaaaaagttgcggcagcccctcaactacattctggtcaacgtgtccttcggaggcttcctcctctgcatcttctctgtcttccctgtcttcgtcgccagctgtaacggatacttcgtcttcggtcgccatgtttgtgctttggagggcttcctgggcactgtagcaggtactgcaggggaaaaaggggttaggggaaggcaaaggttgctactccactggagggggttcctaaagaggagtctgggggaaatgagtctggtgctttttaaaatactggggtacaaagcaacccagactagaagtttggctaaaataggatgtttgagtcttcactccaaatgtcagtccagccctgtctctctgtgcttgctccacccgatctgtttgccactctgccagccaggctgggtggggctctgtctagcccattatcctcacatttcaccacaggtctggttacaggatggtcactggccttcctggcctttgagcgctacattgtcatctgtaagcccttcggcaacttccgcttcagctccaagcatgcactgacggtggtcctggctacctggaccattggtattggcgtctccatcccacccttctttggctggagccggtgagagtgcagggcagtggtgctgagttaactaggagctcaggttgatgtgggtggaaagagagcttgggtataactatttagtctttgacctctacttttaaagagttgcaatatgaggcgaaagggcagtgggagacaagtgctaacgtttactctgcagttggaattgctgtagcttctcccagtcaggacagaaaacccccctgcttgaagccttagggcattccgtgggttctaagtggagaacacaatccaggcatctcagctcccactgcactcttgtggagagtccagtgagcaagtgtttggtcctttgcaggttcatccctgagggcctgcagtgttcctgtggccctgactggtacaccgtgggcaccaaataccgcagcgagtcctatacgtggttcctcttcatcttctgcttcattgtgcctctctccctcatctgcttctcctacactcagctgctgagggccctgaaagctgtgagtggcatttgatagtcagggaagaaggggttcggggctccacatgagaaggaagagtgctctgaaacataagatgcctggaaatgtccatagccagagagggtatctaaaagcagcaaaggaagtaggaggagggagaatgatggagatccaaaggaactaggccaggagatgggacagaaaagaggcaatcagagtggatgccccctcccccatcccacagaaaagcatccagagaccgggcgcagtggctcacgcctataatcccagcactttgggaggccgaagcagacggatcacctgaggtcaatagttccagaccagcctggcctacatggcaaaatgctaaaatgcgaaaattagctgggcatggtggtgtgtgcttgtaaccccagctactcaggaggctgagacaaaagaatcacttgaacccgggaggtggaggttcagtgagccgagactgcaccactgcaactccagcctgggcaacagagcgagactcggtctcaaaaagaaaattaaaaattaaaaattaaaaaaaaaaaaaaaaagcatccagagggccaggaaaaagagagatgtgatgctttccgtgctccaccccaggttgcagctcagcagcaggagtcagctacgacccagaaggctgaacgggaggtgagccgcatggtggttgtgatggtaggatccttctgtgtctgctacgtgccctacgcggccttcgccatgtacatggtcaacaaccgtaaccatgggctggacttacggcttgtcaccattccttcattcttctccaagagtgcttgcatctacaatcccatcatctactgcttcatgaataagcaggtaaagctctttattcacattcctatggtccagaagaccctggttcttttctcaccattgacttttaactcagagcaccctggactctacccaggtttctagtagacgagggaagccacaaaacccccgagtaggttgggaagcctttggtaagcacagggaggaaggcacggttatcaagacgagaaaatagaaccccggaggaaagaacttgagtcaggaaaatgaagttgctccaaagaacaggatgaatgaaagcattttattgaaaaactcgtgcagcaaaccaccatggcacacgtttacctatgtaacaaacctgcacatcctgcacatgtatcctggaacttaaataaaattaaaaaaataaaaataaaaactcagattcctctcaattttcagtccttgcatttaataatttcttaatcatttcccttccaacttttagcctgcacgagcatgtgtgaagcacagaaatcataccacatgcaaaaatctctaaaatatcttatcatctgaaggtactgggggatttcctatcccatctgaaatccgagctaataaacaccaaaccctaagtggcaaaaaccctactttcagatggtattgtttcctcaatcccagaggtagactcaaaactaatttgaaacctccctggatagaagagaattggcagtcctttccagctgggagcacctgctagtaatggaggggcctctgctgacagtgcttttatgaagcaggatggtttgtgaatttaccaacagtgaggtctcagacttgaccagtttaggattaccgtagacccaggagtagttctagactggaatctagatagttttcaggatggggaagatagattcaaaaccacctaagggcattctgggtacaaagcattgtgcaaggctttggtgatacagagaataaggtcttttttcccatacttcctcatctgccaaggttatctccaattgtacctttctctccagttccaagcttgcatcatgaagatggtgtgtgggaaggccatgacagatgaatccgacacatgcagctcccagaaaacagaagtttctactgtctcgtctacccaagttggccccaactgaggacccaatattggcctgtttgcaacagctagaattaaattttacttttaa
    TTAAAAGTAAAATTTAATTCTAGCTGTTGCAAACAGGCCAATATTGGGTCCTCAGTTGGGGCCAACTTGGGTAGACGAGACAGTAGAAACTTCTGTTTTCTGGGAGCTGCATGTGTCGGATTCATCTGTCATGGCCTTCCCACACACCATCTTCATGATGCAAGCTTGGAACTGGAGAGAAAGGTACAATTGGAGATAACCTTGGCAGATGAGGAAGTATGGGAAAAAAGACCTTATTCTCTGTATCACCAAAGCCTTGCACAATGCTTTGTACCCAGAATGCCCTTAGGTGGTTTTGAATCTATCTTCCCCATCCTGAAAACTATCTAGATTCCAGTCTAGAACTACTCCTGGGTCTACGGTAATCCTAAACTGGTCAAGTCTGAGACCTCACTGTTGGTAAATTCACAAACCATCCTGCTTCATAAAAGCACTGTCAGCAGAGGCCCCTCCATTACTAGCAGGTGCTCCCAGCTGGAAAGGACTGCCAATTCTCTTCTATCCAGGGAGGTTTCAAATTAGTTTTGAGTCTACCTCTGGGATTGAGGAAACAATACCATCTGAAAGTAGGGTTTTTGCCACTTAGGGTTTGGTGTTTATTAGCTCGGATTTCAGATGGGATAGGAAATCCCCCAGTACCTTCAGATGATAAGATATTTTAGAGATTTTTGCATGTGGTATGATTTCTGTGCTTCACACATGCTCGTGCAGGCTAAAAGTTGGAAGGGAAATGATTAAGAAATTATTAAATGCAAGGACTGAAAATTGAGAGGAATCTGAGTTTTTATTTTTATTTTTTTAATTTTATTTAAGTTCCAGGATACATGTGCAGGATGTGCAGGTTTGTTACATAGGTAAACGTGTGCCATGGTGGTTTGCTGCACGAGTTTTTCAATAAAATGCTTTCATTCATCCTGTTCTTTGGAGCAACTTCATTTTCCTGACTCAAGTTCTTTCCTCCGGGGTTCTATTTTCTCGTCTTGATAACCGTGCCTTCCTCCCTGTGCTTACCAAAGGCTTCCCAACCTACTCGGGGGTTTTGTGGCTTCCCTCGTCTACTAGAAACCTGGGTAGAGTCCAGGGTGCTCTGAGTTAAAAGTCAATGGTGAGAAAAGAACCAGGGTCTTCTGGACCATAGGAATGTGAATAAAGAGCTTTACCTGCTTATTCATGAAGCAGTAGATGATGGGATTGTAGATGCAAGCACTCTTGGAGAAGAATGAAGGAATGGTGACAAGCCGTAAGTCCAGCCCATGGTTACGGTTGTTGACCATGTACATGGCGAAGGCCGCGTAGGGCACGTAGCAGACACAGAAGGATCCTACCATCACAACCACCATGCGGCTCACCTCCCGTTCAGCCTTCTGGGTCGTAGCTGACTCCTGCTGCTGAGCTGCAACCTGGGGTGGAGCACGGAAAGCATCACATCTCTCTTTTTCCTGGCCCTCTGGATGCTTTTTTTTTTTTTTTTTAATTTTTAATTTTTAATTTTCTTTTTGAGACCGAGTCTCGCTCTGTTGCCCAGGCTGGAGTTGCAGTGGTGCAGTCTCGGCTCACTGAACCTCCACCTCCCGGGTTCAAGTGATTCTTTTGTCTCAGCCTCCTGAGTAGCTGGGGTTACAAGCACACACCACCATGCCCAGCTAATTTTCGCATTTTAGCATTTTGCCATGTAGGCCAGGCTGGTCTGGAACTATTGACCTCAGGTGATCCGTCTGCTTCGGCCTCCCAAAGTGCTGGGATTATAGGCGTGAGCCACTGCGCCCGGTCTCTGGATGCTTTTCTGTGGGATGGGGGAGGGGGCATCCACTCTGATTGCCTCTTTTCTGTCCCATCTCCTGGCCTAGTTCCTTTGGATCTCCATCATTCTCCCTCCTCCTACTTCCTTTGCTGCTTTTAGATACCCTCTCTGGCTATGGACATTTCCAGGCATCTTATGTTTCAGAGCACTCTTCCTTCTCATGTGGAGCCCCGAACCCCTTCTTCCCTGACTATCAAATGCCACTCACAGCTTTCAGGGCCCTCAGCAGCTGAGTGTAGGAGAAGCAGATGAGGGAGAGAGGCACAATGAAGCAGAAGATGAAGAGGAACCACGTATAGGACTCGCTGCGGTATTTGGTGCCCACGGTGTACCAGTCAGGGCCACAGGAACACTGCAGGCCCTCAGGGATGAACCTGCAAAGGACCAAACACTTGCTCACTGGACTCTCCACAAGAGTGCAGTGGGAGCTGAGATGCCTGGATTGTGTTCTCCACTTAGAACCCACGGAATGCCCTAAGGCTTCAAGCAGGGGGGTTTTCTGTCCTGACTGGGAGAAGCTACAGCAATTCCAACTGCAGAGTAAACGTTAGCACTTGTCTCCCACTGCCCTTTCGCCTCATATTGCAACTCTTTAAAAGTAGAGGTCAAAGACTAAATAGTTATACCCAAGCTCTCTTTCCACCCACATCAACCTGAGCTCCTAGTTAACTCAGCACCACTGCCCTGCACTCTCACCGGCTCCAGCCAAAGAAGGGTGGGATGGAGACGCCAATACCAATGGTCCAGGTAGCCAGGACCACCGTCAGTGCATGCTTGGAGCTGAAGCGGAAGTTGCCGAAGGGCTTACAGATGACAATGTAGCGCTCAAAGGCCAGGAAGGCCAGTGACCATCCTGTAACCAGACCTGTGGTGAAATGTGAGGATAATGGGCTAGACAGAGCCCCACCCAGCCTGGCTGGCAGAGTGGCAAACAGATCGGGTGGAGCAAGCACAGAGAGACAGGGCTGGACTGACATTTGGAGTGAAGACTCAAACATCCTATTTTAGCCAAACTTCTAGTCTGGGTTGCTTTGTACCCCAGTATTTTAAAAAGCACCAGACTCATTTCCCCCAGACTCCTCTTTAGGAACCCCCTCCAGTGGAGTAGCAACCTTTGCCTTCCCCTAACCCCTTTTTCCCCTGCAGTACCTGCTACAGTGCCCAGGAAGCCCTCCAAAGCACAAACATGGCGACCGAAGACGAAGTATCCGTTACAGCTGGCGACGAAGACAGGGAAGACAGAGAAGATGCAGAGGAGGAAGCCTCCGAAGGACACGTTGACCAGAATGTAGTTGAGGGGCTGCCGCAACTTTTTGTAGCGCAGTGTGGCCACCAGCACCATGGCATTGAGTGGGAACCCTATAAGGAAGACAGTGCCCATGAAAGCTGCCTGGAGGTAGAAGGCCCAGACAGGGGCAATGTGGTACTGAGGCCCATCCCACGGCCCCACTGAAGAGATATTTTTGAACAGATAAAACTCTTCCTCCGACATTTTTCTCAT


## 6A. Cleaving DNA into Restriction Fragments

#### Given: A DNA string *s*
#### Required: Find the position of the first restriction site "CCATGG" in *s* (i.e., CCA/TGG blunt end cutter)


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Use the find() method to identify the position of the first restriction site
restriction_site = my_gene.find("CCATGG")
print(restriction_site)
```

    157


## 6B. Count the Number of Restriction Fragments

#### Given: A DNA string *s*
#### Required: Count the number of bands that would appear on a gel after digesting *s* with a restriction enzyme that cuts at the same restriction site(s) mentioned in Problem 6A.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Use the count() method to identify all restriction sites
restriction_sites = my_gene.count("CCATGG")
print(restriction_sites)

# For a linear strand of DNA, the number of bands would be equal to the number of restriction sites + 1
number_of_bands = restriction_sites + 1
```

    3


## 7. Storing Exons in Lists

#### Given: A DNA string *s*
#### Required: Store each of the given exon lists extracted from *s* into a single larger list.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Provided start and end postions of 5 different exons within the DNA string
Exon_1 = [1, 352]
Exon_2 = [637, 805]
Exon_3 = [1128, 1293]
Exon_4 = [1903, 2142]
Exon_5 = [3132, 3251]

# Join all given exons into a list of lists, in order of position 1 through 5
exons = [Exon_1,Exon_2,Exon_3,Exon_4,Exon_5]
print(exons)
```

    [[1, 352], [637, 805], [1128, 1293], [1903, 2142], [3132, 3251]]


## 8.  Splicing a gene

#### Given: A DNA string *s*
#### Required: Given the previous exon coordinates for *s*, extract the exons and splice (join) them together to make a protein coding gene.


```python
# Read in DNA string
file = open("my_gene.txt")
my_gene = file.read().rstrip()

# Provided start and end postions of 5 different exons within the DNA string
Exon_1 = [1, 352]
Exon_2 = [637, 805]
Exon_3 = [1128, 1293]
Exon_4 = [1903, 2142]
Exon_5 = [3132, 3251]

# Join all given exons into a list of lists, in order of position 1 through 5
exons = [Exon_1,Exon_2,Exon_3,Exon_4,Exon_5]
print(exons)

# Create a new empty string to store the spliced exons
spliced_gene = ""

# For each exon, extract their start and stop positions
# Given that the exon coordinates were provided in human readable form, apply knowledge of Python indices to include the starting base and excluding the ending base in each exon

for exon in exons:
    start = int(exon[0]-1)
    stop = int(exon[1])
    print(start)
    print(stop)
    
# For each exon, extract and store the new Python notation coordinates from the DNA string
    exon = my_gene[start:stop]
    
# Append each exon to the spliced gene string
    spliced_gene += exon
print(spliced_gene)
```

    [[1, 352], [637, 805], [1128, 1293], [1903, 2142], [3132, 3251]]
    0
    352
    636
    805
    1127
    1293
    1902
    2142
    3131
    3251
    ATGAGAAAAATGTCGGAGGAAGAGTTTTATCTGTTCAAAAATATCTCTTCAGTGGGGCCGTGGGATGGGCCTCAGTACCACATTGCCCCTGTCTGGGCCTTCTACCTCCAGGCAGCTTTCATGGGCACTGTCTTCCTTATAGGGTTCCCACTCAATGCCATGGTGCTGGTGGCCACACTGCGCTACAAAAAGTTGCGGCAGCCCCTCAACTACATTCTGGTCAACGTGTCCTTCGGAGGCTTCCTCCTCTGCATCTTCTCTGTCTTCCCTGTCTTCGTCGCCAGCTGTAACGGATACTTCGTCTTCGGTCGCCATGTTTGTGCTTTGGAGGGCTTCCTGGGCACTGTAGCAGGTCTGGTTACAGGATGGTCACTGGCCTTCCTGGCCTTTGAGCGCTACATTGTCATCTGTAAGCCCTTCGGCAACTTCCGCTTCAGCTCCAAGCATGCACTGACGGTGGTCCTGGCTACCTGGACCATTGGTATTGGCGTCTCCATCCCACCCTTCTTTGGCTGGAGCCGGTTCATCCCTGAGGGCCTGCAGTGTTCCTGTGGCCCTGACTGGTACACCGTGGGCACCAAATACCGCAGCGAGTCCTATACGTGGTTCCTCTTCATCTTCTGCTTCATTGTGCCTCTCTCCCTCATCTGCTTCTCCTACACTCAGCTGCTGAGGGCCCTGAAAGCTGTTGCAGCTCAGCAGCAGGAGTCAGCTACGACCCAGAAGGCTGAACGGGAGGTGAGCCGCATGGTGGTTGTGATGGTAGGATCCTTCTGTGTCTGCTACGTGCCCTACGCGGCCTTCGCCATGTACATGGTCAACAACCGTAACCATGGGCTGGACTTACGGCTTGTCACCATTCCTTCATTCTTCTCCAAGAGTGCTTGCATCTACAATCCCATCATCTACTGCTTCATGAATAAGCAGTTCCAAGCTTGCATCATGAAGATGGTGTGTGGGAAGGCCATGACAGATGAATCCGACACATGCAGCTCCCAGAAAACAGAAGTTTCTACTGTCTCGTCTACCCAAGTTGGCCCCAACTGA

