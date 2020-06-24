---
title: Working with BLAST
linktitle: Topic 2
toc: true
type: docs
date: "2019-05-05T00:00:00+01:00"
draft: false
menu:
  example:
    parent: Python Solutions
    weight: 1

# Prev/next pager order (if `docs_section_pager` enabled in `params.toml`)
weight: 1
---

## 1. Counting amino acids

Given: A protein string *s*.

Required: The number of amino acids in *s*

```python
# Read in protein string
file = open("my_protein.txt")
my_protein = file.read().rstrip()

# Use the count() method to return the number of amino acids in the string
AA_count = len(my_protein)
print(AA_count)
```

## 2. Extracting sequence names from a FASTA file

Given: A protein FASTA file *f*.

Required: Store the protein sequence name from *f* as a string.

```python
# Read in protein FASTA file
file = open('my_protein.fasta')
my_protein = file.readline()

# Use the strip() method to remove line breaks from the end of the header line
# Use Python indices to remove the ">" character from the header line
protein_name = my_protein[1:].strip()
print(protein_name)
```

## 3. Extracting Accession Identifiers for all hits in a BLASTp output

Given: A BLASTp output file *b*.

Required: Extract all accession identifiers from *b* and store them in a list.

### *OPTIONAL STEP*
#### Running BLASTp through a command line call (BASH)

```bash
%%bash
blastp -query my_protein.fasta -db NR -out BLASTp_output.txt -outfmt 6
```

### Where:
    -query is the protein fasta file input
    -db is the protein database used to perform the search (Non Redundant)
    -out is the BLASTp output stored as a txt file
    -outfmt is the flag that determines which BLAST output format to use (6 will output a table)

```python
# Read in BLASTp output
blast_output = open("BLASTp_output.txt")
blast = blast_output.readlines()

# Create an empty list to store accessions
accessions = []

# Iterate over each line of the blast output
# Use the split() method to split each line on the "	" character
for line in blast:
    splits = line.split("	")

# Append the GI encoded in the subject sequence for each line of the blast output to the empty list
    accessions.append(splits[1])
print(accessions)
```

## 4. Calculating the average percent of bases that differ between two species

Given: A BLASTn output file *b*.

Required: Store the percent of sites that differ between Reference 1 and Reference 2 mRNAs from *b* as a float (do not include a percent symbol).

### *Optional Step*
#### Running BLASTn through a command line call (BASH)

```bash
%%bash
#makeblastdb -in ref_2_mRNA.fasta -dbtype nucl -parse_seqids
blastn -query ref_1_mRNA.fasta -db ref_2_mRNA_fix.fasta -out ref1_vs_ref2.mRNA.txt -outfmt 6 -num_alignments 1 -max_hsps 1 -evalue 1e-10
```

    Warning: [blastn] Examining 5 or more matches is recommended



### Additional flags:
    -num_alignments tells BLAST to only report the 1 best alignment for each query
    -max_hsps tells BLAST to only report the top high scoring pair (HSP), or the best local alignment
    -evalue tells BLAST to only output alignments with e-values less than the assigned number

```python
# Read in BLASTn output
blast_output = open("ref1_vs_ref2.mRNA.txt")
blast = blast_output.readlines()

# Create an empty list to store the percent of similar residues between references 1 and 2
percent_similar = []

# Iterate over each line in the blast output
# Use the split() method to split each line on the space character
for line in blast:
    b = line.split()

# Convert each percent of identical residues in the alignment from a number to a float
    c = float(b[2])
    
# Append each percent of identical residues to the empty list
    percent_similar.append(c)

# To find the percent of different residues between the references:
    # Use the sum() function to find the total sum of all identical residues
    # Use the len() function to find the length of the identical residues list
    # Divide the sum of the identical residues by the length of the list
    # Subtract the previous value from 100 to find the percent of residues that differ
    # Ensure that the final result is stored to a new variable as a float
percent_different = 100-(sum(percent_similar)/len(percent_similar))
print(percent_different)
```

## 5. Quantify the ratio of indels to substitutions between two species

Given: A BLASTn output file *b*.

Required: Store the ratio of indels to single nucleotide substitions between Reference 1 and Reference 2 mRNAs from *b* as a float (do not include a percent symbol).

```python
# Read in BLASTn output
blast_output = open("ref1_vs_ref2.mRNA.txt")
blast = blast_output.readlines()

# Create two empty lists to store the values of substitutions and indels 
subs = []
indels = []

# Iterate over each line in the blast output
# Use the split() method to split each line on the space character
for line in blast:
    a = line.split()

# Convert the number of mismatches and the number of gaps in each alignment to a float
    b = float(a[4])
    c = float(a[5])

# Append the number of mismatches and substitutions to their appropriate lists
    subs.append(b)
    indels.append(c)

# To calculate the ratio of indels to substitutions between the references:
    # Use the sum() function to calculate the total number of indels and substitutions in each list
    # Divide the total number of indels by the total number of substitutions
    # Ensure that the final result is stored to a new variable as a float
indels_to_substitions = float(sum(indels)/sum(subs))
print(indels_to_substitions)
```
