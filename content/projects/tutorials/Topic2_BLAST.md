
# Solutions to Bioinformatics Python Problems
## Topic 2: Working with BLAST

## 1. Counting amino acids

#### Given: A protein string *s*
#### Required: The number of amino acids in *s*


```python
# Read in protein string
file = open("my_protein.txt")
my_protein = file.read().rstrip()

# Use the count() method to return the number of amino acids in the string
AA_count = len(my_protein)
print(AA_count)
```

    348


## 2. Extracting sequence names from a FASTA file

#### Given: A protein FASTA file *f*
#### Required: Store the protein sequence name from *f* as a string.


```python
# Read in protein FASTA file
file = open('my_protein.fasta')
my_protein = file.readline()

# Use the strip() method to remove line breaks from the end of the header line
# Use Python indices to remove the ">" character from the header line
protein_name = my_protein[1:].strip()
print(protein_name)
```

    lcl|NG_011705.1_prot_NP_000541.1_1 [gene=TYRP1] [protein=5,6-dihydroxyindole-2-carboxylic acid oxidase precursor] [protein_id=NP_000541.1]


## 3. Extracting Accession Identifiers for all hits in a BLASTp output

#### Given: A BLASTp output file *b*
#### Required: Extract all accession identifiers from *b* and store them in a list.

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

    ['NP_000541.1', 'AKI71692.1', 'AKI71691.1', 'AKI71690.1', 'PNI54501.1', 'XP_004047853.3', 'XP_008963613.1', 'XP_016816075.1', 'EAW58712.1', 'XP_003260485.1', 'XP_031993856.1', 'XP_002819837.1', 'XP_007967434.1', 'CAA35785.1', 'XP_021783195.1', 'XP_025215330.1', 'CAG28611.1', 'XP_011785799.1', 'XP_033093247.1', 'XP_026309656.1', 'XP_011823293.1', 'XP_028691005.1', 'XP_010381456.1', 'XP_011911815.1', 'XP_015292654.1', 'XP_017722683.1', 'XP_011752097.1', 'XP_012308953.1', 'XP_032156115.1', 'XP_003928210.1', 'XP_006163908.1', 'NP_001284424.1', 'XP_014647902.1', 'XP_008063510.1', 'XP_005328038.1', 'XP_026237574.1', 'XP_012625175.1', 'XP_030733024.1', 'XP_010855297.1', 'XP_015354587.1', 'XP_015327928.2', 'XP_027405055.1', 'XP_025784242.1', 'XP_030150001.1', 'XP_019821514.1', 'XP_004313797.3', 'XP_015327929.2', 'XP_016079159.1', 'XP_011229914.1', 'XP_026961659.1', 'DAA26938.1', 'AHZ45543.1', 'VFV18838.1', 'XP_029776228.1', 'XP_004276181.1', 'XP_014929773.1', 'XP_024414130.1', 'KAF0871959.1', 'AHZ45544.1', 'AAK85404.1', 'XP_005891079.1', 'XP_006918424.1', 'XP_003782933.1', 'NP_776905.2', 'NP_001036025.2', 'XP_007098614.1', 'XP_004706584.1', 'TEA24065.1', 'AAY87040.1', 'XP_011380526.1', 'XP_025784165.1', 'XP_019271058.1', 'XP_013849023.1', 'XP_025784089.1', 'MXQ93697.1', 'AAV65117.1', 'ACS66870.1', 'AAV65119.1', 'XP_003407396.1', 'VFV18836.1', 'XP_029776227.1', 'ADB96142.1', 'VFV18837.1', 'XP_029776226.1', 'XP_026361788.1', 'ELR59806.1', 'ALI89102.1', 'XP_012977523.1', 'XP_008688074.1', 'NP_001020397.1', 'XP_030619664.1', 'XP_025857092.1', 'XP_006208354.1', 'XP_027471880.1', 'XP_004765815.1', 'XP_010986480.1', 'ADB96143.1', 'XP_029058587.1', 'ALV13267.1', 'XP_007465698.1']


## 4. Calculating the average percent of bases that differ between two species

#### Given: A BLASTn output file *b*
#### Required: Store the percent of sites that differ between Reference 1 and Reference 2 mRNAs from *b* as a float (do not include a percent symbol).

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

    8.571173946817837


## 5. Quantify the ratio of indels to substitutions between two species

#### Given: A BLASTn output file *b*
#### Required: Store the ratio of indels to single nucleotide substitions between Reference 1 and Reference 2 mRNAs from *b* as a float (do not include a percent symbol).


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

    0.1154652840199936

