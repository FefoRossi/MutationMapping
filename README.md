# MutationMapping
Scripts and protocol for mutation mapping in prokaryotic genomes

This protocol covers the complete workflow from raw sequencing data quality control and read mapping to variant calling and mutation analysis. The pipeline includes tools for read filtering (FastQC/FastP), reference mapping (bowtie2), variant calling (BCFtools), and downstream analysis to identify missense mutations, silent mutations, and frameshift events. The final output provides comprehensive mutation mapping with amino acid substitution analysis for prokaryotic genomes.



# 1. Sequencing data
## 1.1 Read quality check and filtering

### [FastQC](https://github.com/s-andrews/FastQC)
```bash
fastqc read.R1.fastq read.R2.fastq
```
The result is a html page with all the data necessary for filtering the data!  

### [FastP](https://github.com/OpenGene/fastp)
```bash
fastp -i read.R1.fastq -I read.R2.fastq -o trimmed_filtered.R1.fastq -O trimmed_filtered.R2.fastq -e 20 -M 20 -5 -3 -W 4 -l 200
```
>Parameters:  
>**-e** Mean quality on each read  
>**-M** Mean quality over a slide-window  
>**-5** and **-3** Indicates the positions to start the slide-window (5 prime or 3 prime side)  
>**-W** Sliding-window size  
>**-l** Minimal read length to be considered after trimming (depend on the sequencing technology)  

# 2. Read mapping on reference with [bowtie2](https://github.com/BenLangmead/bowtie2)
## 2.1 Selecting reference
The reference must be the best assembled genome of the desired organism (such as RefSeq genomes). The genome must be in .fasta format and the annotation in GFF3 format.  

### Indexing the reference
```bash
bowtie2-build ref_genome.fasta reference
```
The "reference" positional argument is the index name (and prefix for all outputs). Let's keep it like this for now.  

## 2.2 Mapping the reads on the indexed reference
```bash
bowtie2 -x reference -1 trimmed_filtered.R1.fastq -2 trimmed_filtered.R2.fastq -S mapping_output.sam -p 8
```
>Parameters:  
>**-p** Number of threads  

## 2.3 Preparing bam file for next analysis ([samtools](https://www.htslib.org/))
### Converting sam to bam
```bash
samtools view -b -S mapping_output.sam > mapping_output.bam
```
### Sorting bam file
```bash
samtools sort mapping_output.bam > mapping_output.sorted.bam
```

# 3. Creating VCF from BAM [BCFtools](https://samtools.github.io/bcftools/bcftools.html)  
## 3.1 Create BCF file from reference fasta and file.bam
```bash
bcftools mpileup -Ou -f Reference.fasta bam_file | bcftools call -mv -Ob -o output.bcf
```
## 3.2 Variant Call filter by quality
```bash
bcftools view -i 'QUAL>=20' output.bcf > filtered_calls.bcf
```
## 3.3 Converting BCF to VCF
```bash
bcftools view filtered_calls.bcf -o converted.vcf
```

# 4. Starting variant analysis  

### On Jupyter notebook
```python
#Load libraries
import pandas as pd
import os
import re
from Bio import SeqIO

#Prepare codon dictionary
codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}
```

### Reading the reference GFF file

```python
gff = pd.read_table("reference.gff3", comment="#", header=None,
                    names=["contig", "type", "start", "end", "strain", "description"], usecols=[0,2,3,4,6,8])
gff = gff[gff["type"] == "CDS"]

#Usually ORF determination is within the 'locus_tag' delimitation on GFFs (please check your data before running this command)

gff["ORF"] = gff["description"].apply(lambda x: re.findall(r";locus_tag=(.+?);", x)[0])

#The same applies for the annotation, that is usually within the 'product' delimitation

gff["Annotation"] = gff["description"].apply(lambda x: re.findall(r";product=(.+?);", x)[0])
```

### Reading the VCF file

```python
vcf = pd.read_table("file.vcf", comment="#", header=None, usecols=[0,1,3,4],
                    names=["contig", "pos", "ref", "alt"])

#Let's define types of variant
def variant_type(ref, alt):
    if len(ref) > len(alt):
        return "DELETION"
    elif len(ref) < len(alt):
        return "INSERTION"
    else:
        return "SNP"
#And apply that definition in the vcf dataframe
vcf["type"] = vcf.apply(lambda x: variant_type(x["ref"], x["alt"]), axis=1)
```

### Aggregating information from GFF and VCF

```python
#1 Transforming the VCF df into a list
snps = vcf[["type", "pos", "ref", "alt"]].values.tolist()

#2 Iterating through the GFF df to add SNP information
result_data = []

for index, row in gff.iterrows():
    start = row["start"]
    end = row["end"]
    for type,snp,ref,alt in snps:
        if snp in range(start,end):
            results.append([type, snp,start,end,ref,alt,row["ORF"], row["Annotation"]])
results_df = pd.DataFrame.from_records(results, columns=["Type", "SNP_POS", "Start", "End","Ref", "Alt","ORF", "Annotation"])
```
## SNP/SNV analysis
On this step we will compare each affected codon to the reference to check if there is any mutation associated to aminoacid substitution.  

```python
from Bio import SeqIO

# Load sequences only once
ref_seq = str(next(SeqIO.parse("reference.fasta", "fasta")).seq)

for index, row in snps_only.iterrows():
    snp = row["SNP_POS"]
    start = row["Start"]
    end = row["End"]
    ref = row["Ref"]
    alt = row["Alt"]

    # Find codon position that have the SNP
    for codon_start in range(start, end, 3):
        codon_end = codon_start + 3
        if snp >= codon_start and snp < codon_end:
            codon = ref_seq[codon_start-1:codon_end-1]
            # Position of the SNP inside of the Codon
            snp_in_codon = (snp - codon_start)
            # Replace the correct base
            alt_codon = codon[:snp_in_codon] + alt + codon[snp_in_codon+1:]

            snps_only.loc[index, "original codon"] = codon
            snps_only.loc[index, "original AA"] = codontab.get(codon, "X")
            snps_only.loc[index, "mutated codon"] = alt_codon
            snps_only.loc[index, "mutated AA"] = codontab.get(alt_codon, "X")
            break  # Codon found, exit the loop

# Define mutation types
def mutation(before, after):
    if after == "*":
        return "Nonsense"
    elif before != after:
        return "Missense"
    else:
        return "Silence"

snps_only["Mutation"] = snps_only.apply(lambda x: mutation(x["original AA"], x["mutated AA"]), axis=1)
```

After this analysis we can observe if any detected SNP can directly alter aminoacid composition in the subject protein.  

### Summarising all missense (mutations with aminoacid substitution) mutations observed

```python
#Filter all missense mutations
missense = snps_only[snps_only["Mutation"] == "Missense"]

SNP_summary_table = missense.groupby(["ORF", "Annotation"]).size().reset_index(name="Occurrences")
SNP_summary_table = SNP_summary_table.sort_values(by="Occurrences", ascending=False)
```

### Checking for frame shifts

```python
insertions_deletions = results_df[(results_df["Type"] == "INSERTION") | (results_df["Type"] == "DELETION")]

#Creating a definition for the identification of frame shifts
def find_frameshift(ref, alt):
    insertion_size = len(alt) - len(ref)
    if insertion_size%3 != 0:
        return "FRAMESHIFT"
    if insertion_size%3 == 0:
        number_of_codons = insertion_size/3
        return f"INSERTION OF {number_of_codons} CODONS"

insertions_deletions["FrameShift"] = insertions_deletions.apply(lambda x: find_frameshift(x["Ref"],x["Alt"]), axis=1)
```

