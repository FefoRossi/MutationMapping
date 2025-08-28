# Mapping one set of reads to both wild and mutated genomes for analysis

Here the idea is to cerate a mapping that we can then recover the equivalence of mapping for mutation analysis on both the wild type genome and the mutated genome

### 1 Mapping both genomes

>Here it is important to make sure that each genome has a set of fasta headers easy to distinguish!  


```bash
cat genome_1 genome_2 > concat_genomes.fasta

bowtie2-build concat_genomes.fasta concat

bowtie2 -x concat -1 trimmed_filtered.R1.fastq -2 trimmed_filtered.R2.fastq -S mapping_output.sam -p 8 -a

samtools view -b -S mapping_output.sam > mapping_output.bam

#------------ Separating genome 1 and 2 ---------------------

samtools sort mapping_output.bam > mapping_output.sorted.bam

samtools index mapping_output.sorted.bam

samtools idxstats mapping_output.sorted.bam | cut -f 1 | gre
p '^genoma1' > genoma1_contigs.txt

samtools view -b mapping_output.sorted.bam $(cat genoma1_contigs.txt) > genoma1.bam

samtools sort genoma1.bam > genoma1.sorted.bam
```

### 2 Getting SNPs

```bash
bcftools mpileup --ff UNMAP,DUP -Ou -f genoma1.fasta genoma1.sorted.bam | bcftools call -mv -Ob -o genoma1_output.bcf

bcftools view -i 'QUAL>=20' genoma1_output.bcf > genoma1_output.filtered.bcf

bcftools view genoma1_output.filtered.bcf -o genoma1_converted.vc
```

The Python steps for the SNP analysis are the same!  
Then the provided notebooks show how to work with the generated data  