**FINDING TEs IN REF GENOME**  
Earlgrey with publically available data:

got ref genome from NCBI to Roco with
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/762/945/GCA_000762945.1_Bger_1.0/GCA_000762945.1_Bger_1.0_genomic.fna.gz
```
open docker imagine and run earlgrey on docker:
```
docker ps
docker exec -ti containername (the cointerID) sh
earlGrey -g GCA_000762945.2_Bger_2.0_genomic.fna -s B.germanica -o earlGrey/ -c yes -t
```
did not work because ref genome is not of good quality

-----
**LIBRARY PREP**  
aDNA from 3 (150y old) samples, two different libary preps used:  
296004: Shapk (SC1)  
296005: Shapk (SC2)  
296006: Shapk (SC3)  
296007: Dabney (DC1)  
296008: Dabney (DC2)  
296009: Dabney (DC3)  
296010: Library blank (LB)  


concatinate all files for Kraken2 test
```
cat *.fastq.gz > concatenated.fastq.gz
```

concatinate files without LB (library blank) for Kraken2:  
concat_raw_withoutLB.fastq.gz (untrimmed, without LB)  
concat_trimmed_withoutLB.fastq.gz (with less stringent settings in cutadapt, without LB)  
universal_concat_trimmed_withoutLB.fastq.gz (the normal cutadapt settings - more stringent, without LB)

-----
**ADAPTER**  
Illumina pdf: https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2019/03/illumina-adapter-sequences-2019-1000000002694-10.pdf  
its the adapter used for TruSeq Single Indexes  
its TruSeq Universal Adapter (P7 Adapter): AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 

-----
**TRIMMING**
```
#!/bin/bash

# Define input and output directories
INPUT_DIR="/home/vetlinux04/Sarah/cockroach/aDNA/raw"
OUTPUT_DIR="/home/vetlinux04/Sarah/cockroach/aDNA/trimmed"

# Adapter sequences
ADAPTER_A="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

# Loop through each FASTQ file in the input directory
for FILE in $INPUT_DIR/*_R1_001.fastq.gz; do
  # Extract the base name of the file
  BASENAME=$(basename $FILE)
  
  # Define the output file name
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME%.fastq.gz}_trim.fastq.gz"
  
  # Run cutadapt
  cutadapt -a $ADAPTER_A -e 0.1 -O 5 -m 20 -q 20 -o $OUTPUT_FILE $FILE
done
```
-e 0.1: Allows up to 10% errors in the adapter sequence.  
-O 5: Requires at least 5 base pairs of overlap between the adapter and the read.  
-m 20: Discards reads shorter than 20 bases after trimming.  
-q 20: Trims low-quality bases from the 3' end before adapter removal.

this increased the quality from 12 to 24 but raw data had a mean read length of 98nt, after trimming the mean read length is now 30nt


**REDO TRIMMING (more suitable for aDNA)**
```
#!/bin/bash

# Define input and output directories
INPUT_DIR="/home/vetlinux04/Sarah/cockroach/aDNA/raw"
OUTPUT_DIR="/home/vetlinux04/Sarah/cockroach/aDNA/trimmed"

# Adapter sequences
ADAPTER_A="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

# Loop through each FASTQ file in the input directory
for FILE in $INPUT_DIR/*_R1_001.fastq.gz; do
  # Extract the base name of the file
  BASENAME=$(basename $FILE)
  
  # Define the output file name
  OUTPUT_FILE="$OUTPUT_DIR/${BASENAME%.fastq.gz}_trim.fastq.gz"
  
  # Run cutadapt
  cutadapt -a $ADAPTER_A -e 0.1 -O 5 -m 1 -q 5 -o $OUTPUT_FILE $FILE
done
```
----
**QUALITY CONTROL**
```
fastqc -o /home/vetlinux04/Sarah/cockroach/aDNA/fastqc_trimmed /home/vetlinux04/Sarah/cockroach/aDNA/trimmed/*_R1_001_trim.fastq.gz
```

----------------
**INSTALLATION** of Kraken2 on vetlinux01
```
conda init (to initialize conda)
conda config --add channels defaults (add the default channels)
conda config --add channels bioconda (add the bioconda channel)
conda config --add channels conda-forge (add the conda-forge channel)
conda config --set channel_priority strict (set channel priority to strict)
conda create -n SS_Kraken2 kraken2 bracken (install Kraken2 and Bracken)
```

to open environment:
```
conda activate SS_Kraken2
```
it activated the right environment but didn't show that I am in the right environment (I double checked with "conda info"). Therefore I configured the promt by doing:

```
nano ~/.zshrc
```
```
# Set the prompt to show the active conda environment
function set_conda_prompt {
    if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
        PS1="($CONDA_DEFAULT_ENV) $PS1"
    fi
}
precmd_functions+=(set_conda_prompt)
```
then I used
```
source ~/.zshrc
```
----------
**KRAKEN2**
```
kraken2 --db /Volumes/Temp2/KrakenDB/nt --threads 10 --gzip-compressed Blattella_germanica/universal_trimmed.fastq.gz --output Kraken2Bgermanica
```

output using the universal_trimmed.fastq.gz (stringend cutadapt settings): Kraken2Bgermanica.txt  
27280436 sequences (1221.42 Mbp) processed in 43.650s (37499.2 Kseq/m, 1678.94 Mbp/m).  
  905467 sequences classified (3.32%)  
  26374969 sequences unclassified (96.68%)
  
output using the concatenated_samples.fastq.gz (less stringent cutadapt settings): Kraken2Bgermanica_aDNA.txt  
30015710 sequences (1268.15 Mbp) processed in 43.569s (41335.3 Kseq/m, 1746.40 Mbp/m).  
  906698 sequences classified (3.02%)  
  29109012 sequences unclassified (96.98%)

output using the raw_concatenated.fastq.gz (the raw sequences without adapters removed): Kraken2Bgermanica_raw.txt  
30905708 sequences (3090.57 Mbp) processed in 186.839s (9924.8 Kseq/m, 992.48 Mbp/m).
  26787768 sequences classified (86.68%)
  4117940 sequences unclassified (13.32%)


get the lines that are not U (unclassified)
```
grep -v '^U' Kraken2Bgermanica_aDNA.txt > filtered_output.txt
```

```
awk '$1 == "C" {print $3}' Kraken2Bgermanica_aDNA.txt | sort -u > taxon_ids.txt
```

looking only for the B. germanica taxon id:
```
awk '$1 == "C" && $3 == "6973"' Kraken2Bgermanica_aDNA.txt > filtered_6973.txt
```
```
awk '$1 == "C" && $3 == "6973"' Kraken2Bgermanica_aDNA.txt | wc -l
````
finds it 368741 times  
double checking:
```
grep '6973' filtered_6973.txt | wc -l
````
same number

tried the same with the seven-spotted ladybug taxon Id (as a random check how often other species occur)  
```awk '$1 == "C" && $3 == "41139"' Kraken2Bgermanica_aDNA.txt | wc -l```  
did not occur once

tried the same with the human taxon Id  
```awk '$1 == "C" && $3 == "9606"' Kraken2Bgermanica_aDNA.txt | wc -l```  
found it 19899

get the taxon ids with the most hits
```
awk '$1 == "C" {print $3}' Kraken2Bgermanica_aDNA.txt | sort | uniq -c | sort -nr | head -5
```
368741 6973 - B. germanica  
  46668 331104 - Blattabacterium sp. (Blattella germanica) str. Bge   
  26002 2759 - cellular organisms ?  
  23198 543 - Enterobacterales  
  19899 9606 - Homo sapiens

get the taxon ids with the most hits in everything but U
```
awk '$1 != "U" {print $3}' Kraken2Bgermanica_raw.txt | sort | uniq -c | sort -nr | head -n 5
```
12149497 7962 - Cyprinus carpio (common carp) (see this article: https://dgg32.medium.com/carp-in-the-soil-1168818d2191)  
2367994 2759 - Eukaryota  
1968948 1 - ?  
1884534 131567 - Archae  
1716649 29278 - Vectors ?  
1377821 1306438 - Fargesia denudata (plant)  
1163784 29442 - Pseudomonas tolaasii   
 955329 2081755 - Cloning vector pCA-DEST2303  
 470272 33213 - Bilatera  
 429513 6973 - Blattella germanica  

 -----
**MAP aDNA READS TO REF GENOME**  
**indexed and aligned to ref genome**
```
bwa index GCA_000762945.1_Bger_1.0_genomic.fna  
bwa mem ref/GCA_000762945.1_Bger_1.0_genomic.fna raw_concatenated.fastq.gz > raw_concatenated_aligned.sam
```
on vetlinux01 I did:  
GCA_000762945.1_Bger_1.0_genomic.fna  

on vetlinux04 I indexed:   
GCA_003018175.1_Bger_1.1_genomic.fna  
  
on my computer I did:  
GCA_000762945.2_Bger_2.0_genomic.fna  

then
```
samtools view -Sb raw_concatenated_aligned.sam > raw_concatenated_aligned.bam
samtools sort raw_concatenated_aligned.bam -o raw_concatenated_aligned_sorted.bam
samtools index raw_concatenated_aligned_sorted.bam
```

on my computer:
index reference genome
```
samtools faidx GCA_000762945.1_Bger_1.0_genomic.fna
```

**in IGV:**  
Set the Reference Genome:  
Go to Genomes > Load Genome from File....  
Select your reference genome file (GCA_003018175.1_Bger_1.1_genomic.fna).  
Load the BAM File:  
Go to File > Load from File....  
Select your BAM file (raw_concat_Bger1.1_aligned_sorted.bam). IGV will automatically detect the index file (.bai).  
-> super scattered but widely distributed across the genome

----
**GenomeDelta**
```
bash /mnt/data2/sarah/cockrock/GD/GenomeDelta/linux/main.sh --fq universal_trimmed.fastq.gz --fa ref/GCA_000762945.1_Bger_1.0_genomic.fna --of GD --t 8
```
too low coverage, gave an Error (-1) after calculating coverage support for each gap step

-------
**ART_Illumina**  
transform ref genome to fastq so I can use Kraken to see what contaminated the Bgerm ref genome
ART_Illumina simulates Illumina reads from ref genome data
```
conda activate arcitc_env
art_illumina -ss HS25 -sam -i ref/GCA_000762945.2_Bger_2.0_genomic.fna -l 150 -f 30 -m 350 -s 50 -o simulated_reads
```
Generates paired-end reads (simulated_reads1.fq and simulated_reads2.fq) and a SAM file (simulated_reads.sam) with the read alignments to the reference genome.  
-ss HS25: Use the Illumina HiSeq 2500 error profile.  
-sam: Output a SAM file with the read alignments.  
-i ref/GCA_000762945.2_Bger_2.0_genomic.fna: Use the specified reference genome file.  
-l 150: Generate reads of length 150 base pairs.  
-f 30: Achieve 30-fold coverage of the reference genome.  
-m 350: Mean fragment size of 350 base pairs.  
-s 50: Standard deviation of 50 base pairs for fragment size.  
-o simulated_reads: Prefix for the output files (e.g., simulated_reads1.fq, simulated_reads2.fq, simulated_reads.sam, etc.).


--------
**CREATE THE STANDARD KRAKEN" DATABASE**  

```
kraken2-build --standard --threads 50 --db KrakenDB_2024
```
