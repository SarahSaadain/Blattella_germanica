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
compare fastqc outputs (without 296010, which is the library blank)
```
multiqc . -f -i "Cockroach Samples FastQC Report" --ignore "*I1*" --ignore "*I2*" --ignore "*296010*"
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
**ART_Illumina**
transform ref genome to fastq so I can use Kraken to see what contaminated the Bgerm ref genome ART_Illumina simulates Illumina reads from ref genome data (used it on Roco)

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

upload simulated reads to vetlinux01 and run Kraken2 to see what contaminated the ref genome

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

tried the same with the seven-spotted ladybug taxon Id (as a random check how often other species occur)  
```
awk '$1 == "C" && $3 == "41139"' Kraken2Bgermanica_aDNA.txt | wc -l
```  
did not occur once

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

I redid Kraken2 with concat_trimmed_withoutLB.fastq.gz to produce Kraken2Bger_withoutLB.txt  
(in this I removed the library blank file: 296010_S31_R1_001.fastq.gz):  
29897699 sequences (1257.47 Mbp) processed in 465.368s (3854.7 Kseq/m, 162.13 Mbp/m).  
  903641 sequences classified (3.02%)  
  28994058 sequences unclassified (96.98%)  
the 5 most common taxons were the same as in the runs with the LB included  

I did Kraken2 with the ref genome Bger2.0 that I split into simulated reads with Art_illumina:  
```
kraken2 --db /Volumes/Temp2/KrakenDB/nt \  
        --threads 10 \  
        --paired /Volumes/Temp2/ssaadain/Art_Illumina/Bger2.0/simulated_reads1.fq /Volumes/Temp2/ssaadain/Art_Illumina/Bger2.0/simulated_reads2.fq \  
        --output Kraken2Bger2.0_ref.txt
```  
172998887 sequences (51899.67 Mbp) processed in 20628.859s (503.2 Kseq/m, 150.95 Mbp/m).  
  35687906 sequences classified (20.63%)  
  137310981 sequences unclassified (79.37%)  
and those are the 5 most common hits:  
10282418 6973 - Blattella germanica  
1690615 2759 - cellular organisms ?  
1167144 131567 - cellular organisms ?  
1110270 33213 - Bilateria ?  
779169 117571 - Euteleostomi (bony vertebrates) ?  

same with Bger1.0 ref genome  
160581312 sequences (48174.39 Mbp) processed in 1060.463s (9085.5 Kseq/m, 2725.66 Mbp/m).  
  31826567 sequences classified (19.82%)  
  128754745 sequences unclassified (80.18%)  
8473589 6973 - Blattella germanica  
1542698 2759 cellular organisms ?  
1092662 131567 - cellular organisms ? 
1000409 33213 - Bilateria ?  
713602 117571 - Euteleostomi (bony vertebrates) ? 

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

**CREATE THE STANDARD KRAKEN2 DATABASE**  

```
kraken2-build --standard --threads 50 --db KrakenDB_2024
```
-------

**CENTRIFUGE**
checking for contamination with Centrifuge as it should avoid the kmer problem
```
conda activate SS_Centrifuge
nohup centrifuge -x /Volumes/Temp2/Centrifuge/nt -U /Volumes/Temp2/ssaadain/concat_trimmed_withoutLB.fastq.gz -S centrifuge_output.txt --report-file centrifuge_report.tsv -p 10 --verbose --seed 999 & disown
```
get only hits where centrifuge found a match (could classify the reads):  
```
awk '$3 != 0' centrifuge_output.txt > centrifuge_hits.txt
```
look for specific taxon  
```
awk '$3 == 9606' centrifuge_output.txt > human_hits.txt
awk '$3 == 6973' centrifuge_output.txt > cockroach_hits.txt
```
count hits:  
```
wc -l cockroach_hits.txt
wc -l human_hits.txt
```
796083 cockroach_hits.txt  
142314 human_hits.txt  

get hits above a certain threshold:  
```
awk '$4 > 150' centrifuge_output.txt > high_score_hits.txt
```
count hits per TaxID:  
```
awk '{print $3}' centrifuge_hits.txt | sort | uniq -c | sort -nr > taxon_counts.txt
```
Best hits:  
796083 6973 - Blattella germanica  
169081 215358 - Larimichthys crocea (some weird fish)  
142314 9606 - Homo sapiens  
98697 331104 - Blattabacterium sp. (Blattella germanica) str. Bge  
86691 10090 - Mus musculus  

some centrifuge specific tools that needs some indexing before but I didn't do it yet:
```
centrifuge-kreport -x <index> centrifuge_output.txt > centrifuge_report.txt
```
later use Pavian to parse Centrifuge (I can also use it for Kraken2 and compare it)

--------
**get number of endogenous reads**

- align reads to ref genome
- filter aligned reads that have mapped to ref genome with high confidence
- count the aligned reads
- compare against total reads
- quality control (Kraken, Centrifuge etc to identify and quantify non-target reads). Substract these from your total reads
 
```
samtools view -f 4 raw_concat_Bger1.1_aligned_sorted.bam | wc -l
samtools view -c raw_concat_Bger1.1_aligned_sorted.bam
```
number of unmapped reads: 11765670  
number of mapped reads: 19776094  
total reads: 31541764  

Proportion of endogenous reads= Number of mapped reads / Total reads x 100  
​
19,776,094/31,541,764​ × 100 ≈ 62.71%  

or another option is
```
samtools flagstat raw_concat_Bger1.1_aligned_sorted.bam
```
output:  
31541764 + 0 in total (QC-passed reads + QC-failed reads)  
30905708 + 0 primary  
0 + 0 secondary  
636056 + 0 supplementary  
0 + 0 duplicates  
0 + 0 primary duplicates  
19776094 + 0 mapped (62.70% : N/A)  
19140038 + 0 primary mapped (61.93% : N/A)  
0 + 0 paired in sequencing  
0 + 0 read1  
0 + 0 read2  
0 + 0 properly paired (N/A : N/A)  
0 + 0 with itself and mate mapped  
0 + 0 singletons (N/A : N/A)  
0 + 0 with mate mapped to a different chr  
0 + 0 with mate mapped to a different chr (mapQ>=5)  

with Picard
```
picard MarkDuplicates I=raw_concat_Bger1.1_aligned_sorted.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
cat marked_dup_metrics.txt
```
inside the marked_dup_metrics.txt:  
MarkDuplicates INPUT=[raw_concat_Bger1.1_aligned_sorted.bam] OUTPUT=marked_duplicates.bam METRICS_FILE=marked_dup_metrics.txt    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000  
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000  
SORTING_COLLECTION_SIZE_RATIO=0.25  
TAG_DUPLICATE_SET_MEMBERS=false  
REMOVE_SEQUENCING_DUPLICATES=false  
TAGGING_POLICY=DontTag  
CLEAR_DT=true 
DUPLEX_UMI=false 
ADD_PG_TAG_TO_READS=true 
REMOVE_DUPLICATES=false 
ASSUME_SORTED=false 
DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES 
PROGRAM_RECORD_ID=MarkDuplicates 
PROGRAM_GROUP_NAME=MarkDuplicates 
READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 
MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 
VERBOSITY=INFO QUIET=false 
VALIDATION_STRINGENCY=STRICT 
COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 
CREATE_INDEX=false 
CREATE_MD5_FILE=false 
GA4GH_CLIENT_SECRETS=client_secrets.json 
USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
htsjdk.samtools.metrics.StringHeader

METRICS CLASS picard.sam.DuplicationMetrics
LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED SECONDARY_OR_SUPPLEMENTARY_RDS UNMAPPED_READS UNPAIRED_READ_DUPLICATES READ_PAIR_DUPLICATES READ_PAIR_OPTICAL_DUPLICATES PERCENT_DUPLICATION ESTIMATED_LIBRARY_SIZE
Unknown Library 19140038 0 636056 11765670 3320468 0 0 0,173483

-------

Kraken
368 741 germanica
46 668 bacterium
19 899 humans

Centrifuge
796 083 germanica
169 081 bacterium
142 314 human





