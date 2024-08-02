FINDING TEs IN REF GENOME  
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
LIBRARY PREP
aDNA from 3 (150y old) samples, two different libary preps used:  
296004: Shapk (SC1)  
296005: Shapk (SC2)  
296006: Shapk (SC3)  
296007: Dabney (DC1)  
296008: Dabney (DC2)  
296009: Dabney (DC3)  
296010: Library blank (LB)  

ADAPTER
Illumina pdf: https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2019/03/illumina-adapter-sequences-2019-1000000002694-10.pdf  
its the adapter used for TruSeq Single Indexes  
its TruSeq Universal Adapter (P7 Adapter): AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 

TRIMMING
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


REDO TRIMMING (more suitable for aDNA)
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

QUALITY CONTROL
```
fastqc -o /home/vetlinux04/Sarah/cockroach/aDNA/fastqc_trimmed /home/vetlinux04/Sarah/cockroach/aDNA/trimmed/*_R1_001_trim.fastq.gz
```

----------------
Installation of Kraken2 on vetlinux01
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
run Kraken2:
```
kraken2 --db /Volumes/Temp2/KrakenDB/nt --threads 10 --gzip-compressed Blattella_germanica/universal_trimmed.fastq.gz --output Kraken2Bgermanica
```

output:
27280436 sequences (1221.42 Mbp) processed in 43.650s (37499.2 Kseq/m, 1678.94 Mbp/m).  
  905467 sequences classified (3.32%)  
  26374969 sequences unclassified (96.68%)
