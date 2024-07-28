got genomes to Roco with
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/762/945/GCA_000762945.1_Bger_1.0/GCA_000762945.1_Bger_1.0_genomic.fna.gz
```
open docker imagine and run earlgrey on docker:
```
docker ps
docker exec -ti containername (the cointerID) sh
earlGrey -g GCA_000762945.2_Bger_2.0_genomic.fna -s B.germanica -o earlGrey/ -c yes -t
```

-----
trimming script:
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


