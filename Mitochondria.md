**Compare aDNA COI with recent COI**


Idea:  
- Mitochondrial DNA should be more abundant, therefore look for mitochondrial gene and commenly used marker COI
- map known COI sequence to ref genome  
- check if my aDNA mapped there too and look for differences   

Downloaded COI sequences from two samples from Tang 2024 paper downloaded: voucherAC2 and voucherXA1:  
https://www.ncbi.nlm.nih.gov/nuccore/PP692292  
Map it to Bger2.0 ref genome  
```bwa mem ../ref/GCA_000762945.2_Bger_2.0_genomic.fna Bger_voucherAC2_COX1.fasta > COI_voucherAC2_alignment.sam```  
Convert SAM to BAM  
```samtools view -S -b COI_voucherAC2_alignment.sam > COI_voucherAC2_alignment.bam```  
```samtools sort COI_voucherAC2_alignment.bam -o COI_voucherAC2_alignment_sorted.bam```  
```samtools index COI_voucherAC2_alignment_sorted.bam```  
Create summary statistics of the alignment  
```samtools flagstat COI_voucherAC2_alignment_sorted.bam```  
Get the location from 3rd and 4th column to find the read in IGV  
```samtools view COI_voucherAC2_alignment_sorted.bam```  
Scaffold KZ616132.1	position 4294  
Exact region in the ref genome: KZ616132.1:4.294 - KZ616132.1:5.826  
make a bedfile for this region  
``` echo -e "KZ616132.1\t4294\t5826" > region.bed``` 
filter for this region  
```samtools view -b -L region.bed mapBger2_sorted.bam > COI_aDNA_reads.bam```
index
```samtools index COI_aDNA_reads.bam```

COI: KZ616132.1:4,291-5,823

**map modern COI to ref genome and location of COI containing.bed file**
```#!/bin/bash

# Variables
REFERENCE_GENOME="/Users/ssaadain/Documents/cockroach/ref/GCA_000762945.2_Bger_2.0_genomic.fna"  # Path to the reference genome
COI_BED_FILE="../COI_region.bed"  # Path to the COI region BED file (one directory up)
OUTPUT_DIR="./aligned_bam_files"  # Directory to store BAM files

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Step 1: Align the modern FASTA files to the COI region of the reference genome using bwa mem
echo "Aligning modern FASTA files to COI region of reference genome..."
for fasta_file in *.fasta; do
    # Align the FASTA file to the COI region using the -L option to limit to the COI region
    echo "Aligning $fasta_file..."
    bwa mem -L $COI_BED_FILE $REFERENCE_GENOME $fasta_file > $OUTPUT_DIR/${fasta_file%.fasta}.sam
done

# Step 2: Convert SAM files to BAM, sort, and index them
echo "Converting SAM files to BAM, sorting, and indexing..."
for sam_file in $OUTPUT_DIR/*.sam; do
    # Convert SAM to BAM
    samtools view -Sb $sam_file > ${sam_file%.sam}.bam
    
    # Sort the BAM file
    samtools sort ${sam_file%.sam}.bam -o ${sam_file%.sam}_sorted.bam
    
    # Index the BAM file
    samtools index ${sam_file%.sam}_sorted.bam
    
    # Clean up intermediate SAM file
    rm $sam_file
done

echo "Conversion complete. BAM files are located in $OUTPUT_DIR."
```


**map each sample individually (not sure if necessary yet)**  
```
bwa mem -t 8 /Users/ssaadain/Documents/cockroach/ref/GCA_000762945.2_Bger_2.0_genomic.fna /Users/ssaadain/Documents/cockroach/aDNA/trimmed_aDNA/296004_S25_R1_001_trim.fastq.gz > /Users/ssaadain/Documents/cockroach/aDNA/mapped/296004_S25_R1_001_aligned.sam
```

**repeat filtering only for COI but for the individual aDNA files**  
get the COI region  
```
echo -e "KZ616132.1\t4294\t5826" > COI_region.bed
```
```
#!/bin/bash

# Define paths
COI_REGION_BED="COI_region.bed"  # Path to the COI region BED file
OUTPUT_DIR="./COI_filtered_bams"  # Directory to save filtered BAM files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each SAM file in the current directory
for sam_file in *.sam; do
    # Extract the base name of the SAM file (without path and extension)
    base_name=$(basename "$sam_file" .sam)
    
    # Convert the SAM file to a BAM file
    echo "Converting $sam_file to BAM..."
    bam_file="$OUTPUT_DIR/${base_name}.bam"
    samtools view -bS "$sam_file" > "$bam_file"

    # Sort the BAM file
    echo "Sorting $bam_file..."
    sorted_bam="$OUTPUT_DIR/${base_name}_sorted.bam"
    samtools sort "$bam_file" -o "$sorted_bam"

    # Filter the sorted BAM file for reads that map to the COI region
    echo "Filtering $sorted_bam for reads mapped to COI region..."
    coi_bam="$OUTPUT_DIR/${base_name}_COI.bam"
    samtools view -b -L "$COI_REGION_BED" "$sorted_bam" > "$coi_bam"

    # Index the filtered BAM file
    echo "Indexing $coi_bam..."
    samtools index "$coi_bam"

    # Clean up intermediate files (optional)
    rm "$bam_file" "$sorted_bam"

    echo "Filtered and indexed BAM file created: $coi_bam"
done

echo "Processing complete for all files."
```

**map modern .fasta files to ref genome**

```
installed ANGSD  
```

