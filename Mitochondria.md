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

**map each sample individually (not sure if necessary yet)**  
```
bwa mem -t 8 /Users/ssaadain/Documents/cockroach/ref/GCA_000762945.2_Bger_2.0_genomic.fna /Users/ssaadain/Documents/cockroach/aDNA/trimmed_aDNA/296004_S25_R1_001_trim.fastq.gz > /Users/ssaadain/Documents/cockroach/aDNA/mapped/296004_S25_R1_001_aligned.sam
```

**repeat filtering only for COI but for the individual aDNA files**  
get the COI region  
```
echo -e "KZ616132.1\t4294000\t5826000" > COI_region.bed
```
```
#!/bin/bash

# Define paths
COI_REGION_BED="COI_region.bed"  # Path to the COI region BED file
OUTPUT_DIR="./COI_filtered_bams"  # Directory to save filtered BAM files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each SAM file in the current directory
for sam_file in *_aligned.sam; do
    # Extract the base name of the SAM file (without path and extension)
    base_name=$(basename "$sam_file" _aligned.sam)
    
    # Define the output BAM file name
    output_bam="$OUTPUT_DIR/${base_name}_aligned.bam"
    
    # Convert SAM to BAM
    echo "Converting $sam_file to BAM..."
    samtools view -Sb "$sam_file" > "$output_bam"
    
    # Filter the BAM file for reads that map to the COI region
    filtered_bam="$OUTPUT_DIR/${base_name}_COI.bam"
    echo "Filtering $output_bam for reads mapped to COI region..."
    samtools view -b -L "$COI_REGION_BED" "$output_bam" > "$filtered_bam"
    
    # Index the filtered BAM file
    samtools index "$filtered_bam"
    
    echo "Filtered BAM file created: $filtered_bam"
done

echo "Processing complete for all files."
```








