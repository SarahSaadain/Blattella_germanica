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


**map each sample individually**  
```
#!/bin/bash

# Define paths
ref="/Users/ssaadain/Documents/cockroach/ref/GCA_000762945.2_Bger_2.0_genomic.fna"
trimmed_dir="/Users/ssaadain/Documents/cockroach/aDNA/trimmed_aDNA"
mapped_dir="/Users/ssaadain/Documents/cockroach/aDNA/mapped"

# Loop through each trimmed file
for fastq_file in "$trimmed_dir"/*_trim.fastq.gz; do
    # Extract base name without path and extension
    base_name=$(basename "$fastq_file" _trim.fastq.gz)
    
    # Define output SAM file path
    output_sam="$mapped_dir/${base_name}_aligned.sam"
    
    # Run bwa mem with 8 threads
    bwa mem -t 8 "$ref" "$fastq_file" > "$output_sam" &
done

# Wait for all background processes to complete
wait
echo "Alignment complete for all files."
```









