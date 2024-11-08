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
Extract aDNA reads that map to the COI gene with a small script:  
```# Step 1: Create a temporary file with the region coordinates in proper BED format
echo -e "KZ616132.1\t4294000\t5826000" > region.bed

# Step 2: Use samtools to filter the reads in the specified region
samtools view -b -L region.bed mapBger2_sorted.bam > COI_aDNA_reads.bam

# Step 3: Clean up the temporary region file (optional)
rm region.bed```



