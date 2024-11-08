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
