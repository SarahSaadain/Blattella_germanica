Compare aDNA COI with recent COI

Idea:

Mitochondrial DNA should be more abundant, therefore look for mitochondrial gene and commenly used marker COI
map known COI sequence to ref genome
check if my aDNA mapped there too and look for differences

52 modern samples containing COI sequences from Tang 2024 paper: https://www.ncbi.nlm.nih.gov/nuccore/PP692292  
Map it to Bger2.0 ref genome  
```bwa mem ../ref/GCA_000762945.2_Bger_2.0_genomic.fna Bger_voucherAC2_COX1.fasta > COI_voucherAC2_alignment.sam```  
or use map_to_refBger2.py in scripts folder ~/scripts/  
convert sam to bam, sort and index with convert_sam2bam.py in scripts folder ~/scripts/  

6 aDNA files (Shpak1, Shpak2, Shpak3, Dabney1, Dabney2, Dabney3) treated separetely  
also mapped, sorted and indexed to the whole Bger2 with map_to_refBger2.py in scripts folder ~/scripts/  
map only to COI-region KZ616132.1:4,291-5,823 by using script extract_region_from_bam.py in folder ~/scripts/ with COI_region.bed


