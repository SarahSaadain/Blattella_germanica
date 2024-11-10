Compare aDNA COI with recent COI

Idea:

Mitochondrial DNA should be more abundant, therefore look for mitochondrial gene and commenly used marker COI
map known COI sequence to ref genome
check if my aDNA mapped there too and look for differences

**1, map sequences to Bger2**  
52 modern samples containing COI sequences from Tang 2024 paper: https://www.ncbi.nlm.nih.gov/nuccore/PP692292  
Map it to Bger2.0 ref genome with map_to_refBger2.py in folder ~/scripts/  
convert sam to bam, sort and index with convert_sam2bam.py in folder ~/scripts/  

6 aDNA files (Shpak1, Shpak2, Shpak3, Dabney1, Dabney2, Dabney3) treated separetely  
also mapped, sorted and indexed to the whole Bger2 with map_to_refBger2.py in folder ~/scripts/  

**extract COI only**  
extract only COI-region KZ616132.1:4,291-5,823 using extract_region_from_bam.py in folder ~/scripts/ with COI_region.bed

**2, make consensus file**  
first tried with bcftools with script test_make_consensus.py in folder ~/scripts/trials, this instead created a variant-file  
used ANGSD instead (publication: https://www.mdpi.com/1422-0067/23/9/4651)  
this created consensus sequences for all the aDNA reads, the rest is replaced by "N"  
