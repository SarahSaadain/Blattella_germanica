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
extract only COI-region KZ616132.1:4294-5826 using extract_region_from_bam.py in folder ~/scripts/ with COI_region.bed

**2, make consensus file**  
first tried with bcftools with script test_make_consensus.py in folder ~/scripts/trials, this instead created a variant-file  

used ANGSD instead (publication: https://www.mdpi.com/1422-0067/23/9/4651)
using the script create_consensus_ANGSD.py using "-doFasta 2" and "-doCounts 1"  
this created consensus sequences for all the aDNA reads, the rest is replaced by "N"  

the output gets then mapped again (map_to_refBger2.py), sam converted to bam (convert_sam2bam.py) to view it on IGV  

**3, make multiple sequence alignment**  
prep: make sure aDNA consensus files are only the length of the COI gene and reverse complemented:
index aDNA-consensus with  
```samtools faidx 296004_S25_R1_001_aligned_sorted_COI_consensus.fasta.fa```  
cut COI region out of aDNA-consensus and reverse complement it with -i (to match the modern COI)   
```samtools faidx 296004_S25_R1_001_aligned_sorted_COI_consensus.fasta.fa KZ616132.1:4291-5826 -i > output.fasta``` 
created a script that automates this called get_COI_and_reverseComplement.py 

check the .fasta file headers to include all infos

then concatinate all .fasta files  
cat *.fasta > combined_modern_aDNA.fasta  
run mafft ```mafft --auto combined_modern_aDNA.fasta > aligned_all.fasta```

**4, plot haplotype network in R**
```
library(pegas)
library(ggplot2)

# Load and align sequences
sequences <- read.dna("~/Documents/cockroach/aDNA/mapped/merged_files_for_mafft/aligned_all.fasta", format = "fasta")
haplotypes <- haplotype(sequences)

# Construct and plot the network
haplotype_network <- haploNet(haplotypes)
plot(haplotype_network, size = attr(haplotype_network, "freq"))
```

