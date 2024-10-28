**Estimate of the percentage of endogenous reads and their duplication rate**

- align reads to ref genome  
```GCA_000762945.1_Bger_1.0_genomic.fna
bwa mem -t 8 ref/GCA_000762945.1_Bger_1.0_genomic.fna aDNA/trimmed_aDNA/concat_trimmed_withoutLB.fastq.gz > concat_trimmed_withoutLB_aligned.sam &
samtools view -bS concat_trimmed_withoutLB_aligned.sam > concat_trimmed_withoutLB_aligned.bam
```
- filter aligned reads that have mapped to ref genome with high confidence & count aligned reads:
```
samtools view -b -q 30 concat_trimmed_withoutLB_aligned.bam > filtered_concat_trimmed_withoutLB_aligned.bam
samtools view -c filtered_concat_trimmed_withoutLB_aligned.bam
```
aligned reads: 8 169 156


- filter aligned reads & count them  
```
samtools view -c -F 4 concat_trimmed_withoutLB_aligned.bam
```
-F 4 filteres out unmapped reads  
aligned reads: 17 668 122  

- compare against total reads
with gunzip I need to devide the number later by 4
```
gunzip -c concat_trimmed_withoutLB.fastq.gz | wc -l
```

or with this unneccessarily long python script I wrote, that does the same:
```
import gzip

def count_reads_in_fastq_gz(file_path):
    with gzip.open(file_path, 'rt') as f:
        line_count = sum(1 for _ in f)
    return line_count // 4

# Update the file path to the correct location
file_path = '/Users/ssaadain/Documents/cockroach/aDNA/trimmed_aDNA/concat_trimmed_withoutLB.fastq.gz'
num_reads = count_reads_in_fastq_gz(file_path)
print(f"Number of reads: {num_reads}")
```
both times I get:  29 897 699  

Using the reads that have mapped with high confidence (-q 30)
 percentage=(aligned reads / total reads) x 100  
 percentage=(8 169 156 / 29 897 699 ) x 100  
 percentage of endogenous reads is 27.26%

Using all reads that have mapped:
 percentage=(aligned reads / total reads) x 100  
 percentage=(17 668 122 / 29 897 699 ) x 100  
 percentage of endogenous reads is 59.1%

-----
**get their duplication rate**  
-first sort .bam file  
```
samtools sort -o sorted_concat_trimmed_withoutLB_aligned.bam concat_trimmed_withoutLB_aligned.bam
```
-then mark the duplicates (-@ 8 to use 8 threads)  
```
samtools markdup -@ 8 sorted_concat_trimmed_withoutLB_aligned.bam marked_duplicates.bam
```
-then count total reads  
```
samtools view -c marked_duplicates.bam
```
-then count duplicate reads  
```
samtools view -c -f 1024 marked_duplicates.bam
```
(Number of Duplicate Reads / Number of Total Reads) ×100 = Duplication Rate
( 1 594 605 / 30 481 946 ) x100 = 5.223 %








-----
with Picard
```
picard MarkDuplicates I=raw_concat_Bger1.1_aligned_sorted.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
cat marked_dup_metrics.txt
```
inside the marked_dup_metrics.txt:  
MarkDuplicates INPUT=[raw_concat_Bger1.1_aligned_sorted.bam] OUTPUT=marked_duplicates.bam METRICS_FILE=marked_dup_metrics.txt    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000  
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000  
SORTING_COLLECTION_SIZE_RATIO=0.25  
TAG_DUPLICATE_SET_MEMBERS=false  
REMOVE_SEQUENCING_DUPLICATES=false  
TAGGING_POLICY=DontTag  
CLEAR_DT=true 
DUPLEX_UMI=false 
ADD_PG_TAG_TO_READS=true 
REMOVE_DUPLICATES=false 
ASSUME_SORTED=false 
DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES 
PROGRAM_RECORD_ID=MarkDuplicates 
PROGRAM_GROUP_NAME=MarkDuplicates 
READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 
MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 
VERBOSITY=INFO QUIET=false 
VALIDATION_STRINGENCY=STRICT 
COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 
CREATE_INDEX=false 
CREATE_MD5_FILE=false 
GA4GH_CLIENT_SECRETS=client_secrets.json 
USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
htsjdk.samtools.metrics.StringHeader

METRICS CLASS picard.sam.DuplicationMetrics
LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED SECONDARY_OR_SUPPLEMENTARY_RDS UNMAPPED_READS UNPAIRED_READ_DUPLICATES READ_PAIR_DUPLICATES READ_PAIR_OPTICAL_DUPLICATES PERCENT_DUPLICATION ESTIMATED_LIBRARY_SIZE
Unknown Library 19140038 0 636056 11765670 3320468 0 0 0,173483

-------
- quality control (Kraken, Centrifuge etc to identify and quantify non-target reads). Substract these from your total reads

------
get coverage for each region:
I already have
```
concat_trimmed_withoutLB_aligned.sam (I aligned reads to Bger_1.0)
```
convert .sam to .bam and sort it (I think its the same as sorted_concat_trimmed_withoutLB_aligned.bam)
```
samtools view -bS concat_trimmed_withoutLB_aligned.sam | samtools sort -o mapped_reads_sorted.bam
```
index
```
samtools index mapped_reads_sorted.bam
```
get coverage per position
```
samtools depth -a mapped_reads_sorted.bam > coverage.txt
```
put it in bins for easier visualization
```
bedtools makewindows -g ref/GCA_000762945.1_Bger_1.0_genomic.fna.fai -w 10000 | bedtools coverage -a - -b mapped_reads_sorted.bam > coverage_bins.txt
```
outputs 7 columns: chromosome(contig), start, end, num_reads, bases_covered, bin_size, Coverage_fraction  
continue plotting in R  

get regions with highest coverage with this bash script:  
```
# Define the input and output files
input_file="coverage_bins.txt"
output_file="high_coverage_regions.txt"

# Calculate mean coverage, sort, get the top 1% regions, and add header
{
    echo -e "chromosome\tstart\tend\tmean_coverage"  # Print header
    awk '{mean_coverage=$4/$6; print $1, $2, $3, mean_coverage}' OFS="\t" "$input_file" | \
    sort -k4,4nr | \
    head -n $(($(wc -l < "$input_file") / 100))
} > "$output_file"
```
