# **RNA-seq Data Processing Pipeline**

*Disclaimer:* I used the counts data directly from the published source. However, these would be the procedures I would perform to process the RNA-seq data and obtain the quantified expression for each conditions of this study.

## **Step 1: Fetching Data Using SRA Toolkit**

Raw sequencing data is fetched from the **SRA database** using `prefetch` and converted to FASTQ format using `fastq-dump`. The process is automated with a for loop to handle multiple conditions.

### **Command:**
```bash
# Define sample names and their corresponding SRA IDs
declare -A samples=(
    [WT1]=SRR32039580
    [WT2]=SRR32039581
    [WT3]=SRR32039582
    [KO1]=SRR32039583
    [KO2]=SRR32039584
    [KO3]=SRR32039585
)

# Loop through each sample and fetch data
for sample in "${!samples[@]}"; do
    sra_id=${samples[$sample]}
    
    # Fetch SRA data
    prefetch $sra_id
    
    # Convert SRA to FASTQ format
    fastq-dump --split-files --gzip $sra_id -O raw_data/
    
    # Rename files for consistency
    mv raw_data/${sra_id}_1.fastq.gz raw_data/${sample}_R1.fastq.gz
    mv raw_data/${sra_id}_2.fastq.gz raw_data/${sample}_R2.fastq.gz

done
```

---

## **Step 2: Preprocessing Reads Using fastp**

To ensure high-quality reads, raw FASTQ files are cleaned using `fastp`. This step removes adapters, trims low-quality bases, and filters out short reads.

### **Command:**
```bash
mkdir -p clean_data

for sample in "${!samples[@]}"; do
    fastp -i raw_data/${sample}_R1.fastq.gz -I raw_data/${sample}_R2.fastq.gz \
          -o clean_data/${sample}_clean_R1.fastq.gz -O clean_data/${sample}_clean_R2.fastq.gz \
          -q 20 -l 50 --detect_adapter_for_pe --thread 8 --html reports/${sample}_fastp.html

done
```

### **Parameters:**
- `-q 20`: Trims bases with quality score below 20.
- `-l 50`: Discards reads shorter than 50bp.
- `--detect_adapter_for_pe`: Automatically detects and removes adapters.
- `--html reports/sample_fastp.html`: Generates a quality report for each sample.

---

## **Step 3: Aligning Reads to the Mouse Genome Using HISAT2**

The cleaned reads are aligned to the **mm10 (GRCm38)** reference genome using `HISAT2`.

### **Command:**
```bash
mkdir -p alignments

for sample in "${!samples[@]}"; do
    hisat2 -p 8 --dta -x genome_index/mm10 \
           -1 clean_data/${sample}_clean_R1.fastq.gz \
           -2 clean_data/${sample}_clean_R2.fastq.gz \
           -S alignments/${sample}_aligned.sam
    
    # Convert SAM to BAM, sort, and index
    samtools view -bS alignments/${sample}_aligned.sam | samtools sort -o alignments/${sample}_aligned_sorted.bam
    samtools index alignments/${sample}_aligned_sorted.bam

done
```

### **Parameters:**
- `--dta`: Optimizes alignment for transcript assembly.
- `-x genome_index/mm10`: Specifies the genome index for alignment.
- `-S alignments/sample_aligned.sam`: Outputs aligned reads in **SAM** format.
- `samtools view -bS | samtools sort`: Converts SAM to sorted BAM.
- `samtools index`: Indexes BAM files for downstream analysis.

---

## **Step 4: Counting Reads Using featureCounts**

Gene expression levels are quantified using `featureCounts` with a **GTF annotation file**.

### **Command:**
```bash
mkdir -p counts

featureCounts -T 8 -p -t exon -g gene_id -a annotation.gtf \
              -o counts/gene_counts.txt alignments/*_aligned_sorted.bam
```

### **Parameters:**
- `-p`: Specifies paired-end reads.
- `-t exon`: Counts only exon-mapped reads.
- `-g gene_id`: Uses gene IDs for counting.
- `-a annotation.gtf`: Provides gene annotations.
- `-o counts/gene_counts.txt`: Outputs a **tab-delimited** file with raw gene counts.

---

## **Output Files:**
1. `counts/gene_counts.txt`: Raw gene counts for all samples.
2. `alignments/sample_aligned_sorted.bam`: Sorted BAM files with aligned reads.
3. `reports/sample_fastp.html`: Quality control reports for read cleaning.

---

## **Reference Genome:**
- **Assembly:** mm10 (GRCm38)
- **Annotation File:** Ensembl GTF

---

This pipeline ensures high-quality read processing, alignment, and gene expression quantification for RNA-seq analysis.
