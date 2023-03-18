## remapping the bulk RNAseq data to create the ref of antenna
1. fastp processing:
```## PBS configure 
#PBS -N fastp
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=16G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/bulk-RNAseq/
ls ./input/*1P.fastq.gz |cut -d "/" -f 3|cut -d "_" -f 1 >sample
cat sample |while read id;
do echo $id
arr=($id)
samplename=${arr[0]}
fastp_out=./01_fastp/
fastp -i ${samplename}_1P.fastq.gz \
        -I ${samplename}_2P.fastq.gz \
        -o ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
        -O ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
        -h ${fastp_out}${samplename}_fastp.html \
        -j ${fastp_out}${samplename}_fastp.json
       ```

2. STAR mapping 


```

## PBS configure 
#PBS -N align
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=16G
source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/bulk-RNAseq/03_STAR_remapping
STAR_index=/md01/nieyg/ref/STAR/Amel_HAv3_STAR
ls ../input/*1P.fastq.gz |cut -d "/" -f 3|cut -d "_" -f 1 >sample
cat sample |while read id;
do echo $id
arr=($id)
sample=${arr[0]}

STAR --runThreadN 12\
        --genomeDir $STAR_index \
        --readFilesIn ../01_fastp/trimmed_${sample}_R1.fastq.gz ../01_fastp/trimmed_${sample}_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix $sample \
        --outFilterType BySJout \
        --alignIntronMin 10 --alignIntronMax 1000000 --alignMatesGapMax  1000000 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
        --clip3pNbases 1  --outSAMstrandField intronMotif \
        --outSAMattrIHstart 0  --outFilterMultimapNmax 20  --outSAMattributes NH HI AS NM MD\
        --outSAMattrRGline  --outSAMtype BAM SortedByCoordinate 
```

### Need to test (mapping parameters) (changed one at a time):

* Default

```
STAR --runThreadN 12\
        --genomeDir /md01/nieyg/ref/STAR/Amel_HAv3_STAR \
        --readFilesIn ../01_fastp/trimmed_D1_R1.fastq.gz ../01_fastp/trimmed_D1_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix D1-Default \
        --outSAMtype BAM Unsorted
```
* keep only those reads that contain junctions that passed filtering into SJ.out.tab
* keep uniquely reads
* Set max and min intron length

```
STAR --runThreadN 12\
        --genomeDir /md01/nieyg/ref/STAR/Amel_HAv3_STAR \
        --readFilesIn ../01_fastp/trimmed_D1_R1.fastq.gz ../01_fastp/trimmed_D1_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix D1-filter_junction_unique_Intron_range \
        --outFilterType BySJout \
        --alignIntronMin 10 --alignIntronMax 100000 --alignMatesGapMax  100000 \
        --outSAMtype BAM Unsorted
```
* junction_gap max
```
STAR --runThreadN 12\
        --genomeDir /md01/nieyg/ref/STAR/Amel_HAv3_STAR \
        --readFilesIn ../01_fastp/trimmed_D1_R1.fastq.gz ../01_fastp/trimmed_D1_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix D1-filter_junction_unique_Intron_range_filtergap \
        --outFilterType BySJout \
        --alignIntronMin 10 --alignIntronMax 100000 --alignMatesGapMax  100000 \
        --outSJfilterIntronMaxVsReadN 10000 20000 40000 \
        --outSAMtype BAM Unsorted
```

### subtract the mapping quality high reads 
* mapping quality plot 
* subtract reads (MAPQ> 10)

```
samtools sort -@ 12 -O bam -o D1-Default_sorted.bam D1-DefaultAligned.out.bam
samtools index D1-Default_sorted.bam
python MAPQ.py --bam D1-Default_sorted.bam --png D1-Default.png
samtools view -h -b -q 10 D1-Default_sorted.bam -o D1-Default_sorted_mapqm10.bam
```
```
samtools sort -@ 12 -O bam -o D1-filter_junction_unique_Intron_range_sorted.bam D1-filter_junction_unique_Intron_rangeAligned.out.bam
samtools index D1-filter_junction_unique_Intron_range_sorted.bam
python MAPQ.py --bam D1-filter_junction_unique_Intron_range_sorted.bam --png D1-filter_junction_unique_Intron_range.png
samtools view -h -b -q 10 D1-filter_junction_unique_Intron_range_sorted.bam -o D1-filter_junction_unique_Intron_range_sorted_mapqm10.bam
```
```
samtools sort -@ 12 -O bam -o D1-filter_junction_unique_Intron_range_filtergap_sorted.bam D1-filter_junction_unique_Intron_range_filtergapAligned.out.bam
samtools index D1-filter_junction_unique_Intron_range_filtergap_sorted.bam
python MAPQ.py --bam D1-filter_junction_unique_Intron_range_filtergap_sorted.bam --png D1-filter_junction_unique_Intron_range_filtergap.png
samtools view -h -b -q 10 D1-filter_junction_unique_Intron_range_filtergap_sorted.bam -o D1-filter_junction_unique_Intron_range_filtergap_sorted_mapqm10.bam
```

### test the parameter in AaegL5.0
```
STAR --runThreadN 12\
        --genomeDir /md01/nieyg/ref/STAR/Amel_HAv3_STAR \
        --readFilesIn ../01_fastp/trimmed_D1_R1.fastq.gz ../01_fastp/trimmed_D1_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix D1-ref_AaegL \
        --outFilterType BySJout \
        --alignIntronMin 10 --alignIntronMax 1000000 --alignMatesGapMax  1000000 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
        --clip3pNbases 1  --outSAMstrandField intronMotif \
        --outSAMattrIHstart 0  --outFilterMultimapNmax 20  --outSAMattributes NH HI AS NM MD\
        --outSAMtype BAM SortedByCoordinate 
```
