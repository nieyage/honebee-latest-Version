## make ref
1. get fasta and gff3 file from NCBI (Amel_HAv3.1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Apis_mellifera/latest_assembly_versions/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz


2. change the chr name to group 

* the chromsome ID relationship 
```
# trans NCBI_ID to Group123,etc
sed -n '35,121p' GCF_003254395.2_Amel_HAv3.1_assembly_report.txt | awk '{print $1"\t"$8}' > Amel_HAv3.1_assembly_report_1.txt
sed -n '122,211p' GCF_003254395.2_Amel_HAv3.1_assembly_report.txt | awk '{print $1"\t"$7}' > Amel_HAv3.1_assembly_report_2.txt
cat Amel_HAv3.1_assembly_report_1.txt Amel_HAv3.1_assembly_report_2.txt > Amel_HAv3.1_accession_convert.txt
```
```
#!/bin/bash
nLine=`wc -l Amel_HAv3.1_accession_convert.txt | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        name=`awk '{print $1}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        accession=`awk '{print $2}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        sed -i "s/${accession}/${name}/g" GCF_003254395.2_Amel_HAv3.1_genomic.fna
done
```
```
#!/bin/bash
nLine=`wc -l Amel_HAv3.1_accession_convert.txt | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        name=`awk '{print $1}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        accession=`awk '{print $2}' Amel_HAv3.1_accession_convert.txt | sed -n "${i}p"`
        sed -i "s/${accession}/${name}/g" GCF_003254395.2_Amel_HAv3.1_genomic.gtf
done
``` 
I still haven't seen the corresponding file on the website

3. manually annotated chemoreceptors

### back to R
```
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf')
```
* get all chemoreceptor gene # get from NCBI 

```
# unique(gtf[grep("odorant receptor",gtf$product),]$gene_id) # not all 
OR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/odorant receptor gene_result.txt",sep="\t",header=T)
GR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Gustatory receptor gene_result.txt",sep="\t",header=T)
IR_info <- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/Ionotropic receptor gene_result.txt",sep="\t",header=T)
chemoreceptor <- c(OR_info$Symbol,GR_info$Symbol,IR_info$Symbol)
chemoreceptor <- chemoreceptor[which(chemoreceptor%in% gtf$gene_id)]
chemoreceptor_gtf <- gtf[which(gtf$gene_id%in%chemoreceptor),]
chemoreceptor_gtf_exon <- chemoreceptor_gtf[which(chemoreceptor_gtf$type=="exon"),]
chemoreceptor_exon <- as.data.frame(table(chemoreceptor_gtf_exon$gene_id))
chemoreceptor_exon <- chemoreceptor_exon[order(chemoreceptor_exon$Freq,decreasing=T),]
chemoreceptor_exon$type <- "receptor";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% OR_info$Symbol)]="OR";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% IR_info$Symbol)]="IR";
chemoreceptor_exon$type[which(chemoreceptor_exon$Var1%in% GR_info$Symbol)]="GR";

write.csv(chemoreceptor_exon,"chemoreceptor_exon_number.csv")

```
* annotate chemoreceptor manually by bulk RNAseq data and ATACseq data 

4. get the protein sequence of chemoreceptor
gffread GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_fixed.gtf -o- > merged.gff3
agat_convert_sp_gxf2gxf.pl -gff merged.gff3 -o agat_fixed.gff3
agat_sp_add_start_and_stop.pl --gff agat_fixed.gff3 --fasta GCF_003254395.2_Amel_HAv3.1_genomic.fna -o  agat_fixed_add_start_and_stop.gff
gffread agat_fixed_add_start_and_stop.gff -T -o agat_fixed_add_start_and_stop.gtf






grep -P '\btranscript_id\s+"[^"]+"'   GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED.gtf >  GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_fixed.gtf
gffread merged.gff3 -T -o my.gtf

LOC410623

>XP_006561077.2 ionotropic receptor 93a isoform X3 [Apis mellifera]
MISVLLLVWWINYGSSYNNFPSLITSNATMAVIIDKGFFSNKDEYQNATKVIQDLITDAVKKEMNLGSISIRVFRDMNVN
FKDYTILLSVATCYLTWRLHEVAQKEELTHFAITDPDCPRIPDTDGITVPSIVPGEELSQIFLDLRMTDILSWNVINILH
DDTFGDKATSSNDNVTILLSNANTCSRLVSDRDTISRVLKAISNKLPNKRMNLISRSIFSLRYGNTGSGRKSSVKKMLND
FHVEQLGHCFLVIATVDMVADVMSVANSLNMVHPGSQWLYVITNSVSGNLINTSFINLLAEGGNVAFMYNATNLDGFYKI
KLKCYIKDLIEALAKALEYSLKNEIELFKRMNEDEFEMIRLTKSKKRAELLKNVRIHLSRNTSASNSVCEQCLLWRFFSS
ITWGNFFSHDRNMAHLLDIGTWTPIIGVNLTDVIFPHIVHGFRGINLPIATYHNPPWQIISMSKTGKKLYEGLIFDAINY
LSMKLNFTYTVIMPETSQISRSWNTSQFAKLGEKIKEMTMSTTKKVPLEIIDLVRQKKVLLAACALTVNECGNTTFNYTV
PIFVQTYSFLTAKPSQLSRVLLFASPFTKETWACLAVSIIIMGPILYLIHKYSPYSTKASGLNSSWQCVWYVYGALLQQG
GMYLPQNDSARILIGMWWLVVMVLVATYSGSLVAFLTFPRMDTSILSVEDLIAHKDSISWGFPNGSFLEMYLQNAEEPKY
HVLFSRAERHNDTEEERLVERVKEGKHALIDWRSSLRFLMRKDFLLTGSCHFSLSMDEFLDEPIAMIIPYGSPYLSVINA
ELHRMLESGLMNKWITEKMPMKDKCWEAPGSNQMVNKRKVNVTDMQGIFFVLFIGITLAFFFLFCEFYCHRRKIAKERKL
IHPFVS

>XP_016772032.2 ionotropic receptor 93a isoform X1 [Apis mellifera]
MHYSQVSSVIYCANNTRGYFSLSIE

MISVLLLVWWINYGSSYNNFPSLITSNATMAVIIDKGFFSNKDEYQNATKVIQDL
ITDAVKKEMNLGSISIRVFRDMNVNFKDYTILLSVATCYLTWRLHEVAQKEELTHFAITDPDCPRIPDTDGITVPSIVPG
EELSQIFLDLRMTDILSWNVINILHDDTFGDKATSSNDNVTILLSNANTCSRLVSDRDTISRVLKAISNKLPNKRMNLIS
RSIFSLRYGNTGSGRKSSVKKMLNDFHVEQLGHCFLVIATVDMVADVMSVANSLNMVHPGSQWLYVITNSVSGNLINTSF
INLLAEGGNVAFMYNATNLDGFYKIKLKCYIKDLIEALAKALEYSLKNEIELFKRMNEDEFEMIRLTKSKKRAELLKNVR
IHLSRNTSASNSVCEQCLLWRFFSSITWGNFFSHDRNMAHLLDIGTWTPIIGVNLTDVIFPHIVHGFRGINLPIATYHNP
PWQIISMSKTGKKLYEGLIFDAINYLSMKLNFTYTVIMPETSQISRSWNTSQFAKLGEKIKEMTMSTTKKVPLEIIDLVR
QKKVLLAACALTVNECGNTTFNYTVPIFVQTYSFLTAKPSQLSRVLLFASPFTKETWACLAVSIIIMGPILYLIHKYSPY
STKASGLNSSWQCVWYVYGALLQQGGMYLPQNDSARILIGMWWLVVMVLVATYSGSLVAFLTFPRMDTSILSVEDLIAHK
DSISWGFPNGSFLEMYLQNAEEPKYHVLFSRAERHNDTEEERLVERVKEGKHALIDWRSSLRFLMRKDFLLTGSCHFSLS
MDEFLDEPIAMIIPYGSPYLSVINAELHRMLESGLMNKWITEKMPMKDKCWEAPGSNQMVNKRKVNVTDMQGIFFVLFIG
ITLAFFFLFCEFYCHRRKIAKERKLIHPFVS


seqtk subseq GCF_003254395.2_Amel_HAv3.1_genomic.fna test.bed > test.fa
Group1  11241   11244
Group1  9477    9480


#代码有所改变，参考：https://biopython.org/wiki/Alphabet
from Bio import Seq 
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
# read the input sequence
dna = open("test.fa").read().strip() #该文件内容为一条DNA编码序列
dna = Seq.Seq(dna)
#Seq对象为不可更改序列，mutableSeq对象为可变序列对象
# transcribe and translate
mrna = dna.transcribe() #转录
protein = mrna.translate() #翻译
# write the protein sequence to a file
protein_record = SeqRecord(protein, id='test',
                     description="Hemoglobin subunit alpha, Homo sapiens")
outfile = open("HBA_HUMAN.fasta", "w")
#SeqIO.write(protein_record, outfile,"fasta")
#SeqIO.write可将多个SeqRecord对象写入指定文件
outfile.close()


> usage: pfam_scan.py [-out OUT] [-outfmt {csv,json}] [-evalue EVALUE] [-cpu CPU] [-h] [-v]
                    fasta_file pfam_dir

[pfam_scan](https://github.com/aziele/pfam_scan)

pfam_scan.pl -fasta ChemoreceptorPeptides.fasta -dir /data/R02/nieyg/software/pfamdb

## get cell * gene/peak matrix by CellRanger 
1. create a cellranger-arc-compatible reference by cellranger
```
organism: "Apis mellifera"
genome: ["Amel_HAv3_1"]
input_fasta: ["/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna"]
input_gtf: ["/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED.gtf"]
input_motifs: "/md01/nieyg/ref/10X/Amel_HAv3.1/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
non_nuclear_contigs: ["MT"]
```
`nohup cellranger-arc mkref --config=mkref.contig &`

2.  Run cellranger-arc count
* Forager-NCBI-manually.pbs
```
#PBS -N Forager
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/Forager
cellranger-arc count --id=Forager-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/Forager/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```

* NE-NCBI-manually.pbs
```
#PBS -N NE
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/NE
cellranger-arc count --id=NE-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/NE/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```


* Nurse-NCBI-manually.pbs
```
#PBS -N Nurse
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/honeybee/data/cellranger/Nurse
cellranger-arc count --id=Nurse-NCBI-manually\
                       --reference=/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1 \
                       --libraries=/md01/nieyg/project/honeybee/data/cellranger/Nurse/libraries-add.csv \
                       --localcores=24 \
                       --localmem=30

```
