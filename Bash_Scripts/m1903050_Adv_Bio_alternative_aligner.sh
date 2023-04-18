##!/usr/bin/env bash
#First create a project environment with parent and subdirectories
mkdir -p dnaseq/data/aligned/BAM_files  dnaseq/data/bed dnaseq/data/reference/bowtie_index dnaseq/data/trimmed dnaseq/data/untrimmed/compressed dnaseq/results/annotated_data/annovar dnaseq/results/annotated_data/snpEff dnaseq/results/statistics dnaseq/results/VCF dnaseq/results/fastqc/untrimmed dnaseq/results/fastqc/trimmed

#Download the required files
wget -P dnaseq/data/untrimmed/compressed https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget -P dnaseq/data/untrimmed/compressed https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget -P dnaseq/data/bed https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget -P dnaseq/data/reference http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

#Extract fastq.qz files to .fastq
zcat -d dnaseq/data/untrimmed/compressed/NGS0001.R1.fastq.qz > dnaseq/data/untrimmed/NGS0001.R1.fastq
zcat -d dnaseq/data/untrimmed/compressed/NGS0001.R2.fastq.qz > dnaseq/data/untrimmed/NGS0001.R2.fastq

#Perform quality assessment of raw data
for raw_fastq in dnaseq/data/untrimmed/NGS0001.R*;do
fastqc -o dnaseq/results/fastqc/untrimmed -t 4 ${raw_fastq}
done

#Use trimmomatic to trim raw data
trimmomatic PE \
 -threads 4 \
 -phred33 \
dnaseq/data/untrimmed/NGS0001.R1.fastq dnaseq/data/untrimmed/NGS0001.R2.fastq \
-baseout dnaseq/data/trimmed/NGS0001R \
TRAILING:25 MINLEN:50

#Perform quality assessment of raw data
for trimmed_fastq in dnaseq/data/trimmed/*P*;do
fastqc -o dnaseq/results/fastqc/trimmed -t 4 ${trimmed_fastq}
done


#using Bowtie2 index the reference genome
zcat dnaseq/data/reference/hg19.fa.gz > dnaseq/data/reference/hg19.fa
bowtie2-build --threads 4 dnaseq/data/reference/hg19.fa dnaseq/data/reference/bowtie_index/hg19

#Bowtie2 alignment
bowtie2 -p 4 --rg "ID:@RG11V6WR1:111:D1375ACXX:1.NGS0001" --rg "SM:NGS001,PL:ILLUMINA,PU:11V6WR1" -x dnaseq/data/reference/bowtie_index/hg19 -1 dnaseq/data/trimmed/NGS0001R_1P -2 dnaseq/data/trimmed/NGS0001R_2P -S dnaseq/data/aligned/NGS0001.sam

#Converting sam to bam
samtools view -h -b dnaseq/data/aligned/NGS0001.sam > dnaseq/data/aligned/BAM_files/NGS0001.bam
samtools sort dnaseq/data/aligned/BAM_files/NGS0001.bam > dnaseq/data/aligned/BAM_files/NGS0001_sorted.bam
samtools index dnaseq/data/aligned/BAM_files/NGS0001_sorted.bam

#Marking duplicates
picard MarkDuplicates I=dnaseq/data/aligned/BAM_files/NGS0001_sorted.bam O=dnaseq/data/aligned/BAM_files/NGS0001_sorted_marked.bam M=dnaseq/data/aligned/BAM_files/marked_dup_metrics.txt
samtools index dnaseq/data/aligned/BAM_files/NGS0001_sorted_marked.bam

#Filter BAM based on mapping quality and bitwise flags
samtools view -F 1796 -q 20 -o dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam dnaseq/data/aligned/BAM_files/NGS0001_sorted_marked.bam
samtools index dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam


#Apply statistical analysis on the marked and filtered bam files (flagstat, idxstats, depth of coverage, insert size)
#Flagstat
#Before filter
samtools flagstat dnaseq/data/aligned/BAM_files/NGS0001_sorted_marked.bam > dnaseq/results/statistics/NGS0001_sorted_marked_flagstat.txt

#After filter
samtools flagstat dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam > dnaseq/results/statistics/NGS0001_sorted_filtered_flagstat.txt

#idxstats
samtools idxstats dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam > dnaseq/results/statistics/NGS0001_sorted_filtered_idxstats.txt

#CollectInsertSizeMetrics
picard CollectInsertSizeMetrics I=dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam O=dnaseq/results/statistics/NGS0001_sorted_filtered_CollectInsertSizeMetrics.txt H=dnaseq/results/statistics/CollectInsertSize_Metrics.pdf M=0.5

#depth
samtools depth dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam > dnaseq/results/statistics/NGS0001_sorted_filtered_depth.txt


#Variant calling with Freebayes
samtools faidx dnaseq/data/reference/hg19.fa
freebayes --bam dnaseq/data/aligned/BAM_files/NGS0001_sorted_filtered.bam --fasta-reference dnaseq/data/reference/hg19.fa --vcf dnaseq/results/VCF/NGS0001.vcf
bgzip -c dnaseq/results/VCF/NGS0001.vcf > dnaseq/results/VCF/NGS0001.vcf.gz
tabix -p vcf dnaseq/results/VCF/NGS0001.vcf.gz

#Filtering the VCF
vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
dnaseq/results/VCF/NGS0001.vcf.gz > dnaseq/results/VCF/NGS0001_filtered.vcf

bedtools intersect -header -wa -a dnaseq/results/VCF/NGS0001_filtered.vcf \
-b dnaseq/data/bed/annotation.bed \
dnaseq/results/VCF/NGS0001_filtered.vcf > dnaseq/results/VCF/bedtools_intersect.txt
bgzip -c dnaseq/results/VCF/NGS0001_filtered.vcf > dnaseq/results/VCF/NGS0001_filtered.vcf.gz
tabix -p vcf dnaseq/results/VCF/NGS0001_filtered.vcf.gz


#Annotation with annovar
Download annovar.latest.tar.gz
wget -P dnaseq/results/annotated_data/ http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

#Extract annovar.latest.tar.gz
tar -zxvf dnaseq/results/annotated_data/annovar.latest.tar.gz -C dnaseq/results/annotated_data/

#Download annovar databases
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro dnaseq/results/annotated_data/annovar/humandb/
perl dnaseq/results/annotated_data/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp137 dnaseq/results/annotated_data/annovar/humandb/

#VCF to annovar input format
perl dnaseq/results/annotated_data/annovar/convert2annovar.pl -format vcf4 dnaseq/results/VCF/NGS0001_filtered.vcf.gz > dnaseq/results/VCF/NGS0001_filtered.avinput

#CSV output
perl dnaseq/results/annotated_data/annovar/table_annovar.pl dnaseq/results/VCF/NGS0001_filtered.avinput dnaseq/results/annotated_data/annovar/humandb/ -buildver hg19 \
 -out dnaseq/results/VCF/NGS0001_filtered -remove \
 -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,snp137 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout


#Annotate with snpEff
#Download snpEff
wget -P dnaseq/results/annotated_data/ https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
wget -P dnaseq/results/annotated_data/snpEff/ https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz

#Extract snpEff_latest_core.zip
unzip dnaseq/results/annotated_data/snpEff_latest_core.zip -d dnaseq/results/annotated_data/

#Annotate using snpEff
java -jar dnaseq/results/annotated_data/snpEff/snpEff.jar -c dnaseq/results/annotated_data/snpEff/snpEff.config -v GRCh37.75 dnaseq/results/VCF/NGS0001_filtered.vcf > dnaseq/results/VCF/NGS0001_snpeff.ann.vcf

#Annotate the snpeff output file using snpSift
gzip -d dnaseq/results/annotated_data/snpEff/00-common_all.vcf.gz
java -Xmx8g -jar dnaseq/results/annotated_data/snpEff/SnpSift.jar annotate -v dnaseq/results/annotated_data/snpEff/00-common_all.vcf dnaseq/results/VCF/NGS0001_snpeff.ann.vcf > dnaseq/results/VCF/NGS0001_snpsift_dbsnp.ann.vcf

#Filtering exonic variants not seen in dbSNP
java -jar dnaseq/results/annotated_data/snpEff/SnpSift.jar filter "(ANN[*].BIOTYPE = 'protein_coding')" dnaseq/results/VCF/NGS0001_snpsift_dbsnp.ann.vcf > dnaseq/results/VCF/NGS0001_snpsift_dbsnp_exonic_filtering.ann.vcf
java -jar dnaseq/results/annotated_data/snpEff/SnpSift.jar filter "(ID !~ 'rs')" dnaseq/results/VCF/NGS0001_snpsift_dbsnp_exonic_filtering.ann.vcf  > dnaseq/results/VCF/NGS0001_snpsift_dbsnp_exonic_filtering_no_dbSNP.ann.vcf
