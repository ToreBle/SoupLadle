#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=STAR_SNP
#SBATCH --time=24:00:00
#SBATCH --output=output.%J.out
#SBATCH --array=1-5

OUT=$"/data/SNP_Multiplexing/data/bulkRNA/PBMC"
input=$OUT"/Samples/Samples.txt"

IDs=`awk -v id="${SLURM_ARRAY_TASK_ID}" 'NR==id{print $1,$2}' $input`
PATH_FASTQ=$(echo "$IDs" | awk '{print $1}')
SAMPLE=$(echo "$IDs" | awk '{print $2}')

echo "$PATH_FASTQ"
echo "$SAMPLE"
echo "$OUT"


#FASTQC
fastqc $PATH_FASTQ"/"$SAMPLE"_1.fq.gz"
fastqc $PATH_FASTQ"/"$SAMPLE"_2.fq.gz"

gunzip $PATH_FASTQ"/"$SAMPLE"_1.fq.gz"
gunzip $PATH_FASTQ"/"$SAMPLE"_2.fq.gz"


#STAR alignment
STAR --runThreadN 6 \
--genomeDir /Genome/index_GRCh38_98 \
--sjdbGTFfile /Genome/Homo_sapiens.GRCh38.98.chr.gtf \
--readFilesIn $PATH_FASTQ"/"$SAMPLE"_1.fq" $PATH_FASTQ"/"$SAMPLE"_2.fq" \
--outFileNamePrefix $OUT"/"$SAMPLE"/"$SAMPLE \
--outSAMattributes MD NH

#Sorting afgter STAR alignment
samtools view -b $OUT"/"$SAMPLE"/"$SAMPLE"Aligned.out.sam" > $OUT"/"$SAMPLE"/"$SAMPLE".bam"
samtools sort $OUT"/"$SAMPLE"/"$SAMPLE".bam" -o $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam"
samtools index $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bai"
samtools flagstat $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" > $OUT"/"$SAMPLE"/"$SAMPLE"_flagstat_alignment.txt"

#Deduplication and sorting
picard MarkDuplicates INPUT=$OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" OUTPUT=$OUT"/"$SAMPLE"/"$SAMPLE"_deduplicat.bam" METRICS_FILE=$OUT"/"$SAMPLE"/"$SAMPLE".txt" REMOVE_DUPLICATES=true
samtools sort $OUT"/"$SAMPLE"/"$SAMPLE"_deduplicat.bam" -o $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam"
samtools index $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bai"
samtools flagstat $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" > $OUT"/"$SAMPLE"/"$SAMPLE"_flagstat_deduplicate.txt"
rm $OUT"/"$SAMPLE"/"$SAMPLE"_deduplicat.bam"
rm $OUT"/"$SAMPLE"/"$SAMPLE".bam"

#Variant calling with freebayes
source activate souporcell
~/programmes/freebayes -f /Genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -iXu --min-mapping-quality 30 --min-base-quality 10 --read-max-mismatch-fraction 0.5 --read-snp-limit 5 --read-indel-limit 5 --min-coverage 10 --genotype-qualities $OUT"/"$SAMPLE"/"$SAMPLE"_sorted.bam" > $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf"
echo $SAMPLE > $OUT"/"$SAMPLE"/"$SAMPLE".txt"
bcftools reheader -s $OUT"/"$SAMPLE"/"$SAMPLE".txt" $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf" > $OUT"/"$SAMPLE"/"$SAMPLE".vcf"


#bgzip and indexing
bgzip -c $OUT"/"$SAMPLE"/"$SAMPLE".vcf" > $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf.gz" && tabix $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf.gz"
zcat $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf.gz"| awk -F"\t" '{if ($0 !~ /^#/) {print "chr"$0} else{print $0}}' | bgzip -c > $OUT"/"$SAMPLE"/"$SAMPLE".vcf.gz"
tabix $OUT"/"$SAMPLE"/"$SAMPLE".vcf.gz"
rm $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf.gz"
rm $OUT"/"$SAMPLE"/"$SAMPLE"_tmp.vcf.gz.tbi"
