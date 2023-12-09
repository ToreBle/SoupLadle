#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=SNP_p
#SBATCH --time=144:00:00
#SBATCH --output=output.%J.out
#SBATCH --array=1-1

export PATH=/home/uv525372/samtools-1.16.1/bin:$PATH
export PATH=/home/uv525372/programmes/minimap2/:$PATH
export PATH=/home/uv525372/programmes/bedtools2/bin:$PATH
export PATH=/home/uv525372/programmes/:$PATH
source activate souporcell

#input="/data/uv525372/multiplexing/Samples_SouporCell/Samples_PBMC_KH30.txt"
input="/data/uv525372/multiplexing/Samples_SouporCell/Samples_Heart_KH108.txt"


IDs=`awk -v id="${SLURM_ARRAY_TASK_ID}" 'NR==id{print $1,$2,$3,$4,$5}' $input`
SAMPLE=$(echo "$IDs" | awk '{print $1}')
CLUSTER=$(echo "$IDs" | awk '{print $2}')
BULKRNA=$(echo "$IDs" | awk '{print $3}')
WES=$(echo "$IDs" | awk '{print $4}')
REFERENCE=$(echo "$IDs" | awk '{print $5}')

#Working directory
mkdir "/data/SNP_Multiplexing/output/"$SAMPLE"_unfiltered/"
cd "/data/SNP_Multiplexing/output/"$SAMPLE"_unfiltered/"

#Consider the same cell barcode as defiend in the CellRanger pipeline (commented)
#awk -F ',' 'NR>1{print $7}' '/data/SNP_Multiplexing/data/multiplex/'$SAMPLE'/outs/multi/multiplexing_analysis/assignment_confidence_table.csv' > $SAMPLE_'all_cell_barcodes.tsv'

# Merge all pooled sample and the unassigned reads from cellranger multi
#samtools merge merged_cellranger_alignments.bam '/data/SNP_Multiplexing/data/multiplex/'$SAMPLE'/outs/per_sample_outs/P'*'/count/'*'bam' '/data/SNP_Multiplexing/data/multiplex/'$SAMPLE'/outs/multi/count/unassigned_alignments.bam'

# Extract only reads assigned to a cell by CellRanger
#samtools view merged_cellranger_alignments.bam | LC_ALL=C grep -F -f $SAMPLE_'all_cell_barcodes.tsv' | LC_ALL=C grep "xf:i:25" > body_filtered_sam

# Extract the BAM header and write to header_filted_sam
#samtools view -H merged_cellranger_alignments.bam > header_filtered_sam
#cat header_filtered_sam body_filtered_sam > filtered_alignments.sam
#samtools view -b filtered_alignments.sam > filtered_alignments.bam

#samtools sort filtered_alignments.bam -o $SAMPLE'_cellranger_alignments_sorted.bam'
#rm filtered_alignments.sam
#rm filtered_alignments.bam
#rm header_filtered_sam
#rm body_filtered_sam


samtools sort merged_cellranger_alignments.bam -o $SAMPLE'_cellranger_alignments_sorted.bam'
samtools index $SAMPLE'_cellranger_alignments_sorted.bam'

rm merged_cellranger_alignments.bam

#Comparisions of different running modes
mkdir output_souporcell_minimap
python /home/uv525372/programmes/souporcell/souporcell_pipeline.py -i $SAMPLE'_cellranger_alignments_sorted.bam' -b $SAMPLE_'all_cell_barcodes.tsv' -f $REFERENCE -t 8 -o ./output_souporcell_minimap/ -k $CLUSTER #--skip_remap SKIP_REMAP --ignore True

mkdir output_souporcell_SNP_only
python /home/uv525372/programmes/souporcell/souporcell_pipeline.py -i $SAMPLE'_cellranger_alignments_sorted.bam' -b $SAMPLE_'all_cell_barcodes.tsv' -f $REFERENCE -t 8 -o ./output_souporcell_SNP_only/ -k $CLUSTER --skip_remap SKIP_REMAP --ignore True


## Vireo
mkdir ./vireo_WES
bcftools merge -Oz /data/SNP_Multiplexing/data/WES/Heart/*/M3????.vcf.gz > ./vireo_WES/WES_merged.vcf.gz
bcftools index ./vireo_WES/WES_merged.vcf.gz
bcftools view ./vireo_WES/WES_merged.vcf.gz -R ./cellSNP/cellSNP.cells.vcf.gz -Oz -o './vireo_WES/'$SAMPLE'_sc_overlap.vcf.gz'
vireo -c ./cellSNP/cellSNP.cells.vcf.gz -d './vireo_WES/'$SAMPLE'_sc_overlap.vcf.gz' -t GT --forceLearnGT -o ./vireo_WES
GTbarcode -i ./vireo_WES/GT_donors.vireo.vcf.gz -o ./vireo_WES/GT_barcodes.tsv --randSeed 1

#Vireo with sample assignment of bulkRNA-seq data
mkdir ./vireo_bulkRNA
bcftools merge -Oz /data/SNP_Multiplexing/data/bulkRNA/Kidney/*/KH*.vcf.gz > './vireo_bulkRNA/'$SAMPLE'_merged.vcf.gz'
bcftools index './vireo_bulkRNA/'$SAMPLE'_merged.vcf.gz'
#Filter for regions that overlap with the scRNA-seq SNPs
bcftools view './vireo_bulkRNA/'$SAMPLE'_merged.vcf.gz' -R ./cellSNP/cellSNP.cells.vcf.gz -Oz -o './vireo_bulkRNA/'$SAMPLE'_sc_overlap.vcf.gz'
vireo -c ./cellSNP/cellSNP.cells.vcf.gz -d './vireo_bulkRNA/'$SAMPLE'_sc_overlap.vcf.gz' -t GT --forceLearnGT -o ./vireo_bulkRNA
GTbarcode -i ./vireo_bulkRNA/GT_donors.vireo.vcf.gz -o ./vireo_bulkRNA/GT_barcodes.tsv --randSeed 1

#vireo with freebayes and vartrix SNP
vireo --vartrixData=./output_souporcell_SNP_only/alt.mtx,./output_souporcell_SNP_only/ref.mtx,$SAMPLE_'all_cell_barcodes.tsv',./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -N $CLUSTER -o ./vireo_vartrix
GTbarcode -i ./vireo_vartrix/GT_donors.vireo.vcf.gz -o ./vireo_vartrix/GT_barcodes.tsv --randSeed 1

#vireo with SNPs only and cellSNP
mkdir ./vireo_SNP
vireo -c ./cellSNP/cellSNP.cells.vcf.gz -N $CLUSTER -o ./vireo_SNP_only
GTbarcode -i ./vireo_SNP_only/GT_donors.vireo.vcf.gz -o ./vireo_SNP_only/GT_barcodes.tsv --randSeed 1


#vireo with vartrix and bulkRNA
mkdir ./vartrix_vireo_bulkRNA
bcftools view './vireo_bulkRNA/'$SAMPLE'_merged.vcf.gz' -R ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -Oz -o './vartrix_vireo_bulkRNA/'$SAMPLE'_sc_overlap.vcf.gz'
vireo --vartrixData=./output_souporcell_SNP_only/alt.mtx,./output_souporcell_SNP_only/ref.mtx,$SAMPLE_'all_cell_barcodes.tsv',./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -d './vartrix_vireo_bulkRNA/'$SAMPLE'_sc_overlap.vcf.gz' -t GT --forceLearnGT -o ./vartrix_vireo_bulkRNA
GTbarcode -i ./vartrix_vireo_bulkRNA/GT_donors.vireo.vcf.gz -o ./vartrix_vireo_bulkRNA/GT_barcodes.tsv --randSeed 1

