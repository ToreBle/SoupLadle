#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=WES
#SBATCH --time=24:00:00
#SBATCH --output=output.%J.out

source activate cellSNP
export PATH=/home/uv525372/programmes/htslib:$PATH
export PATH=/home/uv525372/programmes/bcftools/:$PATH
module load R

cd /data/SNP_Multiplexing/output/KH30_unfiltered/
Rscript --vanilla /data/SNP_Multiplexing/output/KH30_unfiltered/SNP_assignment/vcf_subsampling.R
bgzip -c KH30_merged_bulkRNA.vcf > KH30_merged_bulkRNA.vcf.gz && tabix KH30_merged_bulkRNA.vcf.gz


cd /data/SNP_Multiplexing/output/KH30_unfiltered/SNP_assignment/
seed=$(ls -1 | grep -E '^SNP_fraction_PBMC_bulk_RNA_.*\.tab$' | sed -E 's/^SNP_fraction_PBMC_bulk_RNA_(.*)\.tab$/\1/')
#cd /data/SNP_Multiplexing/output/KH30_unfiltered/
#fractions=(1_1 0.8 0.6 0.4 0.2 0.1)

fractions=($seed)

for fraction in "${fractions[@]}"; do
  #fraction_parts=(${fraction//_/ })
  filename="SNP_fraction_PBMC_bulk_RNA_${fraction}"
  base_dir="/data/SNP_Multiplexing/output/KH30_unfiltered"
  cd /data/SNP_Multiplexing/output/KH30_unfiltered/SNP_assignment/
  bcftools view "${base_dir}/KH30_merged_bulkRNA.vcf.gz" -R "${filename}.tab" -Oz -o "${filename}.vcf.gz"
  bcftools index "${filename}.vcf.gz"
  bcftools view "${filename}.vcf.gz" -R "${base_dir}/output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz" -Oz -o "${filename}_sc_overlap.vcf.gz"
  mkdir "${fraction}"
  vireo --vartrixData="${base_dir}/output_souporcell_SNP_only/alt.mtx","${base_dir}/output_souporcell_SNP_only/ref.mtx","${base_dir}/all_cell_barcodes.tsv","${base_dir}/output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz" -d "${filename}_sc_overlap.vcf.gz" -t GT --forceLearnGT -o "./${fraction}/vartrix_bulkRNA"
  mkdir "${fraction}"
  vireo -c "${base_dir}/cellSNP/cellSNP.cells.vcf.gz" -d "${filename}_sc_overlap.vcf.gz" -t GT --forceLearnGT -o "./${fraction}/vireo_bulkRNA"
done

#Patient_Assignment_V2
