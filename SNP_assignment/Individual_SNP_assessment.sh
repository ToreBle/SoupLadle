#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=WES
#SBATCH --time=24:00:00
#SBATCH --output=output.%J.out

export PATH=/home/uv525372/programmes/htslib:$PATH
export PATH=/home/uv525372/programmes/bcftools/:$PATH
module load R


#CellSNP vs Vartrix
cd /data/SNP_Multiplexing/output/KH30_unfiltered/
bcftools index ./cellSNP/cellSNP.cells.vcf.gz
bcftools isec -n-2 -O z -p "./SNP_evaluation/cellSNP_Vartrix_isec/" ./cellSNP/cellSNP.cells.vcf.gz ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz
Rscript --vanilla /data/uv525372/multiplexing/SNP_calling_evaluation.R "./SNP_evaluation/cellSNP_Vartrix_isec/" "PBMC"

cd /data/SNP_Multiplexing/output/KH108_unfiltered/
bcftools index ./cellSNP/cellSNP.cells.vcf.gz
bcftools isec -n-2 -O z -p "./SNP_evaluation/cellSNP_Vartrix_isec/" ./cellSNP/cellSNP.cells.vcf.gz ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz
Rscript --vanilla /data/uv525372/multiplexing/SNP_calling_evaluation.R "./SNP_evaluation/cellSNP_Vartrix_isec/" "Heart"


#WES
process_WES_samples() {
    local data_type=$1
    local samples=("${@:2}")
    local base_path="/data/SNP_Multiplexing/data/WES/$data_type"

    for sample in "${samples[@]}"; do
        #bcftools view ./vireo_WES/"$sample"_GATK.vcf.gz -R ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -Oz -o "$base_path/$sample/${sample}_vartrix_overlap.vcf.gz"
        bgzip -c "$base_path/$sample"/*GATK.vcf > "$base_path/$sample/$sample.vcf.gz"
        bcftools index "$base_path/$sample/$sample.vcf.gz"
        bcftools view "$base_path/$sample/$sample.vcf.gz" -R ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -Oz -o "$base_path/$sample/${sample}_vartrix_overlap.vcf.gz"
        bcftools index "$base_path/$sample/${sample}_vartrix_overlap.vcf.gz"
        #bcftools view "$base_path/$sample/$sample.vcf.gz" -R ./cellSNP/cellSNP.cells.vcf.gz -Oz -o "$base_path/$sample/${sample}_cellSNP_overlap.vcf.gz"
        #bcftools index "$base_path/$sample/${sample}_cellSNP_overlap.vcf.gz"
    done

    bcftools isec -n-4 -O z -p "./SNP_evaluation/WES_vartrix_isec/" "$base_path"/*/*_vartrix_overlap.vcf.gz
    bcftools isec -n-4 -O z -p "./SNP_evaluation/WES_isec/" "$base_path"/*/M3????.vcf.gz
    #bcftools isec -n-4 -O z -p "./SNP_evaluation/WES_cellSNP_isec/" "$base_path"/*/*_cellSNP_overlap.vcf.gz
}


cd /data/SNP_Multiplexing/output/KH30_unfiltered/
PBMC_samples=("M34287" "M34288" "M34289" "M34290" "M34291")
process_WES_samples "PBMC" "${PBMC_samples[@]}"

cd /data/SNP_Multiplexing/output/KH108_unfiltered/
Heart_samples=("M36483" "M36484" "M36485" "M36486" "M36487")
process_WES_samples "Heart" "${Heart_samples[@]}"


#bulkRNA
process_bulk_samples() {
    local data_type=$1
    local samples=("${@:2}")
    local base_path="/data/SNP_Multiplexing/data/bulkRNA/$data_type"

    for sample in "${samples[@]}"; do
        bcftools view "$base_path/$sample/$sample.vcf.gz" -R ./output_souporcell_SNP_only/souporcell_merged_sorted_vcf.vcf.gz -Oz -o "$base_path/$sample/${sample}_vartrix_overlap.vcf.gz"
        bcftools index "$base_path/$sample/${sample}_vartrix_overlap.vcf.gz"
        #bcftools view "$base_path/$sample/$sample.vcf.gz" -R ./cellSNP/cellSNP.cells.vcf.gz -Oz -o "$base_path/$sample/${sample}_cellSNP_overlap.vcf.gz"
        #bcftools index "$base_path/$sample/${sample}_cellSNP_overlap.vcf.gz"
    done

    bcftools isec -n-4 -O z -p "./SNP_evaluation/bulk_RNA_vartrix_isec/" "$base_path"/*/*_vartrix_overlap.vcf.gz
    bcftools isec -n-4 -O z -p "./SNP_evaluation/bulk_RNA_isec/" "$base_path"/*/KH???.vcf.gz
    #bcftools isec -n-4 -O z -p "./SNP_evaluation/bulk_RNA_cellSNP_isec/" "$base_path"/*/*_cellSNP_overlap.vcf.gz
}


cd /data/SNP_Multiplexing/output/KH30_unfiltered/
PBMC_samples=("KH_75" "KH_76" "KH_77" "KH_78" "KH_79")
process_bulk_samples "PBMC" "${PBMC_samples[@]}"
Rscript --vanilla /data/uv525372/multiplexing/Transcript_Coverage_WES_bulkRNA.R "/data/Tore/scSLAM_seq/Genome/Homo_sapiens.GRCh38.98.chr.gtf" "/data/SNP_Multiplexing/output/KH30_unfiltered/" "PBMC"
Rscript --vanilla /data/uv525372/multiplexing/SNP_bulk_demultiplexing_pre_evaluation.R "/data/Tore/scSLAM_seq/Genome/Homo_sapiens.GRCh38.98.chr.gtf" "/data/SNP_Multiplexing/output/KH30_unfiltered/" "PBMC"


cd /data/SNP_Multiplexing/output/KH108_unfiltered/
Heart_samples=("KH113" "KH114" "KH115" "KH116" "KH117")
process_bulk_samples "Heart" "${Heart_samples[@]}"
Rscript --vanilla /data/uv525372/multiplexing/Transcript_Coverage_WES_bulkRNA.R "/data/Tore/scSLAM_seq/Genome/Homo_sapiens.GRCh38.98.chr.gtf" "/data/SNP_Multiplexing/output/KH108_unfiltered/" "Heart"
Rscript --vanilla /data/uv525372/multiplexing/SNP_bulk_demultiplexing_pre_evaluation.R "/data/Tore/scSLAM_seq/Genome/Homo_sapiens.GRCh38.98.chr.gtf" "/data/SNP_Multiplexing/output/KH108_unfiltered/" "Heart"
