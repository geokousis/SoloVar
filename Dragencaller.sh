#!/bin/bash
bwa_files="/home/kousis/work/meth/Dra-GATK/mutec/BWA_files"
bwa_files="/home/kousis/work/meth/ILCEM/mapping_normal/bwa_files"
N_T=25
reference="/home/kousis/work/meth/Dra-GATK/mutec/reference/Homo_sapiens_assembly38.fasta"
intervals="/home/kousis/work/meth/ILCEM/S07604514_hg38_1/S07604514_Regions.clean.bed"
# Create chromosome list
output_file="chromosome_list.list"
mkdir vcf
grep "^>" "$reference" | awk '{print substr($1, 2)}' | sort | uniq | grep '^chr' > "$output_file"

pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/*recali*.bam > tmp.txt
popd > /dev/null


Data_Pre_Processing() {
    local reference="$1"
    gatk ComposeSTRTableFile \
         -R "$reference" \
         -O vcf/str_table.tsv
    if [ $? -ne 0 ]; then
        echo "Error during Data pre processing"
        exit 1
    fi
}

Haplotype_caller() {
    local reference="$1"
    local input_bam="$2"
    local base_name=$(basename "$input_bam" | sed 's/^md_//')
    gatk CalibrateDragstrModel \
         -R "$reference" \
         -I "$input_bam" \
         -str vcf/str_table.tsv \
         -O vcf/"${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during Calibration"
        exit 1
    fi

    gatk HaplotypeCaller \
         -R "$reference" \
         -I "$input_bam" \
         -L "$intervals" \
         --interval-padding 50 \
         -O vcf/"${base_name}.g.vcf" \
         -ERC GVCF \
         --tmp-dir tmp_v \
         --create-output-variant-index \
         --dragen-mode true \
         --dragstr-params-path vcf/"${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Uncomment if per-sample filtration is needed
    # gatk VariantFiltration \
    #      -V vcf/"${base_name}.g.vcf" \
    #      --filter-expression "QUAL < 50" \
    #      --filter-name "DRAGENHardQUAL" \
    #      --filter-expression "DP < 60" \
    #      --filter-name "LowDP" \
    #      -O vcf/filtered_"${base_name}.g.vcf" \
    #      --create-output-variant-index \
    #      --tmp-dir tmp_v
}

Merge() {
    local reference="$1"
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
         -L chromosome_list.list \
         --genomicsdb-workspace-path vcf/my_database \
         --tmp-dir tmp_v \
         --sample-name-map vcf/cohort.sample_map \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during GDB"
        exit 1
    fi
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
         -R "$reference" \
         -V gendb://vcf/my_database \
         -O vcf/Merged.vcf \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Select SNPs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V vcf/Merged.vcf \
         --select-type-to-include SNP \
         -O SNPs.vcf
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection SNPs"
        exit 1
    fi

    # Select INDELs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V vcf/Merged.vcf \
         --select-type-to-include INDEL \
         -O INDELs.vcf
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection Indels"
        exit 1
    fi
}

# Create cohort.sample_map file
create_sample_map() {
    while read -r bam_file; do
        base_name=$(basename "$bam_file" | sed 's/^md_//')
        echo -e "${base_name}\tvcf/${base_name}.g.vcf"
    done < "${bwa_files}/tmp.txt" > vcf/cohort.sample_map
}

# Data pre-processing
Data_Pre_Processing "$reference"

# Process each BAM file
count=0
while read -r bam_file; do
    running_jobs=$(jobs -p | wc -l)
    while [ "$running_jobs" -ge "$N_T" ]; do
        sleep 1
        running_jobs=$(jobs -p | wc -l)
    done
    Haplotype_caller "$reference" "$bam_file" &
    count=$((count + 1))
done < "${bwa_files}/tmp.txt"

wait

# Create sample map after all HaplotypeCaller jobs are done
create_sample_map
# Merge the GVCF files
Merge "$reference"

echo "All BAM files have been processed and merged."
