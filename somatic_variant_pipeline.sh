#!/bin/bash

# This integrated pipeline performs somatic variant calling using Mutect2
# followed by variant filtering. It reads configuration from a YAML file.
#
# Expected YAML keys include:
#   bwa_files:         Directory containing analysis-ready BAM files (e.g., f_*recalibrated*.bam)
#   output_directory:  Base directory for outputs (VCFs, logs, filtered VCFs, etc.)
#   threads:           Number of threads to use
#   reference:         Reference genome FASTA file
#   intervals:         Intervals file (BED or interval list)
#   pon:               Panel-of-Normals VCF file
#   gr:                Germline resource VCF file (e.g., gnomAD AF-only)
#   subgr:             Common sites VCF (for GetPileupSummaries)
#   gatk_path:         Path to GATK executable (or "gatk" if in PATH)
#   learn_orientation: yes/no (if yes, run LearnReadOrientationModel to model FFPE artifacts)
#
# The pipeline consists of two parts:
#
# Part 1: Mutect2 Variant Calling
#   - For each analysis-ready BAM file (pattern: f_*recalibrated*.bam in bwa_files),
#     run GATK Mutect2 to produce a VCF and an F1R2 tar file (output to output_directory/vcf).
#
# Part 2: Variant Filtering
#   - Optionally learn the read orientation model (if learn_orientation is "yes").
#   - Run GetPileupSummaries on each tumor BAM (using subgr and intervals).
#   - Run CalculateContamination on each resulting table.
#   - Run FilterMutectCalls to filter the raw calls (using the contamination table and,
#     if available, the orientation model as ob-priors).
#   - Select passed variants from each filtered VCF.
#
# Usage:
#   ./integrated_somatic_variant_pipeline.sh config.yaml

########################################
# YAML Parser (simple key: value pairs)
########################################
parse_yaml() {
    local yaml_file=$1
    while IFS=":" read -r key value; do
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        if [[ -n "$key" && ! "$key" =~ ^# ]]; then
            eval "${key}=\"${value}\""
        fi
    done < "$yaml_file"
}

########################################
# Logging helper
########################################
# Global settings
DRY_RUN=${DRY_RUN:-false}    # Set to 'true' to enable dry-run mode (show commands without executing)
########################################
# Execute helper function (with logging & dry-run)
########################################
execute() {
    local cmd="$1"
    local description="$2"
    local runtime_flag="${3:-}"
    local verbose_flag="${4:-}"
    local tool
    tool=$(echo "$cmd" | awk '{print $1}')

    # Dry-run: show but do not execute
    if [[ "$DRY_RUN" == "true" ]]; then
        echo -e "\033[2m[DRY RUN]\033[0m \033[1m$description\033[0m"
        echo "[CMD] $cmd"
        return 0
    fi

    # Verify tool availability before running
    if ! command -v "$tool" &>/dev/null; then
        echo -e "\033[1mError: Required tool '$tool' not found.\033[0m"
        exit 127
    fi

    # Log start: [INFO] dim, description bold
    echo -e "\033[2m[INFO]\033[0m $(date '+%Y-%m-%d %H:%M:%S') - \033[1m$description\033[0m" | tee -a "$LOG_FILE"
    echo "[CMD] $cmd" >> "$LOG_FILE"

    local start_time end_time runtime
    if [[ "$runtime_flag" == "runtime" ]]; then
        start_time=$(date +%s)
    fi

    if [[ "$verbose_flag" == "verbose" ]]; then
        bash -c "$cmd" 2>&1 | tee -a "$LOG_FILE"
        status=${PIPESTATUS[0]}
    else
        bash -c "$cmd" >> "$LOG_FILE" 2>&1
        status=$?
    fi

    if [[ "$runtime_flag" == "runtime" ]]; then
        end_time=$(date +%s)
        runtime=$((end_time - start_time))
        echo -e "\033[2m[INFO]\033[0m Command completed in $runtime seconds" | tee -a "$LOG_FILE"
    fi

    if [[ $status -ne 0 ]]; then
        echo -e "\033[1mError: '$description' failed (exit code $status). See log: $LOG_FILE\033[0m"
        exit $status
    fi
}

########################################
# Cleanup on exit
########################################
cleanup() {
    echo -e "\033[2m[INFO]\033[0m Exiting at $(date '+%Y-%m-%d %H:%M:%S')"
}
trap cleanup EXIT


########################################
# Part 1: Mutect2 Variant Calling
########################################

# Process a single tumor BAM file with Mutect2.
Mutect2_Calling_single() {
    local tumor_bam="$1"
    local base_name
    base_name=$(basename "$tumor_bam" | sed 's/^f_//; s/\.bam$//')
    echo "Running Mutect2 for $base_name" | tee -a "$LOG_FILE"
    execute "$gatk_path Mutect2 -R \"$reference\" -I \"$tumor_bam\" -L \"$intervals\" --interval-padding \"$padding\" --germline-resource \"$gr\" --panel-of-normals \"$pon\" -O \"${output_directory}/vcf/${base_name}.g.vcf\" --create-output-variant-index --f1r2-tar-gz \"${output_directory}/vcf/f1r2_${base_name}.tar.gz\" --tmp-dir \"${output_directory}/tmp_v\"" "Mutect2 variant calling for $base_name"
}

# Parallel version: uses GNU parallel to process all tumor BAM files.
Mutect2_Calling() {
    echo "Starting Mutect2 variant calling (parallel mode)..." | tee -a $LOG_FILE
    
    # Generate a list of analysis-ready BAM files (pattern: f_*recalibrated*.bam)
    pushd "$bwa_files" > /dev/null
    ls -1 -d "$(pwd)"/f_*recalibrated*.bam > bam_list.txt
    popd > /dev/null

    # Export functions and required variables for GNU parallel.
    export -f Mutect2_Calling_single execute
    export output_directory threads reference intervals pon gr gatk_path padding LOG_FILE

    # Run Mutect2 variant calling in parallel using the list of BAMs.
    parallel --jobs "$threads" Mutect2_Calling_single ::: $(cat "${bwa_files}/bam_list.txt")
    
    echo "Mutect2 variant calling completed (parallel mode)." | tee -a $LOG_FILE
}

# Serial version: loops over each tumor BAM file one by one.
Mutect2_Calling_Serial() {
    echo "Starting Mutect2 variant calling (serial mode)..." | tee -a $LOG_FILE
    
    pushd "$bwa_files" > /dev/null
    ls -1 -d "$(pwd)"/f_*recalibrated*.bam > bam_list.txt
    popd > /dev/null
    
    while read -r tumor_bam; do
        Mutect2_Calling_single "$tumor_bam"
    done < "${bwa_files}/bam_list.txt"
    
    echo "Mutect2 variant calling completed (serial mode)." | tee -a $LOG_FILE
}

# Wrapper: Checks the 'parallel_mutect' variable and calls the appropriate function.
Run_Mutect2_Calling() {
    if [ "$parallel_mutect" == "yes" ]; then
        Mutect2_Calling
    else
        Mutect2_Calling_Serial
    fi
}

########################################
# Part 2: Variant Filtering Functions
########################################
# Optional: LearnReadOrientationModel (to model FFPE artifacts)
LearnReadOrientationModel_Filt() {
    local input_tar="$1"
    local base_name
    base_name=$(basename "$input_tar" | sed 's/^f1r2_//; s/\.tar\.gz$//')
    local output_model="${output_directory}/vcf/${base_name}.tar.gz"
    echo "Running LearnReadOrientationModel for $input_tar -> $(basename "$output_model")" | tee -a $LOG_FILE
    execute "$gatk_path LearnReadOrientationModel -I \"$input_tar\" -O \"$output_model\"" "LearnReadOrientationModel for $input_tar"
}

GetPileupSummaries_Filt() {
    local tumor="$1"
    local base_name
    base_name=$(basename "$tumor" | sed 's/^f_//; s/\.bam$//')
    local output_table="${output_directory}/filter_vcf/${base_name}.table"

    echo "Running GetPileupSummaries for $tumor -> $(basename "$output_table")" | tee -a $LOG_FILE
    execute "$gatk_path GetPileupSummaries -I \"$tumor\" -V \"$subgr\" -L \"$intervals\" -O \"$output_table\"" "GetPileupSummaries for $tumor"
}

CalculateContamination_Filt() {
    local table="$1"
    local base_name
    base_name=$(basename "$table" .table)
    local output_contam="${output_directory}/vcf/${base_name}.txt"

    echo "Running CalculateContamination for $table -> $(basename "$output_contam")" | tee -a $LOG_FILE
    execute "$gatk_path CalculateContamination -I \"$table\" -tumor-segmentation \"${output_contam}.seg\" -O \"$output_contam\"" "CalculateContamination for $table"
}

FilterMutectCalls_Filt() {
    local vcf="$1"
    local cont_table="$2"
    local ob="$3"   # may be blank if learn_orientation not enabled
    local base_name
    base_name=$(basename "$vcf")
    local output_filtered="${output_directory}/filter_vcf/filtered_${base_name}"

    echo "Filtering Mutect calls for $vcf -> $(basename "$output_filtered")" | tee -a $LOG_FILE
    if [ -z "$ob" ]; then
        execute "$gatk_path FilterMutectCalls -R \"$reference\" -V \"$vcf\" --contamination-table \"$cont_table\" --tumor-segmentation \"${cont_table}.seg\" --stats \"${vcf}.stats\" -O \"$output_filtered\"" "FilterMutectCalls for $vcf (no orientation model)"
    else
        execute "$gatk_path FilterMutectCalls -R \"$reference\" -V \"$vcf\" --contamination-table \"$cont_table\" --tumor-segmentation \"${cont_table}.seg\" --stats \"${vcf}.stats\" --ob-priors \"$ob\" -O \"$output_filtered\"" "FilterMutectCalls for $vcf (with orientation model)"
    fi
}


# Select passed variants from the filtered VCF
SelectPassedVariants_Filt() {
    local input_vcf="$1"
    local base_name
    base_name=$(basename "$input_vcf" .vcf)
    local output_vcf="${output_directory}/filtered_pass_vcf/filtered_pass_${base_name}.vcf"

    echo "Selecting passed variants from $input_vcf -> $output_vcf" | tee -a $LOG_FILE
    execute "$gatk_path SelectVariants -R \"$reference\" -V \"$input_vcf\" --exclude-filtered -O \"$output_vcf\"" "SelectVariants for $input_vcf"
}

VariantFiltering() {
    echo "Starting variant filtering..." | tee -a "$LOG_FILE"
    mkdir -p "${output_directory}/filter_vcf" "${output_directory}/filtered_pass_vcf"

    # Optionally learn orientation model if enabled.
    if [ "$learn_orientation" == "yes" ]; then
        echo "Learning read orientation model (FFPE artifacts)..." | tee -a $LOG_FILE
        for tar in "${output_directory}/vcf"/f1r2_*.tar.gz; do
            LearnReadOrientationModel_Filt "$tar"
        done
    fi

    # Run GetPileupSummaries for each tumor BAM (use same pattern as Mutect2)
    for tumor in "$bwa_files"/f_*recalibrated*.bam; do
        GetPileupSummaries_Filt "$tumor"
    done

    # Run CalculateContamination for each table produced.
    for table in "${output_directory}/filter_vcf"/*.table; do
        CalculateContamination_Filt "$table"
    done

    # Prepare a combined input list for filtering:
    # List VCF files from output_directory/vcf (pattern: f_*.g.vcf)
    pushd "${output_directory}/vcf" > /dev/null
    ls -1 -d "$(pwd)"/*.vcf > vcf_list.txt
    popd > /dev/null

    if [ "$learn_orientation" == "yes" ]; then
        paste "${output_directory}/vcf/vcf_list.txt" \
              <(sed 's|^.*f_|'"${output_directory}/vcf/contamination_"'|; s|\.g\.vcf$|.txt|' "${output_directory}/vcf/vcf_list.txt") \
              <(sed 's|^.*f_|'"${output_directory}/vcf/model_"'|; s|\.g\.vcf$|.tar.gz|' "${output_directory}/vcf/vcf_list.txt") \
              > combined_input_list.txt
    else
        paste "${output_directory}/vcf/vcf_list.txt" \
              <(sed 's|^.*f_|'"${output_directory}/vcf/contamination_"'|; s|\.g\.vcf$|.txt|' "${output_directory}/vcf/vcf_list.txt") \
              <(yes "" | head -n $(wc -l < "${output_directory}/vcf/vcf_list.txt")) \
              > combined_input_list.txt
    fi

    export -f FilterMutectCalls_Filt execute
    export output_directory threads reference intervals pon gr subgr gatk_path learn_orientation LOG_FILE
    parallel --jobs "$threads" --colsep '\t' FilterMutectCalls_Filt {1} {2} {3} :::: combined_input_list.txt

    export -f SelectPassedVariants_Filt execute
    parallel --jobs "$threads" SelectPassedVariants_Filt ::: $(ls "${output_directory}/filter_vcf"/filtered_*.vcf)
    rm combined_input_list.txt
    echo "Variant filtering completed." | tee -a "$LOG_FILE"
}

########################################
# Main Pipeline Function (No QC)
########################################
main_pipeline() {
    # Record overall start time.
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "\033[1mStarting Somatic Variant Calling & Filtering Pipeline.\033[0m"
    # Set up master log file and create required directories.
    LOG_FILE="${output_directory}/pipeline_$(date +'%Y%m%d_%H%M%S').log"
    mkdir -p "${output_directory}/vcf" "${output_directory}/tmp_v" \
             "${output_directory}/filter_vcf" "${output_directory}/filtered_pass_vcf"
    
    echo "Pipeline started at $(date)" >> "$LOG_FILE"
    echo "=== Part 1: Mutect2 Variant Calling ===" | tee -a "$LOG_FILE"
    Run_Mutect2_Calling

    echo "=== Part 2: Variant Filtering ===" | tee -a "$LOG_FILE"
    VariantFiltering

    # Clean up
    rm -r "${output_directory}/tmp_v"
    # Report total runtime.
    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    start_seconds=$(date -d "$start_time" +%s)
    end_seconds=$(date -d "$end_time" +%s)
    runtime_seconds=$((end_seconds - start_seconds))
    echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds" | tee -a "$LOG_FILE"
    echo -e "\033[1mSomatic Variant Calling & Filtering Pipeline completed.\033[0m"
}

########################################
# Entry Point
########################################
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

CONFIG_FILE="$1"
parse_yaml "$CONFIG_FILE"

# Ensure directory variables end with a slash.
[[ "${output_directory}" != */ ]] && output_directory="${output_directory}/"
[[ "${bwa_files}" != */ ]] && bwa_files="${bwa_files}/"

# Export necessary variables for functions.
export bwa_files output_directory threads reference intervals pon gr subgr gatk_path learn_orientation

main_pipeline
