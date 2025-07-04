#!/bin/bash
# vcf_annotation_pipeline.sh
#
# This pipeline annotates VCF files using Funcotator and, if enabled,
# further converts annotated VCFs to MAFs, annotates them with OncoKB,
# and finally merges the resulting MAFs with hotspot and Cosmic data
# that can also be fitlered based on known databases
#
# It reads configuration parameters from a YAML file. Expected keys include:
#
#   vcf_dir:                Directory containing VCF files to annotate
#   output_directory:       Base directory for all outputs (logs, annotated MAFs, etc.)
#   threads:                Number of parallel jobs to run
#   reference:              Reference genome FASTA file
#   funcotator_data_sources: Path to Funcotator data sources
#
#   # OncoKB/MAF conversion parameters (if oncokb_enabled is "yes")
#   oncokb_enabled:         yes/no – whether to perform the additional OncoKB/MAF conversion
#   onko_out_dir:           (Optional) Directory where OncoKB-converted MAFs will be saved.
#                           If not provided, it defaults to ${output_directory}/maf_onko
#   maf2vcf:                Path to the vcf2maf directory (which contains vcf2maf.pl)
#   vep_anot:               Path to the reference FASTA used by VEP (for vcf2maf)
#   vep_path:               Path to the VEP installation (or its bin directory)
#
#   # OncoKB annotation parameters:
#   oncokb_path:            Path to OncoKB annotator directory (contains MafAnnotator.py)
#   token:                  API token to use with OncoKB annotator
#
#   # Merging MAFs parameters:
#   merge_output_folder:    Folder where the final merged MAFs will be written
#   merge_mafs_enabled:     yes/no – whether to run the merging step
#
#   # Optional: Hotspot and Cosmic files (for additional annotation)
#   hotspot_file:           Path to a hotspot file (optional)
#   cosmic_file:            Path to a Cosmic TSV file (optional)
#
#   # Additional mapping file for OncoKB conversion
#   map_list:               Chromosome mapping file for bcftools renaming
#
# Usage:
#   ./vcf_annotation_pipeline.sh config.yaml

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
execute() {
    local cmd="$1"
    local description="$2"

    # Extract the tool name (first token of the command)
    local tool
    tool=$(echo "$cmd" | awk '{print $1}')

    # Check if the tool exists
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: Required tool '$tool' not found. Please install it before proceeding." | tee -a "$LOG_FILE"
        exit 127
    fi

    # Log the command to the log file (but not to the terminal)
    echo "Executing: $description" >> "$LOG_FILE"
    echo "Command: $cmd" >> "$LOG_FILE"

    # Execute the command silently and log both stdout and stderr
    eval "$cmd" >> "$LOG_FILE" 2>&1
    local status=$?
    if [ $status -ne 0 ]; then
        echo "Error: '$description' failed. See log: $LOG_FILE" >> "$LOG_FILE"
        exit $status
    fi
}

########################################
# Part 1: VCF Annotation with Funcotator
########################################
AnnotateVCF() {
    local input_vcf="$1"
    local base_name
    base_name=$(basename "$input_vcf" .vcf)
    local output_vcf="${output_directory}/maf_func/annotated_${base_name}.maf"

    echo "Annotating $input_vcf -> $output_vcf" | tee -a $LOG_FILE
    local cmd="gatk Funcotator \
        -R \"$reference\" \
        -V \"$input_vcf\" \
        -O \"$output_vcf\" \
        --output-file-format MAF \
        --data-sources-path \"$funcotator_data_sources\" \
        --ref-version hg38"
    execute "$cmd" "Funcotator annotation for $input_vcf"
}

########################################
# Part 2: OncoKB/MAF Conversion (Optional)
########################################
OncoKB_Conversion() {
    # Input VCF file passed as the first argument.
    local input_vcf="$1"
    
    # Derive base filename by removing .vcf or .vcf.gz extension.
    local base_name
    base_name=$(basename "$input_vcf" | sed 's/\.vcf$//; s/\.vcf\.gz$//')
    
    echo "Processing: $base_name" | tee -a "$LOG_FILE"

    # Set the output directory if not defined.
    if [ -z "$onko_out_dir" ]; then
        onko_out_dir="${output_directory}/maf_onko"
    fi
    mkdir -p "$onko_out_dir"
    mkdir -p "$onko_out_dir/vep"

    # Rename chromosomes using bcftools with the mapping file.
    local renamed_vcf="${onko_out_dir}/kb_${base_name}.vcf"
    execute "bcftools annotate --rename-chrs \"$map_list\" \"$input_vcf\" -Ov -o \"$renamed_vcf\"" "Chromosome renaming for $base_name"

    # Ensure conda is properly set up, then deactivate any current environment.
    eval "$(conda shell.bash hook)"
    conda deactivate

    # Convert VCF to MAF using vcf2maf.
    local maf_file="${onko_out_dir}/kb_${base_name}.maf"
    execute "perl \"${maf2vcf}/vcf2maf.pl\" \
        --ref-fasta \"$vep_anot\" \
        --input-vcf \"$renamed_vcf\" \
        --output-maf \"$maf_file\" \
        --vep-path \"$vep_path\" \
        --vep-data \"$vep_data\" \
        --ncbi-build GRCh38" "VCF to MAF conversion for $base_name"

    # Activate the TOSV conda environment.
    execute "conda activate TOSV" "Activate TOSV conda environment"
    conda activate TOSV

    # Run the OncoKB MafAnnotator Python script.
    execute "python3 \"$oncokb_path/MafAnnotator.py\" -i \"$maf_file\" -o \"${onko_out_dir}/discr_KB_${base_name}.maf\" -b \"$token\" -d -t \"$tissue\"" "Running MafAnnotator on $base_name"

    # Move the original MAF file to the vep subdirectory.
    mv "$maf_file" "$onko_out_dir/vep"
}

########################################
# Part 3: Merge MAFs with Hotspot and Cosmic Annotation (Optional)
########################################
Merge_MAFs() {
    echo "Merging MAF files with hotspot and Cosmic annotation..." | tee -a "$LOG_FILE"
    # Use Funcotator output from maf_t and OncoKB conversion output from kb_out.
    execute "python3 ${scripts}/merge_mafs.py --func_dir \"${output_directory}/maf_func\" --onko_dir \"$onko_out_dir\" --output_dir \"$merge_output_folder\" --hotspot_file \"$hotspot_file\" --num_threads \"$threads\" --cosmic_file \"$cosmic_file\"" "Merge MAFs with hotspot and Cosmic annotation"
}

########################################
# Main Pipeline Function
########################################
main_pipeline() {
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "\033[1mStarting VCF Annotation & MAF Conversion Pipeline.\033[0m"

    # Set up a master log file and create necessary directories.
    LOG_FILE="${output_directory}/pipeline_$(date +'%Y%m%d_%H%M%S').log"
    mkdir -p "${output_directory}/maf_func"
    echo "Pipeline started at $(date)" >> "$LOG_FILE"

    # --- Part 1: VCF Annotation with Funcotator ---
    # echo "=== Part 1: VCF Annotation with Funcotator ===" | tee -a "$LOG_FILE"
    # export -f AnnotateVCF execute
    # find "$vcf_dir" -name "*.vcf" | parallel --jobs "$threads" AnnotateVCF

    # # --- Part 2: Optional OncoKB/MAF Conversion ---
    # if [ "$oncokb_enabled" == "yes" ]; then
    #     echo "=== Part 2: OncoKB/MAF Conversion ===" | tee -a "$LOG_FILE"
    #     find "$vcf_dir" -name "*.vcf" | parallel --jobs "$threads" OncoKB_Conversion

    # else
    #     echo "OncoKB/MAF Conversion is disabled." | tee -a "$LOG_FILE"
    # fi

    # # --- Part 3: Optional MAF Merging with Hotspot and Cosmic Annotation ---
    # if [ "$merge_mafs_enabled" == "yes" ]; then
    #     echo "=== Part 3: Merging MAFs with Hotspot and Cosmic Annotation ===" | tee -a "$LOG_FILE"
    #     Merge_MAFs
    # else
    #     echo "MAF merging is disabled." | tee -a "$LOG_FILE"
    # fi

    # --- Part 4: Optional AF Filtering based on known db ---
    if [ "$AF_filtering_enabled" == "yes" ]; then
	echo "=== Part 4: AF Filtering ===" | tee -a "$LOG_FILE"
	known_sites_cmd=""
	IFS=',' read -ra sites <<< "$known_db"
	for site in "${sites[@]}"; do
            site=$(echo "$site" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
            known_sites_cmd+=" --vcf-file $site"
	done

	execute "python3 \"$scripts/filter_af.py\" --maf-input \"${merge_output_folder}/combined_merged_with_cosmic.maf\" $known_sites_cmd --output \"${merge_output_folder}/af_filtered.maf\"" "Running MafAnnotator on $base_name"
    else
	echo "AF filtering is disabled." | tee -a "$LOG_FILE"
    fi

    
    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    start_seconds=$(date -d "$start_time" +%s)
    end_seconds=$(date -d "$end_time" +%s)
    runtime_seconds=$((end_seconds - start_seconds))
    echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds" | tee -a "$LOG_FILE"
    echo -e "\033[1mPipeline completed.\033[0m"
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
[[ "${vcf_dir}" != */ ]] && vcf_dir="${vcf_dir}/"

# If onko_out_dir is not provided in the YAML, default it to output_directory/kb_out.
if [ -z "$onko_out_dir" ]; then
    onko_out_dir="${output_directory}/maf_onko"
fi

# Export variables and functions needed for parallel execution.
export vcf_dir output_directory threads reference funcotator_data_sources oncokb_enabled onko_out_dir maf2vcf vep_anot vep_path oncokb_path token LOG_FILE tissue vep_data vep_path onko_out_dir map_list maf2vcf
export merge_output_folder gatk_path hotspot_file cosmic_file merge_mafs_enabled map_list
export -f AnnotateVCF OncoKB_Conversion Merge_MAFs execute

main_pipeline
