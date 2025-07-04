#!/bin/bash

# This integrated pipeline performs:
#   0. (Optional) BWA Indexing (if ans == "yes")
#   1. Mapping with bwa-meme from input FASTQ files (paired or single-end)
#   2. Mark–duplicates (using Picard SortSam, GATK MarkDuplicatesSpark, and Picard SetNmMdAndUqTags)
#   3. Optional BQSR:
#         - Split intervals (using GATK SplitIntervals)
#         - Run BaseRecalibrator (BQSR) per interval (in parallel)
#         - Combine recalibration tables and apply BQSR (in parallel)
#   4. Run an external QC routine by calling: python3 qc_pipeline.py <combined_log> <output_directory>
#
# All settings (paths, thread counts, options, known‐sites, etc.) are provided via a YAML file.
#
# Usage:
#   ./integrated_pipeline.sh config.yaml

########################################
# YAML parser (simple key: value pairs)
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
# === BWA INDEXING (Optional) ===
########################################
BWA_Indexing() {
    echo "Indexing started" | tee -a "$LOG_FILE"
    # You can uncomment the following line if you prefer to use bwa-meme indexing.
    # bwa-meme index -a meme "$reference" -t "$threads"
    ./build_rmis_dna.sh "$reference"
    if [ $? -ne 0 ]; then
        echo "Indexing failed" | tee -a "$LOG_FILE"
        exit 1
    fi
    echo "Indexing done" | tee -a "$LOG_FILE"
}

########################################
# === MAPPING SECTION (bwa-meme) ===
########################################
Mapping() {
    echo "Starting bwa-meme mapping..." | tee -a "$LOG_FILE"
    if [ "$paired" == "yes" ]; then
        mapfile -t sample_files < "$samples"
        if (( ${#sample_files[@]} % 2 != 0 )); then
            echo "Error: The samples file must contain an even number of lines for paired-end reads." | tee -a "$LOG_FILE"
            exit 1
        fi
        local total_pairs=$(( ${#sample_files[@]} / 2 ))
        local counter=1
        for (( i=0; i<${#sample_files[@]}; i+=2 )); do
            local file1="${sample_files[$i]}"
            local file2="${sample_files[$((i+1))]}"
            echo "Mapping pair $counter of $total_pairs: $(basename "$file1") and $(basename "$file2")" | tee -a "$LOG_FILE"
            local sample_name
            sample_name=$(basename "$file1" | cut -d'_' -f1)
            execute "bwa-meme mem -7 -M -t \"$threads\" -R \"@RG\tID:$sample_name\tSM:$sample_name\tPL:ILLUMINA\" \"$reference\" \"$file1\" \"$file2\" > \"${mapping_output_dir}${sample_name}.sam\"" "bwa-meme mapping for $sample_name"
            execute "samtools view -@ \"$threads\" -b \"${mapping_output_dir}${sample_name}.sam\" > \"${mapping_output_dir}${sample_name}.bam\"" "samtools view for $sample_name"
            execute "samtools sort -@ \"$threads\" \"${mapping_output_dir}${sample_name}.bam\" -o \"${mapping_output_dir}sorted_${sample_name}.bam\"" "samtools sort for $sample_name"
            execute "samtools flagstat -@ \"$threads\" \"${mapping_output_dir}sorted_${sample_name}.bam\" > \"${mapping_output_dir}LOG_${sample_name}.txt\"" "samtools flagstat for $sample_name"
            rm -f "${mapping_output_dir}${sample_name}.sam" "${mapping_output_dir}${sample_name}.bam"
            counter=$((counter + 1))
        done
    else
        while IFS= read -r file; do
            echo "Mapping: $(basename "$file")" | tee -a "$LOG_FILE"
            local sample_name
            sample_name=$(basename "$file")
            execute "bwa-meme -M -t \"$threads\" -R \"@RG\tID:$sample_name\tSM:$sample_name\tPL:ILLUMINA\" \"$reference\" \"$file\" > \"${mapping_output_dir}${sample_name}.sam\"" "bwa-meme mapping for $sample_name"
            execute "samtools view -@ \"$threads\" -b \"${mapping_output_dir}${sample_name}.sam\" > \"${mapping_output_dir}${sample_name}.bam\"" "samtools view for $sample_name"
            execute "samtools sort -@ \"$threads\" \"${mapping_output_dir}${sample_name}.bam\" -o \"${mapping_output_dir}sorted_${sample_name}.bam\"" "samtools sort for $sample_name"
            execute "samtools flagstat -@ \"$threads\" \"${mapping_output_dir}sorted_${sample_name}.bam\" > \"${mapping_output_dir}LOG_${sample_name}.txt\"" "samtools flagstat for $sample_name"
            rm -f "${mapping_output_dir}${sample_name}.sam" "${mapping_output_dir}${sample_name}.bam"
        done < "$samples"
    fi
    echo "Mapping completed." | tee -a "$LOG_FILE"
}
########################################
# === MARK DUPLICATES SECTION ===
########################################
########################################
# === MARK DUPLICATES SECTION ===
########################################

QuerySort() {
    local input_bam="$1"
    local base_name
    base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local output_bam="${mapping_output_dir}qs_${base_name}"
    local log_file="${output_directory}/logs_mark_dup/qs_${base_name}.log"

    echo "Query Sorting $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"
    # Set the log file for the execute function
    LOG_FILE="$log_file"
    execute "java -jar \"$picard_path\" SortSam I=\"$input_bam\" O=\"$output_bam\" SORT_ORDER=queryname TMP_DIR=\"${output_directory}/qs_tmp\"" "QuerySort for $(basename "$input_bam")"
}

MarkDuplicates() {
    local input_bam="$1"
    local base_name
    base_name=$(basename "$input_bam" | sed 's/^qs_//')
    local output_bam="${mapping_output_dir}md_${base_name}"
    local metrics_file="${mapping_output_dir}marked_dup_metrics_${base_name%.bam}.txt"
    local log_file="${output_directory}/logs_mark_dup/md_${base_name}.log"

    echo "Marking duplicates in $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"
    LOG_FILE="$log_file"
    execute "\"$gatk_path\" MarkDuplicatesSpark -I \"$input_bam\" -O \"$output_bam\" -M \"$metrics_file\" --tmp-dir \"${output_directory}/md_tmp\" --create-output-bam-splitting-index true --spark-runner LOCAL --conf \"spark.executor.cores=${threads}\"" "MarkDuplicatesSpark for $(basename "$input_bam")"
}

SetNmMdAndUqTags() {
    local input_bam="$1"
    local base_name
    base_name=$(basename "$input_bam" | sed 's/^md_//')
    local output_bam="${mapping_output_dir}f_${base_name}"
    local log_file="${output_directory}/logs_mark_dup/f_${base_name}.log"

    echo "Setting NM, MD, and UQ tags in $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"
    LOG_FILE="$log_file"
    execute "java -jar \"$picard_path\" SetNmMdAndUqTags R=\"$reference\" I=\"$input_bam\" O=\"$output_bam\" TMP_DIR=\"${output_directory}/fix_tmp\" CREATE_INDEX=true" "SetNmMdAndUqTags for $(basename "$input_bam")"
}

########################################
# === BQSR SECTION (Optional) ===
########################################

SplitIntervals() {
    local log_file="${output_directory}/logs_bqsr/split_intervals.log"
    echo "Splitting intervals..." | tee "$log_file"
    LOG_FILE="$log_file"
    execute "\"$gatk_path\" SplitIntervals -R \"$reference\" -L \"$intervals\" --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION --scatter-count \"$threads\" -O intervals" "SplitIntervals"
    echo "Interval splitting completed." | tee -a "$log_file"
}

BQSR() {
    local input_bam="$1"
    local interval_split="$2"
    local base_name
    base_name=$(basename "$input_bam" .bam)
    local interval_code
    interval_code=$(basename "$interval_split" | sed 's/\.interval_list$//')
    local recal_table="${base_name}_${interval_code}.table"
    local log_file="${output_directory}/logs_bqsr/bqsr_${base_name}_${interval_code}.log"

    # Process known_sites: split the comma-separated string and add "--known-sites" before each site.
    local known_sites_cmd=""
    IFS=',' read -ra sites <<< "$known_sites"
    for site in "${sites[@]}"; do
        site=$(echo "$site" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
        known_sites_cmd+=" --known-sites $site"
    done

    echo "Running BQSR on $(basename "$input_bam") for interval $interval_code" | tee "$log_file"
    LOG_FILE="$log_file"
    execute "\"$gatk_path\" BaseRecalibrator -I \"$input_bam\" -R \"$reference\" -L \"$interval_split\" --interval-padding \"$padding\" $known_sites_cmd -O \"$recal_table\" --tmp-dir \"${output_directory}/bqsr_tmp\"" "BaseRecalibrator for $(basename "$input_bam") on interval $interval_code"
}

CombineTables() {
    local base_name="$1"
    shift
    local recal_tables=("$@")
    local output_combined_table="${mapping_output_dir}${base_name}_combined.table"
    local log_file="${output_directory}/logs_bqsr/combine_${base_name}.log"
    local arguments_file="${output_directory}/logs_bqsr/${base_name}_recal_tables.args"

    echo "Combining BQSR tables for $base_name" | tee "$log_file"
    # Create or clear the arguments file
    : > "$arguments_file"
    for table in "${recal_tables[@]}"; do
        echo "-I $table" >> "$arguments_file"
    done
    LOG_FILE="$log_file"
    execute "\"$gatk_path\" GatherBQSRReports --arguments_file \"$arguments_file\" -O \"$output_combined_table\"" "GatherBQSRReports for $base_name"
    echo "Combined tables into $output_combined_table" | tee -a "$log_file"
    rm "$arguments_file"
}

ApplyBQSR() {
    local input_bam="$1"
    local base_name
    base_name=$(basename "$input_bam" .bam)
    local combined_table="${mapping_output_dir}${base_name}_combined.table"
    local output_bqsr_bam="${mapping_output_dir}${base_name}_recalibrated.bam"
    local log_file="${output_directory}/logs_bqsr/applybqsr_${base_name}.log"

    echo "Applying BQSR to $(basename "$input_bam") using table $combined_table" | tee "$log_file"
    if [ ! -f "$combined_table" ]; then
        echo "Error: Combined table $combined_table not found!" | tee -a "$log_file"
        exit 1
    fi
    LOG_FILE="$log_file"
    execute "\"$gatk_path\" ApplyBQSR -R \"$reference\" -I \"$input_bam\" --bqsr-recal-file \"$combined_table\" -O \"$output_bqsr_bam\"" "ApplyBQSR for $(basename "$input_bam")"
    echo "Applied BQSR to $(basename "$input_bam"); output: $(basename "$output_bqsr_bam")" | tee -a "$log_file"
}

########################################
# === QC Routine Call (External) ===
########################################
run_qc() {
    echo "Running QC routine..." | tee -a "$LOG_FILE"
    more "$mapping_output_dir"LOG*.txt > "${output_directory}/combine.txt"
    execute "python3 qc_pipeline.py ${output_directory}/combine.txt ${output_directory}" "Python QC routine"
}

########################################
# === MAIN PIPELINE FUNCTION ===
########################################
main_pipeline() {
    # Record the start time.
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "\033[1mStarting Integrated Pipeline.\033[0m"
    
    # Set up master log file and required output directories.
    LOG_FILE="${output_directory}/pipeline_$(date +'%Y%m%d_%H%M%S').log"
    mkdir -p "${output_directory}/logs_mark_dup" "${output_directory}/logs_bqsr" \
             "${output_directory}/qs_tmp" "${output_directory}/md_tmp" "${output_directory}/fix_tmp" \
             "${output_directory}/bqsr_tmp"
    echo "Pipeline started at $(date)" >> "$LOG_FILE"
    
    # Define mapping_output_dir as a subdirectory under output_directory.
    mapping_output_dir="${output_directory}bwa_files/"
    mkdir -p "$mapping_output_dir"
    
    # --- Optional: BWA Indexing ---
    if [ "$ans" == "yes" ]; then
        echo "=== BWA Indexing ===" | tee -a "$LOG_FILE"
        BWA_Indexing
    fi
    
    # --- Part 1: Mapping ---
    Mapping
    
    # --- Part 2: MarkDuplicates ---
    if [ "$markduplicates_enabled" == "yes" ]; then
        echo "=== MarkDuplicates Section ===" | tee -a "$LOG_FILE"
        pushd "$mapping_output_dir" > /dev/null
        ls -1 -d "$(pwd)"/*.bam > bam_list.txt
        popd > /dev/null
        echo "Query Sorting BAMs..." | tee -a "$LOG_FILE"
        while read -r bam; do
            QuerySort "$bam"
        done < "${mapping_output_dir}bam_list.txt"
        pushd "$mapping_output_dir" > /dev/null
        ls -1 -d "$(pwd)"/qs*.bam > qs_list.txt
        popd > /dev/null
        echo "Marking Duplicates..." | tee -a "$LOG_FILE"
        while read -r qs_bam; do
            MarkDuplicates "$qs_bam"
        done < "${mapping_output_dir}qs_list.txt"
	pushd "$mapping_output_dir" > /dev/null
	ls -1 -d "$(pwd)"/md*.bam > md_list.txt
	popd > /dev/null
	
        echo "Fixing BAMs (SetNmMdAndUqTags)..." | tee -a "$LOG_FILE"
        export -f SetNmMdAndUqTags
	export output_directory threads reference picard_path gatk_path intervals padding known_sites mapping_output_dir # Probalby shoul put in the beggingn 
        parallel --jobs "$threads" SetNmMdAndUqTags :::: "${mapping_output_dir}md_list.txt"
        mv "${mapping_output_dir}bam_list.txt" "${mapping_output_dir}md_list.txt" "${mapping_output_dir}qs_list.txt" "${output_directory}/logs_mark_dup/"
    else
        echo "MarkDuplicates section is disabled." | tee -a "$LOG_FILE"
    fi
    
    # --- Part 3: BQSR Section (Optional) ---
    if [ "$bqsr_enabled" == "yes" ]; then
	export -f BQSR
	export output_directory threads reference picard_path gatk_path intervals padding known_sites mapping_output_dir
        if [ -z "$known_sites" ]; then
            echo "BQSR is enabled but no known_sites provided. Exiting." | tee -a "$LOG_FILE"
            exit 1
        fi
        echo "=== BQSR Section ===" | tee -a "$LOG_FILE"
        mkdir -p intervals
        mkdir -p "${output_directory}/logs_bqsr"
        SplitIntervals
        for bam_file in "${mapping_output_dir}"/f_*.bam; do
            echo "Running BQSR for $(basename "$bam_file")" | tee -a "$LOG_FILE"
            parallel --jobs "$threads" BQSR ::: "$bam_file" ::: intervals/*.interval_list
        done
        for bam_file in "${mapping_output_dir}"/f_*.bam; do
            tables=()
            base_name=$(basename "$bam_file" .bam)
            for interval_split in intervals/*.interval_list; do
                interval_code=$(basename "$interval_split" | sed 's/\.interval_list$//')
                tables+=("${base_name}_${interval_code}.table")
            done
            CombineTables "$base_name" "${tables[@]}"
        done
        echo "Applying BQSR in parallel..." | tee -a "$LOG_FILE"
        export -f ApplyBQSR
        parallel --jobs "$threads" ApplyBQSR ::: "${mapping_output_dir}"/f_*.bam
        rm -r "${output_directory}/bqsr_tmp" intervals
        rm f_*scattered.table 2>/dev/null
    else
        echo "BQSR section is disabled." | tee -a "$LOG_FILE"
    fi
    
    # --- Part 4: QC ---
    run_qc
    
    # Report total runtime.
    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    start_seconds=$(date -d "$start_time" +%s)
    end_seconds=$(date -d "$end_time" +%s)
    runtime_seconds=$((end_seconds - start_seconds))
    echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds" | tee -a "$LOG_FILE"
    echo "Pipeline completed." | tee -a "$LOG_FILE"
}


########################################
# === Entry Point ===
########################################
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

CONFIG_FILE="$1"
parse_yaml "$CONFIG_FILE"

# Ensure directory variables end with a slash.
[[ "${output_directory}" != */ ]] && output_directory="${output_directory}/"
[[ "${samples}" != */ ]] && samples="${samples}"

main_pipeline
