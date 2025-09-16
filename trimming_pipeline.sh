#!/bin/bash

# It performs:
#   1. FastQC on raw samples (all at once to utilize parallel processing).
#   2. MultiQC to summarize FastQC results.
#   3. Trimming with fastp (paired-end or single-end based on config).
#   4. FastQC on trimmed samples (all at once).
#   5. MultiQC to summarize post-trimming FastQC results.
#
# The configuration is provided via a YAML file with keys:
#   samples, output_directory, read_type, and threads.
#
# Example YAML (config.yaml):
#
#   samples: /home/kousis/work/meth/01.RawData/trimm_full.txt
#   output_directory: /home/kousis/work/meth/output
#   read_type: paired
#   threads: 4
#

# --- Helper: Simple YAML parser ---
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
# --- Function to run FastQC in parallel on multiple files ---
run_fastqc_parallel() {
    local out_dir="$1"
    shift
    local files=("$@")
    echo "Running FastQC on files: ${files[*]}" | tee -a "$LOG_FILE"
    # Call fastqc with the provided threads and all files as input.
    execute "fastqc -t \"$threads\" ${files[*]} -o \"$out_dir\"" "FastQC on multiple files"
}

# --- Function to run MultiQC ---
run_multiqc() {
    local qc_dir="$1"
    echo "Running MultiQC on $qc_dir..." | tee -a "$LOG_FILE"
    execute "multiqc \"$qc_dir\" -o \"$qc_dir\"" "MultiQC on $qc_dir"
}

# --- Functions for fastp trimming ---
run_fastp_paired() {
    local infile1="$1"
    local infile2="$2"
    local base1 base2
    base1=$(basename "$infile1")
    base2=$(basename "$infile2")
    local outfile1="${TRIMMED_DIR}/${base1%.*}_trim.${base1##*.}"
    local outfile2="${TRIMMED_DIR}/${base2%.*}_trim.${base2##*.}"
    echo "Running fastp on paired files: $infile1 and $infile2" | tee -a "$LOG_FILE"
    execute "fastp -w \"$threads\" -i \"$infile1\" -I \"$infile2\" -o \"$outfile1\" -O \"$outfile2\"" "fastp on paired files $infile1 and $infile2"
}

run_fastp_single() {
    local infile="$1"
    local base
    base=$(basename "$infile")
    local outfile="${TRIMMED_DIR}/${base%.*}_trim.${base##*.}"
    echo "Running fastp on single file: $infile" | tee -a "$LOG_FILE"
    execute "fastp -w \"$threads\" -i \"$infile\" -o \"$outfile\"" "fastp on single file $infile"
}

# --- Main Workflow ---
main() {
    if [ "$#" -ne 1 ]; then
        echo "Usage: $0 <config.yaml>"
        exit 1
    fi

    CONFIG_FILE="$1"
    parse_yaml "$CONFIG_FILE"

    # Check required parameters; default threads to 1 if not provided.
    if [ -z "$samples" ] || [ -z "$output_directory" ] || [ -z "$read_type" ]; then
        echo "Error: The YAML file must provide samples, output_directory, and read_type."
        exit 1
    fi

    if [ -z "$threads" ]; then
        threads=1
    fi

    # Create output directories.
    FASTQC_PRE_DIR="${output_directory}/fastqc_pre"
    FASTQC_POST_DIR="${output_directory}/fastqc_post"
    TRIMMED_DIR="${output_directory}/trimmed"
    LOG_DIR="${output_directory}/Trimming_logs"

    mkdir -p "$FASTQC_PRE_DIR" "$FASTQC_POST_DIR" "$TRIMMED_DIR" "$LOG_DIR"

    # Create a log file with a timestamp.
    LOG_FILE="${LOG_DIR}/pipeline_$(date +'%Y%m%d_%H%M%S').log"
    echo "Pipeline started at $(date)" | tee -a "$LOG_FILE"

    # Read sample file list.
    mapfile -t raw_files < "$samples"
    if [ "$read_type" = "paired" ]; then
        if (( ${#raw_files[@]} % 2 != 0 )); then
            echo "Error: Number of files in $samples is not even for paired-end reads." | tee -a "$LOG_FILE"
            exit 1
        fi
    fi

    echo "Starting pre-trimming FastQC..." | tee -a "$LOG_FILE"
    # Run FastQC on all raw sample files in one go.
    run_fastqc_parallel "$FASTQC_PRE_DIR" "${raw_files[@]}"

    run_multiqc "$FASTQC_PRE_DIR"

    echo "Starting trimming step with fastp..." | tee -a "$LOG_FILE"
    if [ "$read_type" = "paired" ]; then
        # Process paired files two at a time.
        for ((i = 0; i < ${#raw_files[@]}; i += 2)); do
            run_fastp_paired "${raw_files[i]}" "${raw_files[i+1]}"
        done
    else
        for file in "${raw_files[@]}"; do
            run_fastp_single "$file"
        done
    fi

    # Gather all trimmed files.
    mapfile -t trimmed_files < <(find "$TRIMMED_DIR" -maxdepth 1 -type f -name "*_trim.*")
    echo "Starting post-trimming FastQC..." | tee -a "$LOG_FILE"
    run_fastqc_parallel "$FASTQC_POST_DIR" "${trimmed_files[@]}"

    run_multiqc "$FASTQC_POST_DIR"

    echo "Pipeline completed successfully at $(date)" | tee -a "$LOG_FILE"
}

# Call main function with all provided arguments.
main "$@"
