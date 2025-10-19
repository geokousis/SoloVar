#!/bin/bash
export SHELL=/bin/bash

# ====================================
# Integrated Mapping Pipeline
# ====================================
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
#   ./mapping_pipeline.sh config.yaml

########################################
# Simple YAML parser (key: value pairs)
########################################
parse_yaml() {
    local yaml_file="$1"
    while IFS=":" read -r key value; do
        key="$(echo "$key" | xargs)"
        value="$(echo "$value" | xargs)"
        if [[ -n "$key" && ! "$key" =~ ^# ]]; then
            eval "${key}='${value}'"
        fi
    done < "$yaml_file"
}

########################################
# Global Settings
########################################
DRY_RUN=false
VERBOSE_MODE=silent

########################################
# Execute Command Helper
########################################
execute() {
    local cmd="$1"
    local desc="$2"
    local runtime_flag="${3:-}"
    local verbose_flag="${4:-$VERBOSE_MODE}"

    if [[ "$cmd" =~ ^[[:space:]]*\(\( ]]; then
        : # allow bash arithmetic expressions
    else
        local tool
        tool=$(echo "$cmd" | awk '{print $1}')
        if ! command -v "$tool" &>/dev/null; then
            echo -e "\033[1mError:\033[0m Required tool '$tool' not found."
            exit 127
        fi
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo -e "\033[2m[DRY RUN]\033[0m \033[1m$desc\033[0m"
        echo "[CMD] $cmd"
        return 0
    fi

    echo -e "\033[1m[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $desc\033[0m" | tee -a "$LOG_FILE"
    echo "[CMD] $cmd" >> "$LOG_FILE"

    local start_time end_time runtime
    if [[ "$runtime_flag" == "runtime" ]]; then start_time=$(date +%s); fi

    if [[ "$verbose_flag" == "verbose" ]]; then
        ( eval "$cmd" ) 2>&1 | tee -a "$LOG_FILE"
        status=${PIPESTATUS[0]}
    else
        eval "$cmd" >> "$LOG_FILE" 2>&1
        status=$?
    fi

    if [[ "$runtime_flag" == "runtime" ]]; then
        end_time=$(date +%s)
        runtime=$((end_time - start_time))
        echo -e "\033[1m[INFO]\033[0m Command completed in $runtime seconds." | tee -a "$LOG_FILE"
    fi

    if [[ $status -ne 0 ]]; then
        echo -e "\033[1mError:\033[0m Command failed (exit $status). See $LOG_FILE" >&2
        exit $status
    fi
}

########################################
# Cleanup Handler
########################################
cleanup() {
    echo -e "\033[2m[INFO]\033[0m Exiting at $(date '+%Y-%m-%d %H:%M:%S')"
}
trap cleanup EXIT

########################################
# BWA Indexing (Optional)
########################################
BWA_Indexing() {
    execute "./build_rmis_dna.sh '$reference'" "Running BWA indexing"
}

########################################
# Mapping Function
########################################
Mapping() {
    execute "mkdir -p '$mapping_output_dir'" "Creating mapping output directory"

    if [[ "$paired" == "yes" ]]; then
        mapfile -t samples_list < "$samples"
        execute "(( \${#samples_list[@]} % 2 == 0 ))" "Checking paired-end file count"

        for ((i=0; i<${#samples_list[@]}; i+=2)); do
            f1=${samples_list[i]}
            f2=${samples_list[i+1]}
            sample_name=$(basename "$f1" | cut -d'_' -f1)

            execute "bwa-meme mem -M -t '$threads' -R '@RG\\tID:$sample_name\\tSM:$sample_name\\tPL:ILLUMINA' '$reference' '$f1' '$f2' > '$mapping_output_dir/${sample_name}.sam'" "Mapping sample $sample_name"
            execute "samtools view -b -@ '$threads' '$mapping_output_dir/${sample_name}.sam' > '$mapping_output_dir/${sample_name}.bam'" "Converting SAM to BAM for $sample_name"
            execute "samtools sort -@ '$threads' '$mapping_output_dir/${sample_name}.bam' -o '$mapping_output_dir/sorted_${sample_name}.bam'" "Sorting BAM for $sample_name"
            execute "samtools flagstat -@ '$threads' '$mapping_output_dir/sorted_${sample_name}.bam' > '$mapping_output_dir/LOG_${sample_name}.txt'" "Flagstat report for $sample_name"
            execute "rm '$mapping_output_dir/${sample_name}.sam' '$mapping_output_dir/${sample_name}.bam'" "Cleaning intermediate files for $sample_name"
        done

    else
        while read -r f; do
            sample_name=$(basename "$f" | cut -d'_' -f1)

            execute "bwa-meme mem -M -t '$threads' -R '@RG\\tID:$sample_name\\tSM:$sample_name\\tPL:ILLUMINA' '$reference' '$f' > '$mapping_output_dir/${sample_name}.sam'" "Mapping sample $sample_name (single-end)"
            execute "samtools view -b -@ '$threads' '$mapping_output_dir/${sample_name}.sam' > '$mapping_output_dir/${sample_name}.bam'" "Converting SAM to BAM for $sample_name"
            execute "samtools sort -@ '$threads' '$mapping_output_dir/${sample_name}.bam' -o '$mapping_output_dir/sorted_${sample_name}.bam'" "Sorting BAM for $sample_name"
            execute "samtools flagstat -@ '$threads' '$mapping_output_dir/sorted_${sample_name}.bam' > '$mapping_output_dir/LOG_${sample_name}.txt'" "Flagstat report for $sample_name"
            execute "rm '$mapping_output_dir/${sample_name}.sam' '$mapping_output_dir/${sample_name}.bam'" "Cleaning intermediate files for $sample_name"
        done < "$samples"
    fi
}
########################################
# Mark Duplicates
########################################
QuerySort() {
    local bam="$1"
    local output_bam="${mapping_output_dir}qs_$(basename "$bam")"
    execute "java -jar '$picard_path' SortSam I='$bam' O='$output_bam' SORT_ORDER=queryname TMP_DIR='${output_directory}/qs_tmp'" "Sorting BAM $bam by queryname"
}

MarkDuplicates() {
    local qs_bam="$1"
    local output_bam="${mapping_output_dir}md_$(basename "$qs_bam")"
    local metrics_file="${mapping_output_dir}metrics_$(basename "$qs_bam" .bam).txt"
    execute "$gatk_path MarkDuplicatesSpark -I '$qs_bam' -O '$output_bam' -M '$metrics_file' --tmp-dir '${output_directory}/md_tmp'" "Marking duplicates in $qs_bam"
}

SetNmMdAndUqTags() {
    local md_bam="$1"
    local output_bam="${mapping_output_dir}f_$(basename "$md_bam")"
    execute "java -jar '$picard_path' SetNmMdAndUqTags R='$reference' I='$md_bam' O='$output_bam' TMP_DIR='${output_directory}/fix_tmp' CREATE_INDEX=true" "Fixing tags in $md_bam"
}
########################################
# BQSR (Improved)
########################################
SplitIntervals() {
    echo "Splitting intervals..." | tee -a "$LOG_FILE"
    local out_int="${output_directory%/}/intervals"
    mkdir -p "$out_int"
    execute "$GATK_CMD SplitIntervals \
      -R \"$reference\" \
      ${intervals:+-L \"$intervals\"} \
      --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
      --scatter-count \"$threads\" \
      -O \"$out_int\"" "SplitIntervals"
    echo "Interval splitting completed." | tee -a "$LOG_FILE"
}
CombineTables() {
    local base_name="$1"; shift
    local recal_tables=("$@")
    local output_combined_table="${mapping_output_dir%/}/${base_name}_combined.table"
    local tmp_dir="${output_directory%/}/tmp_bqsr"
    local arguments_file="${tmp_dir}/${base_name}_recal_tables.args"
    echo "Combining BQSR tables for $base_name" | tee -a "$LOG_FILE"
    mkdir -p "$tmp_dir"
    : > "$arguments_file"
    for t in "${recal_tables[@]}"; do
        echo "-I $t" >> "$arguments_file"
    done
    execute "$GATK_CMD GatherBQSRReports --arguments_file \"$arguments_file\" -O \"$output_combined_table\"" "GatherBQSRReports for $base_name"
    echo "Combined tables into $output_combined_table" | tee -a "$LOG_FILE"
    rm -f "$arguments_file"
}
ApplyBQSR() {
    local input_bam="$1"
    local base_name; base_name=$(basename "$input_bam" .bam)
    local combined_table="${mapping_output_dir%/}/${base_name}_combined.table"
    local output_bqsr_bam="${mapping_output_dir%/}/${base_name}_recalibrated.bam"
    echo "Applying BQSR to $(basename "$input_bam") using table $combined_table" | tee -a "$LOG_FILE"
    if [[ ! -f "$combined_table" ]]; then
        echo "Error: Combined table $combined_table not found!" | tee -a "$LOG_FILE"
        return 1
    fi
    execute "$GATK_CMD ApplyBQSR -R \"$reference\" -I \"$input_bam\" --bqsr-recal-file \"$combined_table\" -O \"$output_bqsr_bam\"" "ApplyBQSR for $(basename "$input_bam")"
    echo "Applied BQSR to $(basename "$input_bam"); output: $(basename "$output_bqsr_bam")" | tee -a "$LOG_FILE"
}
run_qc() {
    execute "cat '${mapping_output_dir}/LOG_'*.txt > '${output_directory}/combined_log.txt'" "Combining flagstat logs"
    execute "python3 qc_pipeline.py '${output_directory}/combined_log.txt' '${output_directory}'" "Running QC Python script"
}
BQSR() {
    local input_bam="$1"
    local interval_split="$2"
    local base_name; base_name=$(basename "$input_bam" .bam)
    local interval_code; interval_code=$(basename "$interval_split" | sed 's/\.interval_list$//')
    local recal_table="${base_name}_${interval_code}.table"

    local known_sites_cmd=""
    if declare -p known_sites 2>/dev/null | grep -q 'declare \-a'; then
        for site in "${known_sites[@]}"; do
            site=$(echo "$site" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
            [[ -n "$site" ]] && known_sites_cmd+=" --known-sites $site"
        done
    else
        IFS=',' read -ra sites <<< "${known_sites:-}"
        for site in "${sites[@]}"; do
            site=$(echo "$site" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')
            [[ -n "$site" ]] && known_sites_cmd+=" --known-sites $site"
        done
    fi

    local padding_opt=""
    [[ -n "${padding:-}" ]] && padding_opt=" --interval-padding $padding"

    echo "Running BQSR on $(basename "$input_bam") for interval $interval_code" | tee -a "$LOG_FILE"
    execute "$GATK_CMD BaseRecalibrator \
      -I \"$input_bam\" \
      -R \"$reference\" \
      -L \"$interval_split\"$padding_opt \
      $known_sites_cmd \
      -O \"$recal_table\" \
      --tmp-dir \"${output_directory%/}/bqsr_tmp\"" \
      "BaseRecalibrator for $(basename "$input_bam") on interval $interval_code"
}

########################################
# Main Pipeline
########################################
main_pipeline() {

    # Resolve GATK launcher once
    if command -v gatk >/dev/null 2>&1; then
        GATK_CMD="gatk"
    elif [[ -n "${gatk_path:-}" ]]; then
        if [[ "$gatk_path" == *.jar ]]; then
            GATK_CMD="java -jar \"$gatk_path\""
        elif [[ -x "$gatk_path" ]]; then
            GATK_CMD="\"$gatk_path\""
        else
            echo "Error: gatk_path set but not executable: $gatk_path" | tee -a "$LOG_FILE"
            exit 1
        fi
    else
        echo "Error: Required tool 'gatk' not found; set PATH or gatk_path." | tee -a "$LOG_FILE"
        exit 1
    fi
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    LOG_FILE="${output_directory}/pipeline_$(date +'%Y%m%d_%H%M%S').log"

    execute "mkdir -p '${output_directory}/logs_mark_dup' '${output_directory}/logs_bqsr' '${output_directory}/qs_tmp' '${output_directory}/md_tmp' '${output_directory}/fix_tmp' '${output_directory}/bqsr_tmp' '${mapping_output_dir}'" "Setup output directories"

    if [[ "$ans" == "yes" ]]; then
        execute "echo '=== BWA Indexing ==='" "BWA Indexing step"
        BWA_Indexing
    fi

    execute "echo '=== Mapping ==='" "Mapping step"
    Mapping

    if [[ "$markduplicates_enabled" == "yes" ]]; then
	export mapping_output_dir
        execute "echo '=== MarkDuplicates ==='" "MarkDuplicates step"
        pushd "$mapping_output_dir" > /dev/null
        execute "ls sorted_*.bam > bam_list.txt" "Listing sorted BAMs"

        export -f QuerySort MarkDuplicates SetNmMdAndUqTags execute
        export mapping_output_dir output_directory threads picard_path reference gatk_path LOG_FILE threads

        parallel --jobs "$threads" QuerySort :::: bam_list.txt
        execute "ls qs_*.bam > qs_list.txt" "Listing query-sorted BAMs"
        parallel --jobs "$threads" MarkDuplicates :::: qs_list.txt
        execute "ls md_*.bam > md_list.txt" "Listing duplicate-marked BAMs"
        parallel --jobs "$threads" SetNmMdAndUqTags :::: md_list.txt

        execute "mv bam_list.txt qs_list.txt md_list.txt '${output_directory}/logs_mark_dup/'" "Archiving BAM lists"
        popd > /dev/null
    fi
    if [[ "$bqsr_enabled" == "yes" ]]; then
        echo "=== BQSR ===" | tee -a "$LOG_FILE"

        # Split intervals to an absolute directory
        mkdir -p "${output_directory%/}/bqsr_tmp"
        SplitIntervals

        # Expand intervals to concrete files (no wildcards later)
        INTERVALS_DIR="${output_directory%/}/intervals"
        mapfile -t INTERVALS < <(printf '%s
' "$INTERVALS_DIR"/*.interval_list)
        if (( ${#INTERVALS[@]} == 0 )); then
            echo "ERROR: No interval_list files in $INTERVALS_DIR" | tee -a "$LOG_FILE"
            exit 1
        fi

        pushd "$mapping_output_dir" > /dev/null
        export -f BQSR CombineTables ApplyBQSR
        export GATK_CMD output_directory threads reference picard_path gatk_path intervals padding known_sites mapping_output_dir LOG_FILE

        # Per-sample × per-interval shards
        for bam_file in f_*.bam; do
            [[ -e "$bam_file" ]] || continue
            echo "Running BQSR for $(basename "$bam_file")" | tee -a "$LOG_FILE"
            parallel --jobs "${threads:-1}" BQSR {1} {2} ::: "$bam_file" ::: "${INTERVALS[@]}"
        done

        # Combine per-sample tables
        for bam_file in f_*.bam; do
            [[ -e "$bam_file" ]] || continue
            base_name=$(basename "$bam_file" .bam)
            tables=()
            for f in "${INTERVALS[@]}"; do
                icode=$(basename "$f" .interval_list)
                tables+=( "${base_name}_${icode}.table" )
            done
            CombineTables "$base_name" "${tables[@]}"
        done

        # Apply across samples in parallel
        parallel --jobs "${threads:-1}" ApplyBQSR {1} ::: f_*.bam
        popd > /dev/null
    fi
	
        parallel --jobs "$threads" \
		 :::: bam_list.txt \
		 ::: "${INTERVALS[@]}"

	while read -r bam; do
	    base=$(basename "$bam" .bam)
	    tables=( "${base}"_*.table )
	    if (( ${#tables[@]} == 0 )); then
		echo "WARNING: No tables for ${base}; skipping GatherTables/ApplyBQSR"
		continue
	    fi
	    GatherTables "$base" "${tables[@]}"
	    ApplyBQSR "$bam"
	done < bam_list.txt
	
        popd > /dev/null


    execute "echo '=== QC ==='" "QC section"
    run_qc

    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    start_s=$(date -d "$start_time" +%s)
    end_s=$(date -d "$end_time" +%s)
    dur=$((end_s - start_s))
    execute "echo 'Total runtime: $((dur/3600))h $(((dur%3600)/60))m $((dur%60))s'" "Recording total runtime"
}

########################################
# Entry Point
########################################
if [[ "$#" -lt 1 ]]; then
    echo "Usage: $0 <config.yaml> [--dry-run|--no-dry-run] [--verbose|--silent]"
    exit 1
fi

parse_yaml "$1"

[[ "${output_directory}" != */ ]] && output_directory+="/"
mapping_output_dir="${output_directory}bwa_files/"
mkdir -p "$mapping_output_dir"

for arg; do
    case "$arg" in
        --no-dry-run) DRY_RUN=false ;;
        --dry-run) DRY_RUN=true ;;
        --silent) VERBOSE_MODE="silent" ;;
        --verbose) VERBOSE_MODE="verbose" ;;
    esac
done

start=$(date +%s)
main_pipeline
