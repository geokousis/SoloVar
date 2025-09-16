#!/bin/bash
set -euo pipefail

########################################
# Simple YAML parser (key: value pairs)
########################################
parse_yaml() {
    local yaml_file="$1"
    while IFS=":" read -r key value; do
        key="$(echo "$key" | xargs)"
        value="$(echo "$value" | xargs)"
        # ignore empty keys or comments
        if [[ -n "$key" && ! "$key" =~ ^# ]]; then
            eval "${key}='${value}'"
        fi
    done < "$yaml_file"
}

########################################
# Execute helper function (with logging)
########################################
execute() {
    local cmd="$1"
    local description="$2"
    local tool
    tool=$(echo "$cmd" | awk '{print $1}')

    if ! command -v "$tool" &>/dev/null; then
        echo "Error: Required tool '$tool' not found. Please install it before proceeding." | tee -a "$LOG_FILE"
        exit 127
    fi
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $description" | tee -a "$LOG_FILE"
    echo "[CMD] $cmd" >> "$LOG_FILE"
    eval "$cmd" >> "$LOG_FILE" 2>&1
    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "Error: '$description' failed (exit code $status). See log: $LOG_FILE" | tee -a "$LOG_FILE"
        exit $status
    fi
}

########################################
# Main
########################################
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 config.yaml"
    exit 1
fi

CONFIG_FILE="$1"
parse_yaml "$CONFIG_FILE"

# Ensure out_dir ends with slash
[[ "${out_dir}" != */ ]] && out_dir="${out_dir}/"

# Variables from YAML
BED="${bed}"
REF="${ref}"
ACC="${acc}"
OUT_DIR="${out_dir}"
INPUT_SAMPLES="${inputsamples}"
CELLULARITY_FILE="${cellularity_file}"
THREADS="${threads}"

# Prepare output directories and log
RESULT_DIR="${OUT_DIR}CNV_results"
mkdir -p "$RESULT_DIR"

LOG_FILE="${OUT_DIR}pipeline_$(date +'%Y%m%d_%H%M%S').log"
echo "Pipeline started at $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"

# Derive BED basename
BEDNAME=$(basename "$BED" .bed)
echo "[INFO] Using BED: $BED (basename: $BEDNAME)" | tee -a "$LOG_FILE"

echo "[STEP] Autobinning all samples" | tee -a "$LOG_FILE"
pushd "$RESULT_DIR" >/dev/null
SAMPLES=$(tr '\n' ' ' < "$INPUT_SAMPLES")
execute "cnvkit.py autobin $SAMPLES -t '$BED' -g '$ACC' --annotate '$annotation_file' --short-names > '${BEDNAME}.autobin.out'" \
        "Autobinning samples"
execute "cnvkit.py reference -o '${RESULT_DIR}/FlatReference.cnn' -f '$REF' -t '${BEDNAME}.target.bed' -a '${BEDNAME}.antitarget.bed'" \
        "Creating FlatReference"
popd >/dev/null

# Export for parallel
export BED ACC RESULT_DIR REF BEDNAME CELLULARITY_FILE INPUT_SAMPLES LOG_FILE
export -f execute

########################################
# Process each sample (no annotation/gene steps)
########################################
process_sample() {
    local bam="$1"
    local base=$(basename "$bam" .bam)
    echo "[INFO] Processing sample: $base" | tee -a "$LOG_FILE"

    # Coverage
    execute "cnvkit.py coverage '$bam' '${RESULT_DIR}/${BEDNAME}.target.bed' -o '${RESULT_DIR}/${base}.targetcoverage.cnn'" \
            "Coverage (target) for $base"
    execute "cnvkit.py coverage '$bam' '${RESULT_DIR}/${BEDNAME}.antitarget.bed' -o '${RESULT_DIR}/${base}.antitargetcoverage.cnn'" \
            "Coverage (antitarget) for $base"

    # Fix and segment
    execute "cnvkit.py fix '${RESULT_DIR}/${base}.targetcoverage.cnn' '${RESULT_DIR}/${base}.antitargetcoverage.cnn' '${RESULT_DIR}/FlatReference.cnn' -o '${RESULT_DIR}/${base}.cnr'" \
            "Fix for $base"
    execute "cnvkit.py segment '${RESULT_DIR}/${base}.cnr' --drop-low-coverage -m haar -o '${RESULT_DIR}/${base}.cns'" \
            "Segment for $base"
    execute "cnvkit.py segmetrics -s ${RESULT_DIR}/${base}.cn{s,r} --ci -o '${RESULT_DIR}/${base}.segmetrics.cns'" \
            "Segmetrics for $base"

    # Determine purity
    local sample_id
    sample_id=$(echo "$base" | grep -o 'MB[0-9]\+')
    local cell
    cell=$(grep -P "^${sample_id}\s" "$CELLULARITY_FILE" | awk '{print $2}' | tr -d '%') || cell="ND"
    if [[ "$cell" != "ND" ]]; then
        local purity
        purity=$(printf "%.2f" "$(echo "$cell/100" | bc -l)")
        echo "[INFO] Using purity $purity for $base" | tee -a "$LOG_FILE"
        execute "cnvkit.py call '${RESULT_DIR}/${base}.segmetrics.cns' --filter ci --drop-low-coverage -m clonal --gender x --purity $purity -o '${RESULT_DIR}/${base}.call.cns'" \
                "CNV call for $base with purity"
    else
        echo "[WARN] No purity info for $base; using threshold method" | tee -a "$LOG_FILE"
        execute "cnvkit.py call '${RESULT_DIR}/${base}.segmetrics.cns' --filter ci --drop-low-coverage -m threshold --gender x -o '${RESULT_DIR}/${base}.call.cns'" \
                "CNV call for $base (threshold)"
    fi
}
export -f process_sample

########################################
# Run in parallel
########################################
parallel -j "$THREADS" process_sample ::: $(< "$INPUT_SAMPLES")

echo "Pipeline completed at $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG_FILE"
