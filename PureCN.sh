#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

###############################################################################
# CONFIG (your header; unchanged + VCF_PoN added)
###############################################################################
BED="/home/kousis/work/meth/ILCEM/S07604514_hg38_1/S07604514_Regions.clean.bed"
PURECN="/home/kousis/R/x86_64-pc-linux-gnu-library/4.5/PureCN/extdata"
REF="/home/kousis/work/meth/Dra-GATK/mutec/reference/Homo_sapiens_assembly38.fasta"
BW="/home/kousis/work/meth/ILCEM/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
CORES=10

# Intervals / ref outputs
OUT_REF="/home/kousis/work/meth/ILCEM/PureCN/ref"
OUT_INT="$OUT_REF/baits_hg38_intervals.txt"
EXPORT="$OUT_REF/baits_optimized_hg38.bed"

# Coverage inputs/outputs
BAM_LIST="/media/storage/kousis/meth/ILCEM/bam_list.list"           # tumors (.list)
NBAM_LIST="/media/storage/kousis/meth/ILCEM/normal_bams.list"       # normals (.list)
OUT_COV="/home/kousis/work/meth/ILCEM/PureCN/COV"

# Results root used in your examples for PureCN outputs
OUT="/home/kousis/work/meth/ILCEM/PureCN"

# VCF list you’ll provide (as shown)
VCF_LIST="vcf.list"

# hg38 specifics
GENOME="hg38"
SNP_BLACKLIST="/home/kousis/work/meth/ILCEM/hg38_simpleRepeats.bed"

# VCF Panel of Normals for NormalDB step
# VCF_PoN="/home/kousis/work/meth/ILCEM/vcf_normal/vcf/Merged.vcf.gz"   # <-- set this to your PoN VCF (bgzipped + indexed)

# Extra PureCN tuning
MIN_BASE_QUAL=20
INTERVAL_PADDING=50

###############################################################################
# STEP SELECTION (run e.g. --steps 5,6)  defaults to "all"
###############################################################################
STEPS="${STEPS:-all}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps) STEPS="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
run_step() { local s="$1"; [[ "$STEPS" == "all" || ",$STEPS," == *",$s,"* ]]; }

###############################################################################
# LOGGING
###############################################################################
if [[ -t 1 ]]; then
  BOLD="$(tput bold)"; RED="$(tput setaf 1)"; GREEN="$(tput setaf 2)"
  YELLOW="$(tput setaf 3)"; BLUE="$(tput setaf 4)"; RESET="$(tput sgr0)"
else
  BOLD=""; RED=""; GREEN=""; YELLOW=""; BLUE=""; RESET=""
fi
LOG_DIR="${OUT}/logs"; mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/purecn_$(date +%Y%m%d_%H%M%S).log"
ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo -e "[$(ts)] ${BLUE}${BOLD}$*${RESET}" | tee -a "$LOG_FILE"; }
ok()  { echo -e "[$(ts)] ${GREEN}$*${RESET}" | tee -a "$LOG_FILE"; }
warn(){ echo -e "[$(ts)] ${YELLOW}WARNING:${RESET} $*" | tee -a "$LOG_FILE"; }
err() { echo -e "[$(ts)] ${RED}ERROR:${RESET} $*" | tee -a "$LOG_FILE" >&2; }
trap 'err "Pipeline failed at line $LINENO. See log: $LOG_FILE"' ERR

require_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing command: $1"; exit 127; }; }
check_file()  { [[ -s "$1" ]] || { err "Missing/empty: $1"; exit 1; }; }
ensure_dir()  { mkdir -p "$1"; }
count_list_lines(){ grep -v '^\s*$' "$1" | grep -v '^\s*#' | wc -l | tr -d ' '; }

require_cmd Rscript
require_cmd samtools
command -v figlet >/dev/null 2>&1 && figlet "PureCN" | tee -a "$LOG_FILE"
log "Log file: $LOG_FILE"

ensure_dir "$OUT_REF" "$OUT_COV/tumors" "$OUT_COV/normal" "$OUT"

###############################################################################
# [1] INTERVAL CREATION (hg38)
###############################################################################
if run_step 1; then
  log "--------------- [1] Interval Creation ---------------"
  check_file "$PURECN/IntervalFile.R"
  check_file "$BED"; check_file "$REF"; check_file "$BW"
  Rscript "$PURECN/IntervalFile.R" \
    --in-file "$BED" \
    --fasta "$REF" \
    --genome hg38 \
    --export "$EXPORT" \
    --mappability "$BW" \
    --out-file "$OUT_INT" \
    --force 2>&1 | tee -a "$LOG_FILE"
  ok "Intervals created → $OUT_INT and $EXPORT"
fi

###############################################################################
# [2] TUMOR COVERAGE
###############################################################################
if run_step 2; then
  log "--------------- [2] Tumor Coverage ---------------"
  check_file "$PURECN/Coverage.R"; check_file "$BAM_LIST"; check_file "$OUT_INT"
  Rscript "$PURECN/Coverage.R" \
    --out-dir "$OUT_COV/tumors" \
    --bam "$BAM_LIST" \
    --intervals "$OUT_INT" \
    --cores "$CORES" 2>&1 | tee -a "$LOG_FILE"
  ok "Tumor coverage complete → $OUT_COV/tumors"
fi

###############################################################################
# [3] NORMAL INDEXING + COVERAGE  (NormalDB moved to step 4)
###############################################################################
if run_step 3; then
  log "--------------- [3] Normal indexing + coverage ---------------"
  check_file "$NBAM_LIST"

  # index normals if needed
  while read -r bam; do
    [[ -z "$bam" || "$bam" =~ ^# ]] && continue
    if [[ -f "$bam" ]]; then
      if [[ ! -f "${bam}.bai" ]]; then
        log "Indexing $bam"; samtools index "$bam"
      fi
    else
      warn "Normal BAM not found: $bam"
    fi
  done < "$NBAM_LIST"

  # normal coverage
  check_file "$PURECN/Coverage.R"
  Rscript "$PURECN/Coverage.R" \
    --out-dir "$OUT_COV/normal" \
    --bam "$NBAM_LIST" \
    --intervals "$OUT_INT" \
    --cores "$CORES" 2>&1 | tee -a "$LOG_FILE"
  ok "Normal coverage complete → $OUT_COV/normal"
fi

###############################################################################
# [4] NormalDB build (uses normals' *_loess.txt.gz + your VCF_PoN)
###############################################################################
if run_step 4; then
  log "--------------- [4] NormalDB build (hg38) ---------------"
  check_file "$PURECN/NormalDB.R"
  check_file "$SNP_BLACKLIST" || true  # not used here, but often co-located
  # inputs
  NORMAL_COV_LIST="$OUT_REF/example_normal_coverages.list"
  ls -1 "$OUT_COV/normal"/*_loess.txt.gz > "$NORMAL_COV_LIST"
  check_file "$NORMAL_COV_LIST"
  # PoN VCF (bgzipped & indexed)
  # check_file "$VCF_PoN"
  # [[ -s "${VCF_PoN}.tbi" || -s "${VCF_PoN}.csi" ]] || { err "PoN VCF missing index: ${VCF_PoN}.tbi|.csi"; exit 1; }

  Rscript "$PURECN/NormalDB.R" \
    --out-dir "$OUT_REF" \
    --coverage-files "$NORMAL_COV_LIST" \
    --genome hg38 \
    --assay agilent_v6 \
    --force \
    2>&1 | tee -a "$LOG_FILE"

  check_file "$OUT_REF/normalDB_agilent_v6_hg38.rds"
  ok "NormalDB created → $OUT_REF/normalDB_agilent_v6_hg38.rds"
fi

###############################################################################
# HELPERS (for steps 5–6)
###############################################################################
# Derive SAMPLEID from VCF filename to match your tumor coverage prefix:
#   filtered_pass_filtered_MB10_recalibrated.g.vcf  ->  f_MB10_recalibrated
derive_sampleid_from_vcf() {
  local v; v="$(basename "$1")"
  local core=""
  if [[ "$v" =~ (MB[0-9]+_recalibrated) ]]; then
    core="${BASH_REMATCH[1]}"
  else
    core="${v%.g.vcf}"
    core="${core#filtered_pass_filtered_}"
  fi
  printf "f_%s" "$core"
}

# Find the tumor BAM path corresponding to SAMPLEID by fuzzy match against BAM_LIST
# (looks for the MB##_recalibrated core in the BAM basename)
find_tumor_bam_for_sample() {
  local sampleid="$1"
  local suf="${sampleid#f_}"    # MB10_recalibrated
  local bam
  while read -r bam; do
    [[ -z "$bam" || "$bam" =~ ^# ]] && continue
    if [[ "$(basename "$bam")" == *"$suf"* ]]; then
      echo "$bam"; return 0
    fi
  done < "$BAM_LIST"
  return 1
}

# Compute minDepth as 20% of mean on-target coverage using samtools depth
# Requires: BAM indexed, and $EXPORT (BED of targets)
compute_min_depth() {
  local bam="$1"
  local bed="$2"   # typically $EXPORT
  local mean
  mean="$(samtools depth -b "$bed" "$bam" \
          | awk '{sum+=$3; n++} END{if(n>0) printf "%.2f", sum/n; else print 0}')"
  awk -v m="$mean" 'BEGIN{v=int(m*0.2+0.5); if(v<1)v=1; print v}'
}

###############################################################################
# [5] PURECN PER-SAMPLE (iterate VCFs; no mapping-bias-file; hg38)
###############################################################################
if run_step 5; then
  log "--------------- [5] PureCN per-sample (hg38) ---------------"
  check_file "$PURECN/PureCN.R"
  check_file "$VCF_LIST"
  check_file "$OUT_REF/normalDB_agilent_v6_hg38.rds"
  check_file "$OUT_INT"     # hg38 intervals created in step 1

  while read -r VCF; do
    [[ -z "$VCF" || "$VCF" =~ ^# ]] && continue
    check_file "$VCF"

    SAMPLEID="$(derive_sampleid_from_vcf "$VCF")"
    OUT_SAMPLE="$OUT/$SAMPLEID"
    ensure_dir "$OUT_SAMPLE"

    # prefer exact match first
    TUMOR_LOESS="$OUT_COV/tumors/${SAMPLEID}_coverage_loess.txt.gz"
    if [[ ! -f "$TUMOR_LOESS" ]]; then
      # fallback: match by suffix (minus leading "f_"), choose first sorted
      suf="${SAMPLEID#f_}"
      alt="$(ls -1 "$OUT_COV/tumors"/*"${suf}"*_coverage_loess.txt.gz 2>/dev/null | sort | head -n1 || true)"
      [[ -n "$alt" ]] && TUMOR_LOESS="$alt"
    fi
    check_file "$TUMOR_LOESS"

    log "Sample: $SAMPLEID"
    log "  VCF:   $VCF"
    log "  LOESS: $TUMOR_LOESS"

    Rscript "$PURECN/PureCN.R" \
      --out "$OUT_SAMPLE/$SAMPLEID" \
      --tumor "$TUMOR_LOESS" \
      --sampleid "$SAMPLEID" \
      --vcf "$VCF" \
      --fun-segmentation PSCBS \
      --normaldb "$OUT_REF/normalDB_agilent_v6_hg38.rds" \
      --intervals "$OUT_INT" \
      --snp-blacklist "$SNP_BLACKLIST" \
      --genome "$GENOME" \
      --model betabin \
      --min-base-quality "$MIN_BASE_QUAL" \
      --interval-padding "$INTERVAL_PADDING" \
      --force --post-optimize --seed 123 \
      --cores "$CORES" \
      2>&1 | tee -a "$LOG_FILE"

    # sanity check: RDS should exist
    if [[ ! -s "$OUT_SAMPLE/${SAMPLEID}.rds" ]]; then
      err "Missing RDS after PureCN: $OUT_SAMPLE/${SAMPLEID}.rds"
      ls -l "$OUT_SAMPLE" | tee -a "$LOG_FILE"
      exit 1
    fi

    ok "PureCN done → $OUT_SAMPLE"
  done < "$VCF_LIST"
fi

###############################################################################
# [6] CALLABLE (bedtools) → FILTER → Dx (hg38)  [sorted stream + sorted targets]
###############################################################################
if run_step 6; then
  log "--------------- [6] Callable (bedtools) + Dx (hg38) ---------------"
  require_cmd bedtools
  check_file "$PURECN/Dx.R"
  check_file "$PURECN/FilterCallableLoci.R"
  check_file "$REF"
  check_file "$EXPORT"
  check_file "$BAM_LIST"

  # Ensure FASTA index exists for consistent chromosome ordering
  if [[ ! -s "${REF}.fai" ]]; then
    log "Indexing reference FASTA for sorting: ${REF}.fai"
    samtools faidx "$REF"
  fi

  # Create a sorted targets BED (once)
  TARGETS_SORTED="$OUT_REF/baits_optimized_hg38.sorted.bed"
  if [[ ! -s "$TARGETS_SORTED" ]]; then
    log "Sorting targets BED → $TARGETS_SORTED"
    bedtools sort -faidx "${REF}.fai" -i "$EXPORT" > "$TARGETS_SORTED"
  fi

  while read -r VCF; do
    [[ -z "$VCF" || "$VCF" =~ ^# ]] && continue
    SAMPLEID="$(derive_sampleid_from_vcf "$VCF")"
    OUT_SAMPLE="$OUT/$SAMPLEID"
    check_file "$OUT_SAMPLE/${SAMPLEID}.rds"   # produced by step 5

    # Find tumor BAM and ensure it’s indexed
    BAM="$(find_tumor_bam_for_sample "$SAMPLEID" || true)"
    if [[ -z "$BAM" || ! -f "$BAM" ]]; then
      err "Could not match tumor BAM for $SAMPLEID (looked in $BAM_LIST)"; exit 1
    fi
    [[ -f "${BAM}.bai" ]] || { log "Indexing tumor BAM for $SAMPLEID"; samtools index "$BAM"; }

    pushd "$OUT_SAMPLE" >/dev/null

    CALLABLE_RAW="${SAMPLEID}_callable_status.bed"
    if [[ ! -s "$CALLABLE_RAW" ]]; then
      # Compute minDepth = 20% of mean on-target coverage (from BAM over targets)
      MINDEPTH="$(compute_min_depth "$BAM" "$EXPORT")"
      log "Callable (bedtools): $SAMPLEID  BAM=$(basename "$BAM")  minDepth=${MINDEPTH}"

      # Coverage → filter by depth → intersect with sorted targets → sort by faidx → merge
      bedtools genomecov -ibam "$BAM" -bga -split \
        | awk -v md="$MINDEPTH" '$4>=md {print $1"\t"$2"\t"$3"\tCALLABLE"}' \
        | bedtools intersect -a - -b "$TARGETS_SORTED" -sorted \
        | bedtools sort -faidx "${REF}.fai" -i - \
        | bedtools merge -i - -c 4 -o distinct \
        > "$CALLABLE_RAW"
    else
      log "Callable BED exists for $SAMPLEID → $CALLABLE_RAW"
    fi

    # Validate callable bed
    if [[ ! -s "$CALLABLE_RAW" ]]; then
      err "Callable BED is empty for $SAMPLEID: $CALLABLE_RAW"
      ls -l | tee -a "$LOG_FILE"
      exit 1
    fi

    # Filter/normalize to ensure 4 columns with CALLABLE tag
    CALLABLE_FILT="${SAMPLEID}_callable_status_filtered.bed"
    awk 'NF>=3{if(NF==3){print $1"\t"$2"\t"$3"\tCALLABLE"} else if($4=="CALLABLE"){print $1"\t"$2"\t"$3"\t"$4}}' \
      "$CALLABLE_RAW" \
      | bedtools sort -faidx "${REF}.fai" -i - \
      > "$CALLABLE_FILT"

    if [[ ! -s "$CALLABLE_FILT" ]]; then
      err "Filtered callable BED is empty for $SAMPLEID: $CALLABLE_FILT"
      head -n 5 "$CALLABLE_RAW" | tee -a "$LOG_FILE" || true
      exit 1
    fi

    # Dx on all callable (exclude simple repeats; search COSMIC signatures)
    check_file "$SNP_BLACKLIST"
    Rscript "$PURECN/Dx.R" \
      --out "$OUT_SAMPLE/$SAMPLEID" \
      --rds "$OUT_SAMPLE/${SAMPLEID}.rds" \
      --callable "$CALLABLE_FILT" \
      --exclude "$SNP_BLACKLIST" \
      --force \
      --signatures \
      2>&1 | tee -a "$LOG_FILE"

    # Restrict callable to CDS (and drop HLA)
    CALLABLE_CDS="${SAMPLEID}_callable_status_filtered_cds.bed"
    Rscript "$PURECN/FilterCallableLoci.R" \
      --genome "$GENOME" \
      --in-file "$CALLABLE_FILT" \
      --force \      
      --out-file "$CALLABLE_CDS" \
      --exclude '^HLA' \
      2>&1 | tee -a "$LOG_FILE"

    if [[ ! -s "$CALLABLE_CDS" ]]; then
      err "CDS callable BED is empty for $SAMPLEID: $CALLABLE_CDS"
      exit 1
    fi

    # Dx on CDS-only callable
    Rscript "$PURECN/Dx.R" \
      --out "$OUT_SAMPLE/${SAMPLEID}_cds" \
      --rds "$OUT_SAMPLE/${SAMPLEID}.rds" \
      --callable "$CALLABLE_CDS" \
      --force \
      --exclude "$SNP_BLACKLIST" \
      2>&1 | tee -a "$LOG_FILE"

    ok "Callable + Dx finished → $OUT_SAMPLE"
    popd >/dev/null
  done < "$VCF_LIST"
fi

ok "Pipeline finished!   Steps executed: ${STEPS}"
