#!/usr/bin/env python3
"""
Inputs:
  --cnvkit-dir : directory containing *.call.cns (tab-delim; columns incl. chrom/start/end/cn)
  --arms       : BED of chromosome arms (0-based, half-open). Columns: chrom start end name|arm
  --targets    : (optional) BED of capture targets (0-based, half-open)

Outputs (same structure as before):
  outdir/
    segments_classified/
      <sample>.arm_summary.tsv     # per-sample arm calls
    arm_calls_per_sample.tsv       # stack of all per-sample calls
    cohort_arm_summary.tsv         # counts of non-neutral events across cohort (+ sample lists)
    thresholds_used.json           # thresholds/config used

Usage:
  python cnvkit_arm_aggregator_final.py \
      --cnvkit-dir /path/to/call_cns \
      --arms hg38_chrom_arms.bed \
      --outdir cnv_summary_results \
      [--targets wes_targets.bed]
"""

import argparse
import os
import pandas as pd
from collections import defaultdict
import json
import math
from bisect import bisect_left

# ---------------------------
# Configurable thresholds
# ---------------------------
THRESHOLDS = {
    "CN_HOMO_DEL": 0.2,           # <=0.2: homozygous deletion
    "CN_HEMI_DEL": 1.5,           # <2:   hemizygous deletion
    "CN_LOW_AMP": 6,              # >=6:  low-level amplification
    "CN_HIGH_AMP": 10,            # >=10: high-level amplification
    "ARM_FRACTION_THRESHOLD": 0.7,
    "PARTIAL_ARM_MIN_FRAC": 0.3,
    "MIN_SEG_LEN_BP": 1000
}

# ---------------------------
# Helpers
# ---------------------------
def intersect_len(s1, e1, s2, e2) -> int:
    """Length of intersection between [s1,e1) and [s2,e2) in half-open coords."""
    return max(0, min(e1, e2) - max(s1, s2))

def classify_segment(cn: float) -> str:
    """Classify a segment by absolute CN."""
    if cn is None or (isinstance(cn, float) and math.isnan(cn)):
        return "Neutral"  # conservative default
    if cn <= THRESHOLDS["CN_HOMO_DEL"]:
        return "Homozygous Deletion"
    elif cn < THRESHOLDS["CN_HEMI_DEL"]:
        return "Hemizygous Deletion"
    elif cn >= THRESHOLDS["CN_HIGH_AMP"]:
        return "High Amplification"
    elif cn >= THRESHOLDS["CN_LOW_AMP"]:
        return "Low-level Amplification"
    else:
        return "Neutral"

def read_cns(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = [c.strip().lower() for c in df.columns]
    rename_map = {"chromosome": "chrom", "chr": "chrom"}
    df.rename(columns=rename_map, inplace=True)
    required = {"chrom", "start", "end", "cn"}
    if not required.issubset(df.columns):
        raise ValueError(f"{path} missing required columns: {required} (found: {list(df.columns)})")
    return df[["chrom", "start", "end", "cn"]]

def read_bed(path: str, names=("chrom", "start", "end", "name")) -> pd.DataFrame:
    """Read BED file (0-based start, half-open end)."""
    with open(path) as f:
        first = f.readline().strip().split()
    ncols = len(first)
    cols = names[:ncols]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    if len(df.columns) == 3:
        df.columns = ["chrom", "start", "end"]
    return df

def has_chr_prefix(series: pd.Series) -> bool:
    return series.astype(str).str.startswith("chr").any()

def normalize_chrom_style(df: pd.DataFrame, use_chr_prefix: bool) -> pd.DataFrame:
    """Ensure chrom column style matches arms BED. Minimal normalization for 'chr' prefix."""
    def add_chr(c):  # e.g., '1' -> 'chr1'
        c = str(c)
        return c if c.startswith("chr") else "chr" + c
    def strip_chr(c):  # e.g., 'chr1' -> '1'
        c = str(c)
        return c[3:] if c.startswith("chr") else c
    if use_chr_prefix:
        df["chrom"] = df["chrom"].astype(str).map(add_chr)
    else:
        df["chrom"] = df["chrom"].astype(str).map(strip_chr)
    return df

# ----- interval utilities (no external deps) -----
def merge_intervals(intervals):
    """Merge overlapping or touching intervals. 'intervals' is iterable of (s,e), half-open."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:  # overlap or touching
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def clip_intervals_to_window(intervals, win_start, win_end):
    """Intersect a list of (s,e) with [win_start,win_end)."""
    out = []
    # binary search to first interval with end > win_start
    ends = [e for _, e in intervals]
    i = bisect_left(ends, win_start)
    # iterate forward until start >= win_end
    n = len(intervals)
    while i < n and intervals[i][0] < win_end:
        s, e = intervals[i]
        isect = intersect_len(s, e, win_start, win_end)
        if isect > 0:
            out.append((max(s, win_start), min(e, win_end)))
        i += 1
    return out

def overlap_len_with_intervals(intervals, s, e) -> int:
    """Sum of overlaps between [s,e) and a list of non-overlapping (merged) intervals."""
    if not intervals:
        return 0
    ends = [ie for _, ie in intervals]
    i = bisect_left(ends, s)
    total = 0
    n = len(intervals)
    while i < n and intervals[i][0] < e:
        isect = intersect_len(intervals[i][0], intervals[i][1], s, e)
        if isect > 0:
            total += isect
        i += 1
    return total

def get_arm_name(row) -> str:
    if "name" in row.index:
        return str(row["name"])
    if "arm" in row.index:
        return str(row["arm"])
    return f"{row['chrom']}:{row['start']}-{row['end']}"

# ---------------------------
# Core computation
# ---------------------------
def build_arm_and_target_index(arms_df, targets_df=None):
    """
    Pre-index arms by chromosome and build merged target intervals:
      - targets_by_chr[chrom] = merged non-overlapping intervals
      - targets_by_arm[(chrom, arm_name)] = merged intervals clipped to that arm
    """
    arms_by_chr = defaultdict(list)  # chrom -> list of dicts with start,end,name
    for _, row in arms_df.iterrows():
        arms_by_chr[row["chrom"]].append({
            "start": int(row["start"]),
            "end": int(row["end"]),
            "name": get_arm_name(row)
        })
    # sort arms per chrom
    for chrom in list(arms_by_chr.keys()):
        arms_by_chr[chrom].sort(key=lambda r: (r["start"], r["end"]))

    targets_by_chr = {}
    targets_by_arm = {}

    if targets_df is not None:
        # merge per chromosome
        tmp = defaultdict(list)
        for _, r in targets_df.iterrows():
            tmp[r["chrom"]].append((int(r["start"]), int(r["end"])))
        for chrom, lst in tmp.items():
            targets_by_chr[chrom] = merge_intervals(lst)

        # pre-clip targets into each arm (still merged)
        for chrom, arms in arms_by_chr.items():
            merged_targets = targets_by_chr.get(chrom, [])
            for arm in arms:
                clipped = clip_intervals_to_window(merged_targets, arm["start"], arm["end"])
                targets_by_arm[(chrom, arm["name"])] = clipped

    return arms_by_chr, targets_by_arm  # targets_by_arm empty dict if no targets

def compute_arm_summary(cns_df, arms_by_chr, sample_name, targets_by_arm=None):
    """
    For each arm, aggregate coverage by event class using:
      - base overlap: seg ∩ arm (if no targets)
      - weighted overlap: seg ∩ arm ∩ merged(targets) (if targets provided)
    Then call the arm using non-neutral event fractions only.
    """
    arm_cov = defaultdict(lambda: defaultdict(float))
    targets_given = bool(targets_by_arm)

    # iterate segments
    for _, seg in cns_df.iterrows():
        s_chr = seg["chrom"]
        s_start = int(seg["start"])
        s_end = int(seg["end"])
        if s_end - s_start < THRESHOLDS["MIN_SEG_LEN_BP"]:
            continue

        seg_class = classify_segment(seg["cn"])
        if s_chr not in arms_by_chr:
            continue

        for arm in arms_by_chr[s_chr]:
            # quick reject
            base_overlap = intersect_len(s_start, s_end, arm["start"], arm["end"])
            if base_overlap == 0:
                continue

            if targets_given:
                arm_targets = targets_by_arm.get((s_chr, arm["name"]), [])
                if not arm_targets:
                    continue
                # seg ∩ arm ∩ targets (targets already clipped to arm)
                eff = overlap_len_with_intervals(arm_targets, s_start, s_end)
                if eff == 0:
                    continue
                eff_overlap = eff
            else:
                eff_overlap = base_overlap

            arm_cov[arm["name"]][seg_class] += eff_overlap

    # summarize per arm
    rows = []
    for arm_name, cov in arm_cov.items():
        total = float(sum(cov.values()))
        if total == 0:
            continue
        frac = {k: v / total for k, v in cov.items()}

        non_neutral = {k: v for k, v in frac.items() if k != "Neutral"}
        if non_neutral:
            top_event = max(non_neutral, key=non_neutral.get)
            top_frac = non_neutral[top_event]
            if top_frac >= THRESHOLDS["ARM_FRACTION_THRESHOLD"]:
                call = f"Arm-level {top_event}"
            elif top_frac >= THRESHOLDS["PARTIAL_ARM_MIN_FRAC"]:
                call = f"Partial Arm {top_event}"
            else:
                call = "Neutral"
        else:
            top_event = "Neutral"
            top_frac = frac.get("Neutral", 0.0)
            call = "Neutral"

        rows.append({
            "sample": sample_name,
            "arm": arm_name,
            "event": call,
            "fraction": round(float(top_frac if call != "Neutral" else frac.get("Neutral", 0.0)), 3)
        })

    return pd.DataFrame(rows)

# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser(description="Aggregate CNVkit call.cns files to arm-level CNVs using CN values (WES-ready).")
    ap.add_argument("--cnvkit-dir", required=True, help="Directory with *.call.cns files")
    ap.add_argument("--arms", required=True, help="Chromosome arm BED (e.g., hg38_chrom_arms.bed)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--targets", required=False, help="Optional BED file with WES capture targets")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seg_dir = os.path.join(args.outdir, "segments_classified")
    os.makedirs(seg_dir, exist_ok=True)

    # Read arms and (optional) targets
    arms_df = read_bed(args.arms, names=["chrom", "start", "end", "name"])
    targets_df = read_bed(args.targets, names=["chrom", "start", "end"]) if args.targets else None

    # Normalize chrom style across inputs based on arms BED style
    arms_use_chr = has_chr_prefix(arms_df["chrom"])
    arms_df = normalize_chrom_style(arms_df, arms_use_chr)
    if targets_df is not None:
        targets_df = normalize_chrom_style(targets_df, arms_use_chr)

    # Build indices (pre-clip merged targets into each arm for performance)
    arms_by_chr, targets_by_arm = build_arm_and_target_index(arms_df, targets_df)

    all_arm_calls = []

    for fn in sorted(os.listdir(args.cnvkit_dir)):
        if not fn.endswith(".call.cns"):
            continue
        path = os.path.join(args.cnvkit_dir, fn)
        sample_name = os.path.basename(fn).replace(".call.cns", "")
        print(f"[INFO] Processing {sample_name}")

        cns = read_cns(path)
        cns = normalize_chrom_style(cns, arms_use_chr)

        df_arms = compute_arm_summary(cns, arms_by_chr, sample_name, targets_by_arm if targets_df is not None else None)
        if df_arms.empty:
            print(f"[WARN] No arm data for {sample_name} (after filtering).")
            continue

        df_arms.to_csv(os.path.join(seg_dir, f"{sample_name}.arm_summary.tsv"), sep="\t", index=False)
        all_arm_calls.append(df_arms)

    if not all_arm_calls:
        print("No CNV data processed.")
        return

    all_arms = pd.concat(all_arm_calls, ignore_index=True)
    all_arms.to_csv(os.path.join(args.outdir, "arm_calls_per_sample.tsv"), sep="\t", index=False)

    # ---- Cohort summary (non-neutral only) ----
    non_neutral = all_arms[~all_arms["event"].str.contains("Neutral", case=False)]
    if non_neutral.empty:
        print("No non-neutral events found across cohort.")
    else:
        summary = (
            non_neutral.groupby(["arm", "event"])
            .agg(
                count=("sample", "nunique"),
                samples=("sample", lambda x: ",".join(sorted(set(x))))
            )
            .reset_index()
            .sort_values("count", ascending=False)
        )
        summary.to_csv(
            os.path.join(args.outdir, "cohort_arm_summary.tsv"),
            sep="\t",
            index=False
        )
        print(f"[INFO] Saved cohort_arm_summary.tsv with {len(summary)} non-neutral events.")

    # Save config for reproducibility
    with open(os.path.join(args.outdir, "thresholds_used.json"), "w") as f:
        json.dump(THRESHOLDS, f, indent=2)

    print(f"[DONE] Results in: {args.outdir}")

if __name__ == "__main__":
    main()
