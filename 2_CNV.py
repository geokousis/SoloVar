#!/usr/bin/env python3
"""
Summarize the most affected chromosome arms from cohort_arm_summary.tsv
produced by cnvkit_arm_aggregator_final.py

Example:
  python summarize_cohort_arms.py cohort_arm_summary.tsv arm_summary_report.tsv
"""

import sys
import pandas as pd
from collections import defaultdict

if len(sys.argv) < 2:
    print("Usage: python summarize_cohort_arms.py cohort_arm_summary.tsv [output.tsv]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else "arm_event_summary.tsv"

# ---- Load data ----
df = pd.read_csv(input_file, sep="\t")
required_cols = {"arm", "event", "count", "samples"}
if not required_cols.issubset(df.columns):
    raise ValueError(f"Input file must have columns: {required_cols}")

# ---- Parse samples ----
df["samples"] = df["samples"].fillna("").apply(lambda x: [s.strip() for s in str(x).split(",") if s.strip()])

# ---- Summarize per arm ----
arm_summary = []
for arm, subdf in df.groupby("arm"):
    # unique sample set across all events
    all_samples = sorted({s for lst in subdf["samples"] for s in lst})
    total = len(all_samples)

    # event-specific details
    details = []
    for _, row in subdf.iterrows():
        ev = row["event"]
        smps = ",".join(sorted(row["samples"]))
        details.append(f"{ev} ({smps})")

    arm_summary.append({
        "arm": arm,
        "total_samples": total,
        "events": "; ".join(subdf["event"].unique()),
        "details": " | ".join(details),
        "samples": ",".join(all_samples)
    })

# ---- Sort by most affected ----
arm_summary = sorted(arm_summary, key=lambda x: x["total_samples"], reverse=True)
out_df = pd.DataFrame(arm_summary)

# ---- Save ----
out_df.to_csv(output_file, sep="\t", index=False)

# ---- Display top 10 ----
print("\nTop affected chromosome arms:")
print(out_df.head(15).to_string(index=False))
print(f"\n[DONE] Summary written to: {output_file}")
