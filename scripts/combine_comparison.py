#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Count Mucinous vs Lobular matches (ignore ILCEM) from combined SNV+CNV file.")
    ap.add_argument("--combined", required=True, help="Combined SNV+CNV TSV (e.g., SNV_CNV_comp.list)")
    ap.add_argument("--out-per-sample", help="Write per-sample counts TSV")
    ap.add_argument("--out-overall", help="Write overall counts TSV")
    args = ap.parse_args()

    # Read combined table
    df = pd.read_csv(args.combined, sep="\t", dtype=str)

    # Normalize columns we'll use
    for col in ["GeneGroup", "Sample", "Status"]:
        if col not in df.columns:
            raise SystemExit(f"[ERROR] Missing required column '{col}' in {args.combined}")
    df["GeneGroup_norm"] = df["GeneGroup"].astype(str).str.strip().str.lower()
    df["Status"] = df["Status"].astype(str)

    # Keep only mucinous & lobular groups; ignore ilcem and others
    df = df[df["GeneGroup_norm"].isin({"mucinous", "lobular"})].copy()

    # Keep only rows with an actual event
    df = df[df["Status"].isin({"SNV only", "CNV only", "Both"})].copy()

    if df.empty:
        print("[INFO] No Mucinous/Lobular matches found.")
        # still emit empty outputs if requested
        if args.out_overall:
            pd.DataFrame(columns=["Group", "Count"]).to_csv(args.out_overall, sep="\t", index=False)
        if args.out_per_sample:
            pd.DataFrame(columns=["Sample", "Lobular", "Mucinous", "Total"]).to_csv(args.out_per_sample, sep="\t", index=False)
        return

    # Overall counts
    overall = (
        df.groupby("GeneGroup_norm")
          .size()
          .rename_axis("Group")
          .reset_index(name="Count")
    )
    # make title-case group names in output
    overall["Group"] = overall["Group"].str.title()

    # Per-sample counts
    per_sample = (
        df.groupby(["Sample", "GeneGroup_norm"])
          .size()
          .unstack(fill_value=0)
          .rename(columns={"lobular": "Lobular", "mucinous": "Mucinous"})
          .reset_index()
    )
    if "Lobular" not in per_sample.columns: per_sample["Lobular"] = 0
    if "Mucinous" not in per_sample.columns: per_sample["Mucinous"] = 0
    per_sample["Total"] = per_sample["Lobular"] + per_sample["Mucinous"]
    per_sample = per_sample[["Sample", "Lobular", "Mucinous", "Total"]].sort_values("Sample")

    # Print a quick summary to stdout
    print("== Overall matches ==")
    for _, r in overall.sort_values("Group").iterrows():
        print(f"{r['Group']}: {r['Count']}")
    print("\n== Per-sample matches ==")
    for _, r in per_sample.iterrows():
        print(f"{r['Sample']}: Lobular={r['Lobular']}  Mucinous={r['Mucinous']}  Total={r['Total']}")

    # Optional outputs
    if args.out_overall:
        overall.to_csv(args.out_overall, sep="\t", index=False)
    if args.out_per_sample:
        per_sample.to_csv(args.out_per_sample, sep="\t", index=False)

if __name__ == "__main__":
    main()
