#!/usr/bin/env python3
"""
Summarize gene-level CNV status from CNVkit *.call.cns files
based on a custom gene list (main gene + related genes), using
the AVERAGE copy number across all overlapping segments for each
(gene, sample) before classifying.

Each line of the gene list file looks like:
  CDH1    JUP,KLRG1,CTNND1,CTNNB1,CTNNA1,CTNNA2,CTNNA3,RAPGEF1,FLOT1,FLOT2,SLC6A3
  CIC     No-String
  ASAP3   ASAP3,ARF6,ASAP2,ASAP1

Output columns:
  MainGene   GeneGroup   Gene   Sample   EventType

Usage:
  python gene_level_cnv_summary_avg.py \
    --cnvkit-dir /path/to/calls \
    --gene-list gene.list \
    --out gene_cnv_summary.tsv
"""

import argparse
import os
from collections import defaultdict
import pandas as pd

# ---------------------------
# CNV thresholds (inclusive at boundaries)
# ---------------------------
CN_THRESHOLDS = {
    "HOMO_DEL": 0.5,  # CN <= 0.5
    "HEMI_DEL": 1.0,  # 0.5 < CN <= 1.0
    "LOW_AMP": 6,     # 6 <= CN < 10
    "HIGH_AMP": 10,   # CN >= 10
    "GAIN": 5         # 5 <= CN < 6
}

# ---------------------------
# Helpers
# ---------------------------
def classify_cn(cn: float):
    """Return CNV event label from an averaged CN value using inclusive bounds."""
    if pd.isna(cn):
        return None
    if cn <= CN_THRESHOLDS["HOMO_DEL"]:
        return "Homozygous Deletion"
    elif cn <= CN_THRESHOLDS["HEMI_DEL"]:
        return "Hemizygous Deletion"
    elif cn >= CN_THRESHOLDS["HIGH_AMP"]:
        return "High Amplification"
    elif cn >= CN_THRESHOLDS["LOW_AMP"]:
        return "Low-level Amplification"
    elif cn >= CN_THRESHOLDS["GAIN"]:
        return "Gain"
    else:
        return None  # neutral

def load_gene_list(path):
    """Load gene list mapping main genes to related genes."""
    gene_groups = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            main, partners = parts
            if partners == "No-String":
                gene_groups[main] = [main]
            else:
                genes = [g.strip() for g in partners.split(",") if g.strip()]
                # always include main gene itself
                if main not in genes:
                    genes.append(main)
                gene_groups[main] = genes
    return gene_groups

def read_cns(path):
    """Read CNVkit call.cns file."""
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = [c.strip().lower() for c in df.columns]
    rename_map = {"chromosome": "chrom", "chr": "chrom"}
    df.rename(columns=rename_map, inplace=True)
    return df

# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser(description="Summarize gene-level CNVs from CNVkit calls (averaged per gene).")
    ap.add_argument("--cnvkit-dir", required=True, help="Directory with *.call.cns files")
    ap.add_argument("--gene-list", required=True, help="Gene list file (2-column: main_gene partners)")
    ap.add_argument("--out", required=True, help="Output TSV file")
    args = ap.parse_args()

    gene_groups = load_gene_list(args.gene_list)
    print(f"[INFO] Loaded {len(gene_groups)} main genes from {args.gene_list}")

    # Build reverse index once: gene -> [main genes it belongs to]
    gene_to_mains = defaultdict(list)
    for main, targets in gene_groups.items():
        for g in targets:
            gene_to_mains[g].append(main)

    results = []

    for fn in sorted(os.listdir(args.cnvkit_dir)):
        if not fn.endswith(".call.cns"):
            continue
        sample = fn.replace(".call.cns", "")
        path = os.path.join(args.cnvkit_dir, fn)
        df = read_cns(path)

        if "gene" not in df.columns or "cn" not in df.columns:
            print(f"[WARN] Missing 'gene' or 'cn' columns in {fn}, skipping.")
            continue

        # Collect all CNs per gene for this sample (simple mean of all segments that list the gene)
        gene_cn_values = defaultdict(list)

        for _, row in df.iterrows():
            if pd.isna(row.get("gene")) or pd.isna(row.get("cn")):
                continue
            # Normalize CN to float; skip if it fails
            try:
                cn = float(row["cn"])
            except Exception:
                continue

            # CNVkit may list multiple comma-separated genes per segment
            genes = [g.strip() for g in str(row["gene"]).split(",") if g.strip()]
            for g in genes:
                if g in gene_to_mains:  # only track genes we care about
                    gene_cn_values[g].append(cn)

        # Average and classify once per gene
        for g, vals in gene_cn_values.items():
            if not vals:
                continue
            avg_cn = sum(vals) / len(vals)
            event = classify_cn(avg_cn)
            if not event:
                continue  # neutral after averaging

            for main_gene in gene_to_mains[g]:
                results.append({
                    "MainGene": main_gene,
                    "GeneGroup": main_gene,  # same for clarity
                    "Gene": g,
                    "Sample": sample,
                    "EventType": event
                })

    if not results:
        print("[WARN] No CNV events found for target genes.")
        return

    out_df = pd.DataFrame(results)
    out_df.sort_values(["MainGene", "Gene", "Sample"], inplace=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[DONE] Saved gene-level CNV summary to {args.out}")

if __name__ == "__main__":
    main()
