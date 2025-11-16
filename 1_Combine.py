#!/usr/bin/env python3
"""
Summarize gene-level SNVs from a MAF and (optionally) merge with CNV summary.

Inputs
------
- Gene list file: 2 columns (whitespace-delimited)
    <MainGene>   <comma_separated_related_genes_or_No-String>
  Same format as your CNV script. The main gene is always included in its group.

- MAF file: tab-delimited; key columns (case-insensitive):
    Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, HGVSp_Short (optional)

Outputs
-------
1) SNV summary TSV with columns:
   MainGene  GeneGroup  Gene  Sample  SNV_EventTypes  Protein_Changes

2) (Optional) Combined SNV+CNV TSV with columns:
   MainGene  GeneGroup  Gene  Sample  SNV_EventTypes  Protein_Changes  CNV_EventType  Status

   where Status ∈ {"Both","SNV only","CNV only"}.

Usage
-----
  # SNV-only summary
  python gene_level_snv_and_merge.py \
    --maf variants.maf \
    --gene-list gene.list \
    --snv-out gene_snv_summary.tsv \
    --snv-sample-regex '(?i)^(?:f_)?(?P<id>MB\\d+).*?$'

  # SNV + CNV combined summary (normalize CNV sample names too)
  python gene_level_snv_and_merge.py \
    --maf variants.maf \
    --gene-list gene.list \
    --snv-out gene_snv_summary.tsv \
    --cnv-summary gene_cnv_summary.tsv \
    --combined-out gene_snv_cnv_combined.tsv \
    --snv-sample-regex '(?i)^(?:f_)?(?P<id>MB\\d+).*?$' \
    --cnv-sample-regex '(?i)^(?:f_)?(?P<id>MB\\d+)(?:_recalibrated)?$'
"""

import argparse
import gzip
import io
import re
from collections import defaultdict
import pandas as pd


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
                if main not in genes:
                    genes.append(main)
                gene_groups[main] = genes
    return gene_groups


def _open_maybe_gzip(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return io.BytesIO(fh.read())
    else:
        return open(path, "rb")


def read_maf(path):
    """Read a (possibly gzipped) MAF and normalize column names to lowercase."""
    with _open_maybe_gzip(path) as fh:
        df = pd.read_csv(fh, sep="\t", comment="#", dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def choose_first_present(df, candidates):
    """Return the first present column in df from candidates (lowercased names)."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _build_normalizer(regex: str | None, fallback_strip_cn_wrappers: bool = False):
    """
    Return a function that normalizes a sample label using an optional regex.
    If regex has a named 'id' group, that is used; otherwise group(1).
    When fallback_strip_cn_wrappers=True, also strip leading 'f_' and trailing '_recalibrated' if regex doesn't match.
    """
    def norm(x):
        if pd.isna(x):
            return x
        s = str(x)
        if regex:
            m = re.search(regex, s)
            if m:
                return m.group("id") if "id" in m.groupdict() else (m.group(1) if m.groups() else s)
        if fallback_strip_cn_wrappers:
            s2 = re.sub(r"(?i)^f_", "", s)
            s2 = re.sub(r"(?i)_recalibrated$", "", s2)
            return s2
        return s
    return norm


def make_snv_summary(maf_df, gene_groups, snv_sample_regex: str | None):
    """
    Build per-(MainGene, Gene, Sample) SNV summary with event types and protein changes.
    """
    # Build reverse index: gene -> [main genes it belongs to]
    gene_to_mains = defaultdict(list)
    for main, targets in gene_groups.items():
        for g in targets:
            gene_to_mains[g].append(main)

    # Pick columns regardless of exact header casing
    col_sample = choose_first_present(maf_df, ["tumor_sample_barcode", "tumor_sample", "sample", "tumor_sample_id"])
    col_gene = choose_first_present(maf_df, ["hugo_symbol", "gene"])
    col_class = choose_first_present(maf_df, ["variant_classification", "vc"])
    col_prot = choose_first_present(maf_df, ["hgvsp_short", "hgvsp", "protein_change", "amino_acids"])

    if col_sample is None or col_gene is None or col_class is None:
        missing = [("Tumor_Sample_Barcode", col_sample),
                   ("Hugo_Symbol", col_gene),
                   ("Variant_Classification", col_class)]
        miss_txt = ", ".join([name for name, val in missing if val is None])
        raise RuntimeError(f"MAF is missing required column(s): {miss_txt}")

    normalize_sample = _build_normalizer(snv_sample_regex, fallback_strip_cn_wrappers=False)

    # Aggregate per (sample, gene)
    agg = {}
    for _, row in maf_df.iterrows():
        g = row[col_gene]
        if g not in gene_to_mains:
            continue
        s = normalize_sample(row[col_sample])
        vclass = row[col_class]
        prot = row[col_prot] if col_prot and pd.notna(row[col_prot]) else None

        key = (s, g)
        if key not in agg:
            agg[key] = {"types": set(), "proteins": set()}
        if pd.notna(vclass) and vclass and vclass != ".":
            agg[key]["types"].add(vclass)
        if prot and prot != ".":
            agg[key]["proteins"].add(prot)

    # Materialize rows
    rows = []
    for (sample, gene), vals in agg.items():
        snv_types = ",".join(sorted(vals["types"])) if vals["types"] else None
        prot_changes = ",".join(sorted(vals["proteins"])) if vals["proteins"] else None
        for main_gene in gene_to_mains[gene]:
            rows.append({
                "MainGene": main_gene,
                "GeneGroup": main_gene,
                "Gene": gene,
                "Sample": sample,
                "SNV_EventTypes": snv_types,
                "Protein_Changes": prot_changes
            })

    out = pd.DataFrame(rows)
    if not out.empty:
        out.sort_values(["MainGene", "Gene", "Sample"], inplace=True)
    return out


def merge_with_cnv(snv_df, cnv_path, cnv_sample_regex: str | None):
    """
    Merge SNV summary with CNV summary from the original script.

    Expects CNV TSV columns: MainGene, GeneGroup, Gene, Sample, EventType (or CNV_EventType)
    """
    cnv_df = pd.read_csv(cnv_path, sep="\t", dtype=str)

    # Normalize expected columns
    rename = {}
    for c in cnv_df.columns:
        lc = c.lower()
        if lc == "eventtype":
            rename[c] = "CNV_EventType"
        elif lc == "maingene":
            rename[c] = "MainGene"
        elif lc == "genegroup":
            rename[c] = "GeneGroup"
        elif lc == "gene":
            rename[c] = "Gene"
        elif lc == "sample":
            rename[c] = "Sample"
    if rename:
        cnv_df.rename(columns=rename, inplace=True)

    need = {"MainGene", "GeneGroup", "Gene", "Sample"}
    if not need.issubset(set(cnv_df.columns)):
        raise RuntimeError("CNV summary missing required columns: "
                           + ", ".join(sorted(need - set(cnv_df.columns))))

    if "CNV_EventType" not in cnv_df.columns:
        if "EventType" in cnv_df.columns:
            cnv_df.rename(columns={"EventType": "CNV_EventType"}, inplace=True)
        else:
            cnv_df["CNV_EventType"] = None

    # Normalize CNV 'Sample' (use regex; if not matched, also try f_/ _recalibrated stripping)
    cnv_df["Sample"] = cnv_df["Sample"].map(
        _build_normalizer(cnv_sample_regex, fallback_strip_cn_wrappers=True)
    )

    # Outer join so we see SNV-only / CNV-only / Both
    combined = pd.merge(
        snv_df, cnv_df,
        on=["MainGene", "GeneGroup", "Gene", "Sample"],
        how="outer"
    )

    # Derive Status
    has_snv = combined["SNV_EventTypes"].notna()
    has_cnv = combined["CNV_EventType"].notna()
    combined["Status"] = None
    combined.loc[ has_snv &  has_cnv, "Status"] = "Both"
    combined.loc[ has_snv & ~has_cnv, "Status"] = "SNV only"
    combined.loc[~has_snv &  has_cnv, "Status"] = "CNV only"

    combined.sort_values(["MainGene", "Gene", "Sample"], inplace=True)
    return combined


def main():
    ap = argparse.ArgumentParser(description="Summarize gene-level SNVs from a MAF and optionally merge with CNV summary.")
    ap.add_argument("--maf", required=True, help="Input MAF (tsv, optionally .gz)")
    ap.add_argument("--gene-list", required=True, help="Gene list file (same format as CNV script)")
    ap.add_argument("--snv-out", required=True, help="Output TSV for SNV summary")
    ap.add_argument("--cnv-summary", help="CNV summary TSV from CNV script (optional)")
    ap.add_argument("--combined-out", help="Output TSV for combined SNV+CNV (requires --cnv-summary)")
    ap.add_argument("--snv-sample-regex", help="Regex to extract sample id from MAF Tumor_Sample_Barcode (use named 'id' or group 1)")
    ap.add_argument("--cnv-sample-regex", help="Regex to extract sample id from CNV 'Sample' (use named 'id' or group 1)")
    args = ap.parse_args()

    gene_groups = load_gene_list(args.gene_list)
    print(f"[INFO] Loaded {len(gene_groups)} main genes from {args.gene_list}")

    maf_df = read_maf(args.maf)
    print(f"[INFO] Read MAF with {len(maf_df):,} rows from {args.maf}")

    snv_df = make_snv_summary(maf_df, gene_groups, snv_sample_regex=args.snv_sample_regex)
    if snv_df.empty:
        print("[WARN] No SNVs found for target genes.")
        # Emit empty file to keep pipelines happy
        pd.DataFrame(columns=["MainGene","GeneGroup","Gene","Sample","SNV_EventTypes","Protein_Changes"]).to_csv(args.snv_out, sep="\t", index=False)
    else:
        snv_df.to_csv(args.snv_out, sep="\t", index=False)
        print(f"[DONE] Saved SNV summary to {args.snv_out} ({len(snv_df):,} rows)")

    if args.cnv_summary and args.combined_out:
        combined = merge_with_cnv(snv_df, args.cnv_summary, cnv_sample_regex=args.cnv_sample_regex)
        combined.to_csv(args.combined_out, sep="\t", index=False)
        print(f"[DONE] Saved combined SNV+CNV summary to {args.combined_out} ({len(combined):,} rows)")
    elif args.cnv_summary and not args.combined_out:
        print("[INFO] --cnv-summary provided but no --combined-out; skipping merge.")
    elif args.combined_out and not args.cnv_summary:
        print("[WARN] --combined-out requested but --cnv-summary not provided; skipping merge.")


if __name__ == "__main__":
    main()
