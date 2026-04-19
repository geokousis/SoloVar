#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd


def load_genes(gene_list_path: str) -> set[str]:
    """
    Parse genes from a mixed-format list file.
    Supports tokens separated by commas, tabs, or spaces.
    """
    text = Path(gene_list_path).read_text(encoding="utf-8", errors="ignore")
    tokens = re.split(r"[,\t\r\n ]+", text)
    genes = {t.strip().upper() for t in tokens if t.strip()}
    return genes


def pick_existing(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [c for c in cols if c in df.columns]


def main():
    ap = argparse.ArgumentParser(
        description="Filter variants by genes in all.list and include IMPACT/SIFT/PolyPhen."
    )
    ap.add_argument("--infile", required=True, help="Input TSV (e.g., 2F_slim_pure.tsv)")
    ap.add_argument("--gene_list", default="all.list", help="Gene list file path")
    ap.add_argument("--outfile", required=True, help="Output TSV")
    ap.add_argument(
        "--keep_all_columns",
        action="store_true",
        help="Keep all input columns (still filtered by genes).",
    )
    args = ap.parse_args()

    genes = load_genes(args.gene_list)
    if not genes:
        raise SystemExit(f"No genes parsed from gene list: {args.gene_list}")

    df = pd.read_csv(args.infile, sep="\t", low_memory=False)
    if "Hugo_Symbol" not in df.columns:
        raise SystemExit("Input TSV is missing required column: Hugo_Symbol")

    mask = df["Hugo_Symbol"].astype(str).str.upper().isin(genes)
    out = df[mask].copy()

    if not args.keep_all_columns:
        requested = [
            "Tumor_Sample_Barcode",
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Genome_Change",
            "Codon_Change",
            "Protein_Change",
            "Variant_Classification",
            "Variant_Type",
            "HGVSc",
            "HGVSp",
            "HGVSp_Short",
            "Transcript_ID",
            "Exon_Number",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "t_ref_count",
            "t_alt_count",
            "VAF",
            "ONCOGENIC",
            "hotspot",
            "IMPACT",
            "SIFT",
            "PolyPhen",
        ]
        keep_cols = pick_existing(out, requested)
        missing_cols = [c for c in requested if c not in out.columns]
        out = out[keep_cols]

    out.to_csv(args.outfile, sep="\t", index=False)
    print(f"parsed_genes={len(genes)}")
    print(f"input_rows={len(df)}")
    print(f"kept_rows={len(out)}")
    if not args.keep_all_columns:
        print(f"kept_columns={','.join(keep_cols)}")
        print(f"missing_requested_columns={','.join(missing_cols) if missing_cols else 'NONE'}")
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
