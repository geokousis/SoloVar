#!/usr/bin/env python3
import pandas as pd
import argparse
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="Path to filtered MAF/TSV file (e.g., af_filtered_f_pure.filtered.tsv)")
    ap.add_argument("--outdir", required=True, help="Output directory for results (TSVs)")
    ap.add_argument("--gene-column", default="Hugo_Symbol", help="Column with gene symbols (default: Hugo_Symbol)")
    ap.add_argument("--sample-column", default="Tumor_Sample_Barcode", help="Column with sample names")
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t", low_memory=False)

    # --- define sample → patient mapping ---
    sample_to_patient = {
        "MB1": "Patient_1", "MB7": "Patient_1",
        "MB4": "Patient_4", "MB9": "Patient_4",
        "MB2": "Patient_2", "MB3": "Patient_3",
        "MB5": "Patient_5", "MB6": "Patient_6",
        "MB8": "Patient_8", "MB10": "Patient_10",
    }
    sample_to_patient = {
        "MB1": "Patient_1", "MB7": "Patient_7",
        "MB4": "Patient_4", "MB9": "Patient_9",
        "MB2": "Patient_2", "MB3": "Patient_3",
        "MB5": "Patient_5", "MB6": "Patient_6",
        "MB8": "Patient_8", "MB10": "Patient_10",
    }

    df = df.copy()
    df[args.sample_column] = df[args.sample_column].astype(str)
    df[args.gene_column] = df[args.gene_column].astype(str)

    # map sample → patient
    df["Patient_ID"] = df[args.sample_column].map(sample_to_patient)
    df = df[df["Patient_ID"].notna()]  # drop any unmatched samples

    # --- collapse multiple hits per gene per patient ---
    unique_hits = df[[args.gene_column, "Patient_ID"]].drop_duplicates()

    # --- count mutated patients per gene ---
    gene_counts = (
        unique_hits.groupby(args.gene_column)["Patient_ID"]
        .nunique()
        .reset_index(name="Num_Patients_Mutated")
        .sort_values("Num_Patients_Mutated", ascending=False)
    )

    # --- per-patient summary: top mutated genes ---
    patient_gene_counts = (
        unique_hits.groupby("Patient_ID")[args.gene_column]
        .value_counts()
        .rename("Variant_Count")
        .reset_index()
    )
    top_genes_per_patient = (
        patient_gene_counts.groupby("Patient_ID")
        .apply(lambda x: x.nlargest(10, "Variant_Count"))
        .reset_index(drop=True)
    )

    # --- output files ---
    import os
    os.makedirs(args.outdir, exist_ok=True)
    path1 = f"{args.outdir}/mutated_gene_counts_overall.tsv"
    path2 = f"{args.outdir}/top10_genes_per_patient.tsv"

    gene_counts.to_csv(path1, sep="\t", index=False)
    top_genes_per_patient.to_csv(path2, sep="\t", index=False)

    print(f"Saved overall gene counts → {path1}")
    print(f"Saved top 10 mutated genes per patient → {path2}")
    print("\nTop overall mutated genes:")
    print(gene_counts.head(20).to_string(index=False))

if __name__ == "__main__":
    main()
