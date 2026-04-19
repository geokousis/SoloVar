#!/usr/bin/env python3
"""
Filter a MAF file by excluding variants that are present in one or more known VCF files.
The script accepts command-line arguments:
    --maf-input  : Input MAF file (tab-delimited).
    --vcf-file   : One or more VCF file paths (can be provided multiple times).
    --output     : Name/path for the filtered output MAF file.

Filtering is done in two steps:
  1. Remove low-confidence/common variants based on the "dbSNP_Val_Status_func" field.
  2. Remove variants that match any entry in the known VCF files, using a variant key
     constructed as "Chromosome_Start_Position_Reference_Allele_func_Tumor_Seq_Allele2_func".
"""

import argparse
import sys
import pandas as pd
import pysam

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Filter MAF file by removing variants found in provided VCF files."
    )
    parser.add_argument(
        "--maf-input",
        required=True,
        help="Input MAF file (tab-delimited, e.g., combined_merged_with_cosmic.maf)"
    )
    parser.add_argument(
        "--vcf-file",
        action="append",
        required=True,
        help="Path to a known VCF file. Can be specified multiple times for multiple VCF files."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path for the output filtered MAF file."
    )
    return parser.parse_args()

def load_known_variants(vcf_files):
    """
    Loads known variants from the provided VCF files. Each variant is stored as a string key:
    "contig_pos_ref_alt"
    """
    known_variants = set()
    for vcf_file in vcf_files:
        print(f"Loading known variants from: {vcf_file}")
        try:
            vcf = pysam.VariantFile(vcf_file)
        except Exception as e:
            print(f"Error loading VCF file {vcf_file}: {e}", file=sys.stderr)
            continue

        for record in vcf:
            if not record.alts:  # Skip records without alternate alleles
                continue
            for alt in record.alts:
                # Create a string key that combines contig, pos, ref, and alt.
                key = f"{record.contig}_{record.pos}_{record.ref}_{alt}"
                known_variants.add(key)
    print(f"Total known variants loaded from VCFs: {len(known_variants)}")
    return known_variants

def filter_maf(maf_file, known_variants):
    """
    Loads the MAF file, filters out low-confidence/common variants as well as any variant
    that is present in the known_variants set.
    """
    # Load the MAF file (ignoring any header lines beginning with '#').
    try:
        df = pd.read_csv(maf_file, sep="\t", dtype=str, low_memory=False)
    except Exception as e:
        print(f"Error reading MAF file {maf_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Filter out low-confidence/common variants, if the column exists.
    if "dbSNP_Val_Status_func" in df.columns:
        df = df[~df["dbSNP_Val_Status_func"].astype(str).str.contains("byFrequency|by1000genomes", na=False)]
    else:
        print("Warning: 'dbSNP_Val_Status_func' column not found in MAF. Skipping frequency filtering.")

    # Ensure required columns are present.
    required_cols = ["Chromosome", "Start_Position", "Reference_Allele_func", "Tumor_Seq_Allele2_func"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in MAF file: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    # Drop rows with missing data in any of the required columns.
    df = df.dropna(subset=required_cols)

    # Create a vectorized variant key for each row:
    # Concatenate the values of the required columns into a single key string.
    df["variant_key"] = (
        df["Chromosome"].astype(str) + "_" +
        df["Start_Position"].astype(int).astype(str) + "_" +
        df["Reference_Allele_func"].astype(str) + "_" +
        df["Tumor_Seq_Allele2_func"].astype(str)
    )

    before_count = len(df)
    # Remove rows where the variant key is found in known_variants.
    df = df[~df["variant_key"].isin(known_variants)]
    after_count = len(df)

    print(f"Variants before VCF filtering: {before_count}")
    print(f"Variants after VCF filtering: {after_count}")

    # Drop the helper column if not needed in the output.
    df = df.drop(columns=["variant_key"])
    return df

def main():
    args = parse_arguments()

    # Load variants from VCF files.
    known_variants = load_known_variants(args.vcf_file)

    # Filter the MAF file.
    filtered_df = filter_maf(args.maf_input, known_variants)

    # Save the filtered MAF file.
    try:
        filtered_df.to_csv(args.output, sep="\t", index=False)
        print(f"Filtered MAF file saved as: {args.output}")
    except Exception as e:
        print(f"Error saving filtered MAF file to {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
