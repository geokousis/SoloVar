#!/usr/bin/env python3
"""
slim_maf.py — Coalesce _onko/_func columns, compute VAF, extract COSMIC IDs,
               parse tumor-type summary, and write a slim TSV.

Usage:
    python3 slim_maf.py --input pf_af_filtered.maf --output af_filtered_f_pure.slim.tsv \
        [--tumor-types-out af_filtered_f_pure.tumor_types_long.tsv] \
        [--summary-dir slim_summaries/]
"""
import argparse
import re
import os
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description="Coalesce MAF columns and produce a slim TSV.")
    ap.add_argument("--input",  required=True, help="Input MAF (tab-separated)")
    ap.add_argument("--output", required=True, help="Output slim TSV")
    ap.add_argument("--tumor-types-out", default=None,
                    help="Optional: long-form tumor-type table (default: <output_stem>.tumor_types_long.tsv)")
    ap.add_argument("--summary-dir", default=None,
                    help="Directory for unique-value summary files (default: alongside --output)")
    return ap.parse_args()


# -------------------------
# Utilities
# -------------------------
def coalesce(df: pd.DataFrame, base: str):
    """Prefer *_onko > *_func > base. Merge into base column, drop alternates."""
    candidates = [f"{base}_onko", f"{base}_func", base]
    present = [c for c in candidates if c in df.columns]
    if not present:
        return None
    s = pd.Series(pd.NA, index=df.index, dtype="object")
    for c in [base, f"{base}_func", f"{base}_onko"]:
        if c in df.columns:
            s = s.where(s.notna(), df[c])
    df[base] = s
    for c in candidates:
        if c in df.columns and c != base:
            df.drop(columns=c, inplace=True)
    return base


COSMIC_PAT = re.compile(r"(?:COSM|COSV)\d+", flags=re.IGNORECASE)

def extract_cosmic_ids_from_str(x):
    if pd.isna(x):
        return []
    vals = [v.upper() for v in COSMIC_PAT.findall(str(x))]
    seen, uniq = set(), []
    for v in vals:
        if v not in seen:
            uniq.append(v); seen.add(v)
    return uniq


def parse_tumor_type_summary(s: str):
    """'biliary_tract(2)|breast(11)' -> [(type, count), ...]"""
    if pd.isna(s):
        return []
    out = []
    for p in str(s).split("|"):
        m = re.match(r"\s*([A-Za-z0-9_\-/\s]+)\((\d+)\)\s*$", p)
        if m:
            out.append((m.group(1).strip(), int(m.group(2))))
    return out


def sanitize_filename(name: str) -> str:
    s = name.strip().lower()
    s = re.sub(r"[^a-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "column"


def write_unique_column(df, col, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.join(out_dir, f"{sanitize_filename(col)}.txt")
    with open(fname, "w") as f:
        if col not in df.columns:
            f.write(f"[{col}] not found\n")
            return fname
        s = (df[col].dropna().astype(str).str.strip()
               .replace("", pd.NA).dropna().unique().tolist())
        s_sorted = sorted(set(s))
        f.write(f"# {col} — {len(s_sorted)} unique values\n")
        for v in s_sorted:
            f.write(f"{v}\n")
    return fname


# -------------------------
# Main
# -------------------------
def main():
    args = parse_args()

    MAF_IN       = args.input
    SLIM_OUT     = args.output
    stem         = os.path.splitext(SLIM_OUT)[0]
    TT_LONG_OUT  = args.tumor_types_out or f"{stem}.tumor_types_long.tsv"
    SUMMARY_DIR  = args.summary_dir or os.path.dirname(os.path.abspath(SLIM_OUT))

    # 1) Load
    try:
        maf = pd.read_csv(MAF_IN, sep="\t", low_memory=False)
    except Exception:
        maf = pd.read_csv(MAF_IN, sep="\t", engine="python", low_memory=False)

    # 2) Coalesce core fields
    for b in [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Strand",
        "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
        "Protein_Change", "Codon_Change",
        "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
        "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "AF",
        "ONCOGENIC", "HIGHEST_LEVEL", "HIGHEST_SENSITIVE_LEVEL", "HIGHEST_RESISTANCE_LEVEL",
        "COSMIC_tissue_types_affected", "MUTATION_SOMATIC_STATUS",
    ]:
        coalesce(maf, b)

    maf.rename(columns={
        "hgvsc": "HGVSc", "hgvsp": "HGVSp", "hgvsp_short": "HGVSp_Short",
        "transcript_id": "Transcript_ID", "exon_number": "Exon_Number",
        "filter": "FILTER", "vcf_qual": "vcf_qual", "vcf_id": "vcf_id",
    }, inplace=True)

    # 3) VAF
    maf["t_alt_count"] = pd.to_numeric(maf.get("t_alt_count"), errors="coerce")
    maf["t_ref_count"] = pd.to_numeric(maf.get("t_ref_count"), errors="coerce")
    denom = (maf["t_alt_count"].fillna(0) + maf["t_ref_count"].fillna(0)).replace(0, pd.NA)
    maf["VAF"] = maf["t_alt_count"] / denom

    # 4) COSMIC IDs
    cosmic_columns = [c for c in maf.columns if "COSMIC" in c.upper()]
    if "Existing_variation" in maf.columns:
        cosmic_columns.append("Existing_variation")
    cosmic_columns = list(dict.fromkeys(cosmic_columns))

    def gather_cosmic_ids(row):
        ids = []
        for c in cosmic_columns:
            if c in maf.columns:
                ids += extract_cosmic_ids_from_str(row.get(c, pd.NA))
        seen, uniq = set(), []
        for v in ids:
            if v not in seen:
                uniq.append(v); seen.add(v)
        return uniq

    if cosmic_columns:
        cosmic_list = maf.apply(gather_cosmic_ids, axis=1)
        maf["COSMIC_IDS"]       = cosmic_list.apply(lambda xs: ";".join(xs) if xs else pd.NA)
        maf["COSMIC_ID_PRIMARY"] = cosmic_list.apply(lambda xs: xs[0] if xs else pd.NA)

    # 5) Tumor-type summary
    tt_candidates = [c for c in maf.columns if "TUMOR_TYPE_SUMMARY" in c.upper()]
    tumor_type_col = tt_candidates[0] if tt_candidates else None

    if tumor_type_col:
        if tumor_type_col != "TUMOR_TYPE_SUMMARY":
            maf.rename(columns={tumor_type_col: "TUMOR_TYPE_SUMMARY"}, inplace=True)
        parsed = maf["TUMOR_TYPE_SUMMARY"].apply(parse_tumor_type_summary)

        def top_tt(lst):
            if not lst:
                return pd.Series({"TUMOR_TYPE_TOP": pd.NA, "TUMOR_TYPE_TOP_COUNT": pd.NA})
            lst_sorted = sorted(lst, key=lambda x: (-x[1], x[0]))
            return pd.Series({"TUMOR_TYPE_TOP": lst_sorted[0][0], "TUMOR_TYPE_TOP_COUNT": lst_sorted[0][1]})

        top_df = parsed.apply(top_tt)
        maf = pd.concat([maf, top_df], axis=1)

        keys = [k for k in ["Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short", "COSMIC_ID_PRIMARY"]
                if k in maf.columns]
        long_rows = []
        for i, lst in enumerate(parsed):
            if not lst:
                continue
            base = {k: maf.loc[i, k] for k in keys}
            for (tt, cnt) in lst:
                long_rows.append({**base, "Tumor_Type": tt, "Tumor_Type_Count": cnt})
        if long_rows:
            tt_long = pd.DataFrame(long_rows)
            tt_long.to_csv(TT_LONG_OUT, sep="\t", index=False)
            print(f"Saved tumor-type long table: {TT_LONG_OUT}  rows={len(tt_long)}")

    # 6) Build slim table
    keep_cols = [
        "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Genome_Change", "Codon_Change", "Protein_Change",
        "Variant_Classification", "Variant_Type", "HGVSc", "HGVSp", "HGVSp_Short",
        "Transcript_ID", "Exon_Number",
        "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
        "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "VAF",
        "FILTER", "vcf_id", "vcf_qual", "dbSNP_RS",
        "gnomADe_AF", "AF",
        "ONCOGENIC", "HIGHEST_LEVEL", "HIGHEST_SENSITIVE_LEVEL", "HIGHEST_RESISTANCE_LEVEL",
        "MUTATION_EFFECT", "MUTATION_EFFECT_DESCRIPTION",
        "COSMIC_ID_PRIMARY", "COSMIC_IDS", "COSMIC_overlapping_mutations", "MUTATION_SOMATIC_STATUS",
        "TUMOR_TYPE_SUMMARY", "TUMOR_TYPE_TOP", "TUMOR_TYPE_TOP_COUNT",
        "COSMIC_tissue_types_affected",
        "hotspot",
        "IMPACT", "SIFT", "PolyPhen",
        "CGC_Mutation_Type", "CGC_Tumor_Types_Somatic",
        "COSMIC_total_alterations_in_gene",
    ]
    keep_cols = [c for c in keep_cols if c in maf.columns]
    slim = maf[keep_cols].copy()

    for poscol in ["Start_Position", "End_Position"]:
        if poscol in slim.columns:
            slim[poscol] = pd.to_numeric(slim[poscol], errors="coerce").astype("Int64")

    slim.to_csv(SLIM_OUT, sep="\t", index=False)
    print(f"Saved slim table: {SLIM_OUT}  rows={len(slim)}  cols={len(slim.columns)}")

    # 7) Unique-value summaries
    qc_dir   = os.path.join(SUMMARY_DIR, "qc_filtering")
    db_dir   = os.path.join(SUMMARY_DIR, "common_dbs")
    onk_dir  = os.path.join(SUMMARY_DIR, "oncokb")

    written = []
    for c in ["FILTER", "vcf_id", "vcf_qual"]:
        written.append(write_unique_column(slim, c, qc_dir))
    for c in ["dbSNP_RS", "gnomADe_AF", "AF"]:
        written.append(write_unique_column(slim, c, db_dir))
    for c in ["ONCOGENIC", "HIGHEST_LEVEL", "HIGHEST_SENSITIVE_LEVEL", "HIGHEST_RESISTANCE_LEVEL"]:
        written.append(write_unique_column(slim, c, onk_dir))

    print("Wrote unique-value summaries:")
    for p in written:
        print("  -", p)


if __name__ == "__main__":
    main()
