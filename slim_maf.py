import re
import pandas as pd
import os

MAF_IN = "pf_af_filtered.maf"
SLIM_OUT = "af_filtered_f_pure.slim.tsv"
TT_LONG_OUT = "af_filtered_f_pure.tumor_types_long.tsv"

# -------------------------
# 1) Load
# -------------------------
try:
    maf = pd.read_csv(MAF_IN, sep="\t", low_memory=False)
except Exception:
    maf = pd.read_csv(MAF_IN, sep="\t", engine="python", low_memory=False)

# -------------------------
# 2) Utilities
# -------------------------
def coalesce(df: pd.DataFrame, base: str):
    """
    Prefer *_onko > *_func > base. Create/overwrite `base` with best-available,
    drop alternates.
    """
    candidates = [f"{base}_onko", f"{base}_func", base]
    present = [c for c in candidates if c in df.columns]
    if not present:
        return None
    out_col = base
    # start with base, overlay func, then onko
    s = pd.Series(pd.NA, index=df.index, dtype="object")
    for c in [base, f"{base}_func", f"{base}_onko"]:
        if c in df.columns:
            s = s.where(s.notna(), df[c])
    df[out_col] = s
    for c in candidates:
        if c in df.columns and c != out_col:
            df.drop(columns=c, inplace=True)
    return out_col

COSMIC_PAT = re.compile(r"(?:COSM|COSV)\d+", flags=re.IGNORECASE)
def extract_cosmic_ids_from_str(x):
    if pd.isna(x):
        return []
    vals = COSMIC_PAT.findall(str(x))
    # normalize case
    vals = [v.upper() for v in vals]
    # de-dup preserving order
    seen = set()
    uniq = []
    for v in vals:
        if v not in seen:
            uniq.append(v)
            seen.add(v)
    return uniq

def parse_tumor_type_summary(s: str):
    """
    Parse strings like 'biliary_tract(2)|breast(11)|central_nervous_system(44)'
    -> list of (type, count:int)
    """
    if pd.isna(s):
        return []
    parts = str(s).split("|")
    out = []
    for p in parts:
        m = re.match(r"\s*([A-Za-z0-9_\-/\s]+)\((\d+)\)\s*$", p)
        if m:
            out.append((m.group(1).strip(), int(m.group(2))))
    return out

# -------------------------
# 3) Coalesce core fields
# -------------------------
for b in [
    # core fields
    "Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Strand",
    "Variant_Classification","Variant_Type","Reference_Allele",
    "Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","dbSNP_Val_Status",
        "Protein_Change","Codon_Change",
    "Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode",
    "t_alt_count","t_ref_count","n_alt_count","n_ref_count","AF",
    "ONCOGENIC","HIGHEST_LEVEL","HIGHEST_SENSITIVE_LEVEL","HIGHEST_RESISTANCE_LEVEL",
        "COSMIC_tissue_types_affected","MUTATION_SOMATIC_STATUS"
]:
    coalesce(maf, b)

# Some MAFs use different casing for these
maf.rename(columns={
    "hgvsc":"HGVSc", "hgvsp":"HGVSp", "hgvsp_short":"HGVSp_Short",
    "transcript_id":"Transcript_ID","exon_number":"Exon_Number",
    "filter":"FILTER","vcf_qual":"vcf_qual","vcf_id":"vcf_id"
}, inplace=True)

# -------------------------
# 4) Compute VAF (tumor)
# -------------------------
# Make sure counts are numeric
maf["t_alt_count"] = pd.to_numeric(maf.get("t_alt_count"), errors="coerce")
maf["t_ref_count"] = pd.to_numeric(maf.get("t_ref_count"), errors="coerce")

denom = maf["t_alt_count"].fillna(0) + maf["t_ref_count"].fillna(0)
# avoid division warnings on zero denominator
denom = denom.replace(0, pd.NA)
maf["VAF"] = maf["t_alt_count"] / denom

# -------------------------
# 5) COSMIC IDs (COSM/COSV)
#    Try common columns + scan 'Existing_variation' if present
# -------------------------
cosmic_columns = [c for c in maf.columns if "COSMIC" in c.upper()]  # e.g., COSMIC_overlapping_mutations
if "Existing_variation" in maf.columns:
    cosmic_columns.append("Existing_variation")
cosmic_columns = list(dict.fromkeys(cosmic_columns))  # unique, keep order

def gather_cosmic_ids(row):
    ids = []
    for c in cosmic_columns:
        if c in maf.columns:
            ids += extract_cosmic_ids_from_str(row.get(c, pd.NA))
    # de-dup, keep order
    seen = set()
    uniq = []
    for v in ids:
        if v not in seen:
            uniq.append(v)
            seen.add(v)
    return uniq

if cosmic_columns:
    cosmic_list = maf.apply(gather_cosmic_ids, axis=1)
    maf["COSMIC_IDS"] = cosmic_list.apply(lambda xs: ";".join(xs) if xs else pd.NA)
    maf["COSMIC_ID_PRIMARY"] = cosmic_list.apply(lambda xs: xs[0] if xs else pd.NA)

# -------------------------
# 6) Tumor-type summary: detect + parse + derive
# -------------------------
tt_candidates = [c for c in maf.columns if "TUMOR_TYPE_SUMMARY" in c.upper()]
tumor_type_col = tt_candidates[0] if tt_candidates else None

if tumor_type_col:
    # Keep raw column under canonical name
    if tumor_type_col != "TUMOR_TYPE_SUMMARY":
        maf.rename(columns={tumor_type_col: "TUMOR_TYPE_SUMMARY"}, inplace=True)

    parsed = maf["TUMOR_TYPE_SUMMARY"].apply(parse_tumor_type_summary)

    # top tissue & count
    def top_tt(lst):
        if not lst:
            return pd.Series({"TUMOR_TYPE_TOP": pd.NA, "TUMOR_TYPE_TOP_COUNT": pd.NA})
        # sort by count desc, then keep first
        lst_sorted = sorted(lst, key=lambda x: (-x[1], x[0]))
        return pd.Series({"TUMOR_TYPE_TOP": lst_sorted[0][0], "TUMOR_TYPE_TOP_COUNT": lst_sorted[0][1]})

    top_df = parsed.apply(top_tt)
    maf = pd.concat([maf, top_df], axis=1)

    # Long/exploded table for plotting/aggregation
    long_rows = []
    # pull minimal keys to carry over
    keys = ["Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","COSMIC_ID_PRIMARY"]
    keys = [k for k in keys if k in maf.columns]
    for i, lst in enumerate(parsed):
        if not lst:
            continue
        base = {k: maf.loc[i, k] for k in keys}
        for (tt, cnt) in lst:
            row = dict(base)
            row["Tumor_Type"] = tt
            row["Tumor_Type_Count"] = cnt
            long_rows.append(row)
    if long_rows:
        tt_long = pd.DataFrame(long_rows)
        tt_long.to_csv(TT_LONG_OUT, sep="\t", index=False)
        print(f"Saved tumor-type long table: {TT_LONG_OUT}  rows={len(tt_long)}")

# -------------------------
# 7) Build SLIM table to save
# -------------------------
keep_cols = [
    # identifiers
    "Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Genome_Change","Codon_Change","Protein_Change",
    # variant impact/type
    "Variant_Classification","Variant_Type","HGVSc","HGVSp","HGVSp_Short","Transcript_ID","Exon_Number",
    # alleles/depth
    "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
    "t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","VAF",
    # QC / filters
    "FILTER","vcf_id","vcf_qual","dbSNP_RS",
    # population AF
    "gnomADe_AF","AF",
    # OncoKB minimal
    "ONCOGENIC","HIGHEST_LEVEL","HIGHEST_SENSITIVE_LEVEL","HIGHEST_RESISTANCE_LEVEL","MUTATION_EFFECT", "MUTATION_EFFECT_DESCRIPTION",
    # COSMIC
    "COSMIC_ID_PRIMARY","COSMIC_IDS","COSMIC_overlapping_mutations","MUTATION_SOMATIC_STATUS",
    # Tumor-type summary (raw + derived)
    "TUMOR_TYPE_SUMMARY","TUMOR_TYPE_TOP","TUMOR_TYPE_TOP_COUNT",
    # Tissue
    "COSMIC_tissue_types_affected",
    # HOTSPOT
    "hotspot",
    # impact / prediction
    "IMPACT", "SIFT", "PolyPhen",
    # gene-level cancer evidence
    "CGC_Mutation_Type", "CGC_Tumor_Types_Somatic",
    "COSMIC_total_alterations_in_gene",
    
    
]

keep_cols = [c for c in keep_cols if c in maf.columns]
slim = maf[keep_cols].copy()

# light cleanup
for poscol in ["Start_Position","End_Position"]:
    if poscol in slim.columns:
        slim[poscol] = pd.to_numeric(slim[poscol], errors="coerce").astype("Int64")

slim.to_csv(SLIM_OUT, sep="\t", index=False)
print(f"Saved slim table: {SLIM_OUT}  rows={len(slim)}  cols={len(slim.columns)}")



# Summary

# Use `maf` (full) or `slim` (your slimmed table). Pick one:
df_for_summary = slim  # or: slim

# ---------- helpers ----------
def sanitize_filename(name: str) -> str:
    # safe filename from column name
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
        # collect unique, non-null, trimmed strings
        s = (
            df[col]
            .dropna()
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .unique()
            .tolist()
        )
        s_sorted = sorted(set(s))
        f.write(f"# {col} — {len(s_sorted)} unique values\n")
        for v in s_sorted:
            f.write(f"{v}\n")
    return fname

# ---------- categories & columns ----------
qc_dir = "qc_filtering"
db_dir = "common_dbs"
onkob_dir = "oncokb"

qc_cols   = ["FILTER", "vcf_id", "vcf_qual"]
db_cols   = ["dbSNP_RS", "gnomADe_AF", "AF"]  # AF kept as fallback if gnomADe_AF missing
onk_cols  = ["ONCOGENIC", "HIGHEST_LEVEL", "HIGHEST_SENSITIVE_LEVEL", "HIGHEST_RESISTANCE_LEVEL"]

# ---------- write files ----------
written = []
for c in qc_cols:
    written.append(write_unique_column(df_for_summary, c, qc_dir))
for c in db_cols:
    written.append(write_unique_column(df_for_summary, c, db_dir))
for c in onk_cols:
    written.append(write_unique_column(df_for_summary, c, onkob_dir))

print("Wrote unique-value summaries:")
for p in written:
    print("  -", p)
