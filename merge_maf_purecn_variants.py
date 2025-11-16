#!/usr/bin/env python3
import os, re, glob, argparse
import pandas as pd

# -----------------------
# helpers
# -----------------------
def norm_sample_name(x: str) -> str:
    if x is None:
        return x
    # f_MB10_recalibrated -> MB10
    m = re.search(r"f_(.+?)_recalibrated", x)
    if m:
        return m.group(1)
    x = re.sub(r"^f_", "", x)
    x = re.sub(r"_recalibrated$", "", x)
    return x

def list_sample_dirs(root: str):
    out = {}
    for d in glob.glob(os.path.join(root, "*")):
        if os.path.isdir(d):
            base = os.path.basename(d)
            out[norm_sample_name(base)] = d
    return out

def load_purecn_recalibrated_variants(sample_dir: str, sample_norm: str) -> pd.DataFrame:
    """
    Load *_recalibrated_variants.csv and:
      - normalize columns
      - add Sample_norm = sample_norm (from folder)
      - keep key PureCN fields (renamed to PureCN_*).
    """
    paths = glob.glob(os.path.join(sample_dir, "*_recalibrated_variants.csv"))
    if not paths:
        return pd.DataFrame()
    df = pd.read_csv(paths[0])

    # lower-case headers for detection
    lower = {c: c.lower() for c in df.columns}
    df.rename(columns=lower, inplace=True)

    # map essentials
    ren = {}
    # coordinates
    if "chr" in df.columns: ren["chr"] = "Chromosome"
    elif "chromosome" in df.columns: ren["chromosome"] = "Chromosome"
    elif "seqnames" in df.columns: ren["seqnames"] = "Chromosome"
    if "start" in df.columns: ren["start"] = "Start"
    if "end" in df.columns: ren["end"] = "End"

    # useful per-variant fields
    mapping = [
        ("ml.m", "PureCN_minorCN"),
        ("ml.c", "PureCN_majorCN"),
        ("ml.loh", "PureCN_LOH"),
        ("ml.somatic", "PureCN_ML_SOMATIC"),
        ("posterior.somatic", "PureCN_POSTERIOR_SOMATIC"),
        ("cellfraction", "PureCN_cell_fraction"),
        ("cellfraction.95.lower", "PureCN_cellfrac_95L"),
        ("cellfraction.95.upper", "PureCN_cellfrac_95U"),
        ("log.ratio", "PureCN_logR"),
        ("depth", "PureCN_depth"),
        ("ar", "PureCN_AR"),
        ("ar.adjusted", "PureCN_AR_adjusted"),
        ("mapping.bias", "PureCN_mapping_bias"),
        ("allelic.imbalance", "PureCN_allelic_imbalance"),
        ("flagged", "PureCN_flagged"),
        ("on.target", "PureCN_on_target"),
        ("gene.symbol", "Hugo_Symbol_PureCN"),
        ("ref", "REF_PureCN"),
        ("alt", "ALT_PureCN"),
        ("id", "VAR_ID_PureCN"),
        ("seg.id", "SEG_ID_PureCN"),
        ("sampleid", "Sampleid_original"),
    ]
    for src, dst in mapping:
        if src in df.columns:
            ren[src] = dst

    df = df.rename(columns=ren)

    # normalize values
    if "Chromosome" not in df.columns: df["Chromosome"] = pd.NA
    df["Chromosome"] = df["Chromosome"].astype(str).str.replace("^chr", "", regex=True)

    if "Start" not in df.columns:
        raise RuntimeError(f"No 'start' column in {paths[0]}")
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce").astype("Int64")

    # keep lean set for merge
    keep_cols = ["Chromosome","Start"] + [c for _, c in mapping if c in df.columns]
    df = df[[c for c in keep_cols if c in df.columns]].copy()

    # attach normalized sample from folder (authoritative for matching)
    df["Sample_norm"] = sample_norm
    return df

def sanitize_filename(name: str) -> str:
    s = name.strip().lower()
    s = re.sub(r"[^a-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "column"

def write_unique_values(series: pd.Series, path: str, float_round: int = 4):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        if series is None:
            f.write("[missing]\n"); return
        s = series.dropna()
        if s.empty:
            f.write("[no values]\n"); return
        # try to coerce numeric for nicer dedup
        try:
            s_num = pd.to_numeric(s, errors="coerce")
            if s_num.notna().any():
                s = s_num.round(float_round).astype(object).where(s_num.notna(), s)
        except Exception:
            pass
        vals = (
            s.astype(str).str.strip().replace("", pd.NA).dropna().unique().tolist()
        )
        vals = sorted(set(vals))
        f.write(f"# {os.path.basename(path)} — {len(vals)} unique values\n")
        for v in vals:
            f.write(f"{v}\n")

# -----------------------
# main
# -----------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slim-maf", required=True, help="Slim MAF TSV")
    ap.add_argument("--purecn-root", required=True, help="Root with f_*_recalibrated folders")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--uniq-dir", required=True, help="Directory to write PureCN unique-value summaries")
    args = ap.parse_args()

    # load slim maf
    maf = pd.read_csv(args.slim_maf, sep="\t", low_memory=False)
    maf["Chromosome"] = maf["Chromosome"].astype(str).str.replace("^chr", "", regex=True)
    maf["Start_Position"] = pd.to_numeric(maf["Start_Position"], errors="coerce").astype("Int64")
    maf["Tumor_Sample_Barcode"] = maf["Tumor_Sample_Barcode"].astype(str)

    # gather all PureCN variant rows across samples
    sample_dirs = list_sample_dirs(args.purecn_root)
    pcv_list = []
    for sample_norm, sdir in sample_dirs.items():
        df = load_purecn_recalibrated_variants(sdir, sample_norm)
        if not df.empty:
            pcv_list.append(df)
    if pcv_list:
        pcv_all = pd.concat(pcv_list, ignore_index=True)
    else:
        pcv_all = pd.DataFrame(columns=["Sample_norm","Chromosome","Start"])

    # do a single left merge on (sample, chromosome, start)
    purecn_cols = [c for c in pcv_all.columns if c not in ("Sample_norm","Chromosome","Start")]
    merged = maf.merge(
        pcv_all,
        how="left",
        left_on=["Tumor_Sample_Barcode","Chromosome","Start_Position"],
        right_on=["Sample_norm","Chromosome","Start"]
    )

    # drop helper keys if you don’t want them in the final table
    merged.drop(columns=[c for c in ["Sample_norm","Start"] if c in merged.columns], inplace=True)

    # save
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Saved merged table: {args.out}  rows={len(merged)}  cols={len(merged.columns)}")

    # ---------------------------
    # unique values per PureCN category
    # ---------------------------
    uniq_root = args.uniq_dir
    os.makedirs(uniq_root, exist_ok=True)

    # categories
    cn_dir   = os.path.join(uniq_root, "purecn_cn")
    prob_dir = os.path.join(uniq_root, "purecn_probabilities")
    loh_dir  = os.path.join(uniq_root, "purecn_loh")
    misc_dir = os.path.join(uniq_root, "purecn_misc")

    for col in ["PureCN_minorCN","PureCN_majorCN","PureCN_logR","PureCN_depth"]:
        if col in merged.columns:
            write_unique_values(merged[col], os.path.join(cn_dir, f"{sanitize_filename(col)}.txt"))
    for col in ["PureCN_ML_SOMATIC","PureCN_POSTERIOR_SOMATIC","PureCN_mapping_bias","PureCN_AR","PureCN_AR_adjusted"]:
        if col in merged.columns:
            write_unique_values(merged[col], os.path.join(prob_dir, f"{sanitize_filename(col)}.txt"))
    for col in ["PureCN_LOH","PureCN_allelic_imbalance"]:
        if col in merged.columns:
            write_unique_values(merged[col], os.path.join(loh_dir, f"{sanitize_filename(col)}.txt"))
    for col in ["PureCN_flagged","PureCN_on_target","SEG_ID_PureCN","Hugo_Symbol_PureCN"]:
        if col in merged.columns:
            write_unique_values(merged[col], os.path.join(misc_dir, f"{sanitize_filename(col)}.txt"))

    print("Wrote PureCN unique-value summaries under:", uniq_root)

if __name__ == "__main__":
    main()
