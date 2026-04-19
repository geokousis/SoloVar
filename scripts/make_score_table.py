#!/usr/bin/env python3
import argparse
import re
import numpy as np
import pandas as pd


def pick(df, *cands):
    for c in cands:
        if c in df.columns:
            return c
    return None


def num_in_parens_or_self(x):
    if pd.isna(x):
        return pd.NA
    s = str(x)
    m = re.search(r"\(([-+]?\d*\.?\d+)\)", s)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return pd.NA
    try:
        return float(s)
    except Exception:
        return pd.NA


def main():
    ap = argparse.ArgumentParser(
        description="Build IMPACT/SIFT/PolyPhen score table from a TSV."
    )
    ap.add_argument("--infile", required=True, help="Input TSV")
    ap.add_argument(
        "--outfile",
        default="impact_sift_polyphen_scores.tsv",
        help="Output scored TSV",
    )
    ap.add_argument(
        "--outsummary",
        default="impact_sift_polyphen_scores_summary.tsv",
        help="Output summary TSV",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t", low_memory=False)

    col_imp = pick(df, "IMPACT")
    col_sift = pick(df, "SIFT")
    col_poly = pick(df, "PolyPhen")

    # IMPACT: HIGH=+2, MODERATE=+1
    if col_imp:
        imp = df[col_imp].astype(str).str.upper()
        w_impact = np.where(imp.eq("HIGH"), 2.0, np.where(imp.eq("MODERATE"), 1.0, 0.0))
    else:
        w_impact = np.zeros(len(df), dtype=float)

    # SIFT: text-only scoring
    if col_sift:
        sift_txt = df[col_sift].astype(str)
        sift_low = sift_txt.str.contains(r"deleterious_low_confidence", case=False, na=False)
        sift_del = sift_txt.str.contains(r"\bdeleterious\b", case=False, na=False) & ~sift_low
        w_sift = (sift_del.astype(float) * 1.0) + (sift_low.astype(float) * 0.5)
    else:
        w_sift = np.zeros(len(df), dtype=float)

    # PolyPhen scoring
    if col_poly:
        poly_txt = df[col_poly].astype(str)
        poly_sc = pd.to_numeric(poly_txt.map(num_in_parens_or_self), errors="coerce")
        poly_prob = poly_txt.str.contains(r"\bprobably_damaging\b", case=False, na=False) | (
            poly_sc.fillna(0.0) >= 0.85
        )
        poly_poss = poly_txt.str.contains(r"\bpossibly_damaging\b", case=False, na=False) | (
            poly_sc.between(0.15, 0.85, inclusive="both").fillna(False)
        )
        w_poly = (poly_prob.astype(float) * 1.0) + ((~poly_prob & poly_poss).astype(float) * 0.5)
    else:
        w_poly = np.zeros(len(df), dtype=float)

    score3_total = w_impact + w_sift + w_poly
    score3_vs_filter_threshold = score3_total / 2.0
    score3_norm_max4 = score3_total / 4.0

    mut_cols_pref = [
        "Tumor_Sample_Barcode",
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
    ]
    mut_cols = [c for c in mut_cols_pref if c in df.columns]

    out = df[mut_cols].copy()
    if col_imp:
        out["IMPACT"] = df[col_imp]
    if col_sift:
        out["SIFT"] = df[col_sift]
    if col_poly:
        out["PolyPhen"] = df[col_poly]

    out["w_impact"] = w_impact
    out["w_sift"] = w_sift
    out["w_polyphen"] = w_poly
    out["score3_total"] = score3_total
    out["score3_vs_filter_threshold"] = score3_vs_filter_threshold
    out["score3_norm_max4"] = score3_norm_max4
    out.to_csv(args.outfile, sep="\t", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "n_mutations",
                "mean_score3_total",
                "median_score3_total",
                "mean_score3_vs_filter_threshold",
                "median_score3_vs_filter_threshold",
                "mean_score3_norm_max4",
                "median_score3_norm_max4",
            ],
            "value": [
                len(out),
                float(np.mean(score3_total)),
                float(np.median(score3_total)),
                float(np.mean(score3_vs_filter_threshold)),
                float(np.median(score3_vs_filter_threshold)),
                float(np.mean(score3_norm_max4)),
                float(np.median(score3_norm_max4)),
            ],
        }
    )
    summary.to_csv(args.outsummary, sep="\t", index=False)

    print(f"wrote {args.outfile}")
    print(f"wrote {args.outsummary}")


if __name__ == "__main__":
    main()
