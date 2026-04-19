#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

def to_num(s):
    return pd.to_numeric(s, errors="coerce")

def coerce_bool(series):
    """
    Normalize boolean-like columns (True/False/1/0/yes/no/t/f -> boolean dtype).
    Returns a BooleanDtype series with True/False/NA.
    """
    if series is None:
        return None
    s = series.astype(str).str.strip().str.lower()
    true_vals  = {"true","t","1","yes","y"}
    false_vals = {"false","f","0","no","n"}
    out = pd.Series(pd.NA, index=series.index, dtype="boolean")
    out[s.isin(true_vals)]  = True
    out[s.isin(false_vals)] = False
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input MAF/TSV")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--report", required=True, help="Filter report TSV with counts")
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t", low_memory=False)
    start_n = len(df)

    # --- bookkeeping for the report
    rows = []
    def add_step(name, colname, before_df, after_df):
        rows.append({
            "step": name,
            "column": colname or "",
            "before": len(before_df),
            "removed": len(before_df) - len(after_df),
            "kept": len(after_df)
        })

    rows.append({"step": "initial", "column": "", "before": start_n, "removed": 0, "kept": start_n})
    cur = df

    # ---- Step 1: AF <= 0.001 (if present)
    af_cols = [c for c in ["AF","af","AF_func","AF_onko"] if c in cur.columns]
    if af_cols:
        c = af_cols[0]
        vals = to_num(cur[c])
        before = cur
        cur = cur[(vals.isna()) | (vals <= 0.001)]
        add_step("max_AF_0.001", c, before, cur)

    # ---- Step 2: gnomADe_AF <= 0.001 (if present)
    gnom_cols = [c for c in ["gnomADe_AF","gnomad_af","gnomad_exome_af","gnomad_genome_af"] if c in cur.columns]
    if gnom_cols:
        c = gnom_cols[0]
        vals = to_num(cur[c])
        before = cur
        cur = cur[(vals.isna()) | (vals <= 0.001)]
        add_step("max_gnomAD_AF_0.001", c, before, cur)

    # ---- Step 3: FILTER/ML presence logic + CONFLICT annotation
    # FILTER normalized
    if "FILTER" in cur.columns:
        s_filter = cur["FILTER"].astype(str).str.strip().str.upper()
        filter_is_pass = s_filter.isin({"PASS", ".", ""})
    else:
        # if FILTER missing, treat all as pass-through for conflict calc
        s_filter = pd.Series("", index=cur.index)
        filter_is_pass = pd.Series(True, index=cur.index)

    # ML column (PureCN) detection
    som_cols = [c for c in ["PureCN_ML_SOMATIC","ML.SOMATIC","ml.somatic"] if c in cur.columns]
    if som_cols:
        ml_col = som_cols[0]
        som_bool = coerce_bool(cur[ml_col])
    else:
        ml_col = None
        som_bool = pd.Series(pd.NA, index=cur.index, dtype="boolean")

    # NEW: drop rows with failing FILTER but no ML info
    before = cur
    if ml_col is None:
        # no ML column at all → drop any non-PASS FILTER
        cur = cur[filter_is_pass].copy()
        add_step("drop_nonPASS_without_ML", "FILTER (no ML column)", before, cur)
        # recompute s_filter/filter_is_pass to match new cur
        s_filter = cur["FILTER"].astype(str).str.strip().str.upper() if "FILTER" in cur.columns else pd.Series("", index=cur.index)
        filter_is_pass = s_filter.isin({"PASS", ".", ""})
        som_bool = pd.Series(pd.NA, index=cur.index, dtype="boolean")
    else:
        # ML present → drop non-PASS where ML is NA
        mask_keep = filter_is_pass | som_bool.notna()
        cur = cur[mask_keep].copy()
        add_step("drop_nonPASS_without_ML", f"FILTER & {ml_col}", before, cur)
        # refresh aligned views
        s_filter = cur["FILTER"].astype(str).str.strip().str.upper() if "FILTER" in cur.columns else pd.Series("", index=cur.index)
        filter_is_pass = s_filter.isin({"PASS", ".", ""})
        som_bool = coerce_bool(cur[ml_col])

    # CONFLICT = (ML True & FILTER non-PASS) OR (ML False & FILTER PASS)
    conflict = (
        (som_bool.fillna(False) & (~filter_is_pass)) |
        ((som_bool == False).fillna(False) & filter_is_pass)
    )
    cur = cur.copy()
    cur["CONFLICT"] = conflict.astype(bool)

    rows.append({
        "step": "annotate_CONFLICT",
        "column": f"{ml_col if ml_col else 'NA'} vs FILTER",
        "before": len(cur),
        "removed": 0,
        "kept": len(cur)
    })

    # ---- Save outputs
    cur.to_csv(args.out, sep="\t", index=False)

    # Final summary rows
    rows.append({
        "step": "final",
        "column": "",
        "before": start_n,
        "removed": start_n - len(cur),
        "kept": len(cur)
    })
    rep = pd.DataFrame(rows, columns=["step","column","before","removed","kept"])

    # add conflict count
    rep.loc[len(rep)] = {
        "step": "CONFLICT_count",
        "column": "CONFLICT==True",
        "before": len(cur),
        "removed": int(cur["CONFLICT"].sum()) if "CONFLICT" in cur.columns else 0,
        "kept": len(cur) - (int(cur["CONFLICT"].sum()) if "CONFLICT" in cur.columns else 0)
    }
    rep.to_csv(args.report, sep="\t", index=False)

    print(f"Input:  {args.inp}")
    print(f"Output: {args.out}")
    print(f"Report: {args.report}")
    print(f"Conflicts flagged: {int(cur['CONFLICT'].sum()) if 'CONFLICT' in cur.columns else 0}")

if __name__ == "__main__":
    main()
