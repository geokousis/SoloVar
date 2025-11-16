#!/usr/bin/env python3
import argparse, re
import pandas as pd
import numpy as np

# =========================
# Config
# =========================
# Extra strict gate for noisy samples; rescue (OncoKB/hotspot) bypasses this.
STRICT_SAMPLES = {
    "MB5": (10, 40, 0.05),  # (min_alt, min_depth, min_vaf)
    "MB6": (10, 40, 0.05),
}

# Conflict override: allow strong base support alone to clear conflicts
CONFLICT_STRONG_VAF   = 0.30
CONFLICT_STRONG_DEPTH = 50

# =========================
# Helpers
# =========================
def pick(df, *cands):
    for c in cands:
        if c in df.columns:
            return c
    return None

def to_num(s):
    return pd.to_numeric(s, errors="coerce")

def coerce_bool(series):
    s = series.astype(str).str.strip().str.lower()
    true_vals  = {"true","t","1","yes","y"}
    false_vals = {"false","f","0","no","n"}
    out = pd.Series(pd.NA, index=series.index, dtype="boolean")
    out[s.isin(true_vals)]  = True
    out[s.isin(false_vals)] = False
    return out

# parse strings like "biliary_tract(2)|breast(11)|kidney(92)"
_pattern_token = re.compile(r"\s*([^()|]+)\((\d+)\)\s*")
def parse_cosmic_tissues(s: str):
    if not isinstance(s, str) or not s.strip():
        return {}, None, None
    counts = {}
    for token in s.split("|"):
        m = _pattern_token.fullmatch(token.strip())
        if not m:
            continue
        tissue = m.group(1).strip().lower()
        cnt = int(m.group(2))
        counts[tissue] = cnt
    if not counts:
        return {}, None, None
    max_tissue = max(counts, key=lambda k: counts[k])
    max_count = counts[max_tissue]
    # unique max?
    n_max = sum(1 for v in counts.values() if v == max_count)
    if n_max > 1:
        return counts, None, max_count
    return counts, max_tissue, max_count

def round_to_step_half_up_series(s: pd.Series, step: float = 0.005) -> pd.Series:
    """
    Round to nearest `step` using HALF-UP (not banker's), preserving NaNs.
    Example with step=0.005:
      0.0222 -> 0.020
      0.0230 -> 0.025
      0.0276 -> 0.030
    """
    s = pd.to_numeric(s, errors="coerce")
    factor = int(round(1.0 / step))  # e.g., 200 for 0.005
    return (np.floor(s * factor + 0.5) / factor)

# =========================
# Main
# =========================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--outfile_breast_top", required=True, help="Variants where breast is dominant among tissues")
    ap.add_argument("--outfile_germline", required=False, help="Variants indicated as germline with no conflict (separate output)")
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t", low_memory=False)

    # --- column handles (support _func or base names) ---
    col_sample = pick(df, "Tumor_Sample_Barcode", "Tumor_Sample_Barcode_func")
    col_gene   = pick(df, "Hugo_Symbol", "Hugo_Symbol_func")
    col_class  = pick(df, "Variant_Classification", "Variant_Classification_func")
    col_vtype  = pick(df, "Variant_Type", "Variant_Type_func")
    col_ref    = pick(df, "Reference_Allele", "Reference_Allele_func")
    col_t1     = pick(df, "Tumor_Seq_Allele1", "Tumor_Seq_Allele1_func")
    col_t2     = pick(df, "Tumor_Seq_Allele2", "Tumor_Seq_Allele2_func")
    col_chr    = "Chromosome" if "Chromosome" in df.columns else None
    col_start  = "Start_Position" if "Start_Position" in df.columns else None
    col_end    = "End_Position" if "End_Position" in df.columns else None

    col_tref   = pick(df, "t_ref_count", "t_ref_count_func")
    col_talt   = pick(df, "t_alt_count", "t_alt_count_func")
    col_tdepth = pick(df, "t_depth")          # optional
    col_vaf    = pick(df, "VAF")              # optional

    col_onco   = pick(df, "ONCOGENIC")
    col_cos_tt = pick(df, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")
    col_hot    = pick(df, "hotspot")

    # for conflict step
    col_mutsom = pick(df, "MUTATION_SOMATIC_STATUS")
    col_filter = pick(df, "FILTER")
    col_ml     = pick(df, "PureCN_ML_SOMATIC", "ML.SOMATIC", "ml.somatic")

    essentials = [col_sample, col_gene, col_class, col_tref, col_talt]
    if any(c is None for c in essentials):
        missing = [n for (n,c) in zip(["sample","gene","class","t_ref","t_alt"], essentials) if c is None]
        raise SystemExit(f"Missing required columns for second filter: {missing}")

    # --- numeric & VAF ---
    df[col_tref] = to_num(df[col_tref])
    df[col_talt] = to_num(df[col_talt])
    if col_tdepth and df[col_tdepth].notna().any():
        depth = to_num(df[col_tdepth])
    else:
        depth = df[col_tref].fillna(0) + df[col_talt].fillna(0)

    if col_vaf and df[col_vaf].notna().any():
        vaf = to_num(df[col_vaf])
    else:
        denom = (df[col_tref].fillna(0) + df[col_talt].fillna(0)).replace(0, pd.NA)
        vaf = df[col_talt] / denom

    # Use rounded VAF for thresholds (half-up to nearest 0.005)
    vaf_cmp = round_to_step_half_up_series(vaf, 0.005)

    # --- STEP 1: Global Filtering (valid variant call) ---
    valid1 = (((vaf_cmp >= 0.02) | (df[col_talt] >= 10)) & (depth >= 20))
    valid2 = (vaf_cmp >= 0.10)
    df1 = df[valid1 | valid2].copy()

    # --- STEP 1.5: Drop biologically irrelevant ---
    cls = df1[col_class].astype(str)
    bad = (
        cls.str.contains(r"\bSilent\b", case=False) |
        cls.str.contains(r"\bRNA\b", case=False) |
        cls.str.contains(r"\bIGR\b", case=False) |
        cls.str.contains(r"\bIntron\b", case=False) |
        cls.str.contains(r"5.?UTR", case=False) |
        cls.str.contains(r"3.?UTR", case=False)
    )
    df2 = df1[~bad].copy()

    # --- STEP 2: evidence keepers (LOF / OncoKB / hotspot / pathogenicity + optional breast weight) ---
    col_class = pick(df2, "Variant_Classification", "Variant_Classification_func")
    col_onco  = pick(df2, "ONCOGENIC")
    col_hot   = pick(df2, "hotspot")
    col_imp   = pick(df2, "IMPACT")
    col_sift  = pick(df2, "SIFT")
    col_poly  = pick(df2, "PolyPhen")
    col_cadd  = pick(df2, "CADD")  # optional
    col_cos_tt= pick(df2, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")

    # 1) LOF keeper
    lof_types = {"Splice_Site","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins"}
    keep_lof = df2[col_class].isin(lof_types)

    # 2) OncoKB keeper
    keep_oncokb = False
    if col_onco:
        onc = df2[col_onco].astype(str).str.strip().str.lower()
        keep_oncokb = onc.isin({"oncogenic","likely oncogenic"})

    # 3) Hotspot keeper: hotspot column uses "+" for true; restrict to missense
    keep_hotspot = False
    if col_hot:
        hot = df2[col_hot].astype(str).str.strip()
        keep_hotspot = (hot.eq("+") & df2[col_class].eq("Missense_Mutation"))

    # 4) Pathogenicity score (IMPACT + SIFT + PolyPhen + optional CADD)
    def _num_in_parens_or_self(x):
        if pd.isna(x): return pd.NA
        s = str(x)
        m = re.search(r"\(([-+]?\d*\.?\d+)\)", s)
        if m:
            try: return float(m.group(1))
            except: return pd.NA
        try: return float(s)
        except: return pd.NA

    # IMPACT: HIGH=+2, MODERATE=+1
    if col_imp:
        imp = df2[col_imp].astype(str).str.upper()
        w_impact = np.where(imp.eq("HIGH"), 2.0, np.where(imp.eq("MODERATE"), 1.0, 0.0))
    else:
        w_impact = np.zeros(len(df2), dtype=float)

    # SIFT
    if col_sift:
        sift_txt  = df2[col_sift].astype(str)
        sift_sc   = pd.to_numeric(sift_txt.map(_num_in_parens_or_self), errors="coerce")
        sift_low  = sift_txt.str.contains(r"deleterious_low_confidence", case=False, na=False)
        sift_del  = sift_txt.str.contains(r"\bdeleterious\b", case=False, na=False) & ~sift_low
        sift_strong = sift_del | (sift_sc.fillna(1.0) <= 0.05)
        w_sift = (sift_strong.astype(float) * 1.0) + (sift_low.astype(float) * 0.5)
    else:
        w_sift = np.zeros(len(df2), dtype=float)

    # PolyPhen
    if col_poly:
        poly_txt  = df2[col_poly].astype(str)
        poly_sc   = pd.to_numeric(poly_txt.map(_num_in_parens_or_self), errors="coerce")
        poly_prob = poly_txt.str.contains(r"\bprobably_damaging\b", case=False, na=False) | (poly_sc.fillna(0.0) >= 0.85)
        poly_poss = poly_txt.str.contains(r"\bpossibly_damaging\b", case=False, na=False) | poly_sc.between(0.15, 0.85, inclusive="both").fillna(False)
        w_poly = (poly_prob.astype(float) * 1.0) + ((~poly_prob & poly_poss).astype(float) * 0.5)
    else:
        w_poly = np.zeros(len(df2), dtype=float)

    # CADD (optional): >=25 => +1
    if col_cadd and df2[col_cadd].notna().any():
        cadd = pd.to_numeric(df2[col_cadd], errors="coerce")
        w_cadd = (cadd.fillna(-1) >= 25).astype(float) * 1.0
    else:
        w_cadd = np.zeros(len(df2), dtype=float)

    # COSMIC breast soft bonus
    w_breast = np.zeros(len(df2), dtype=float)
    if col_cos_tt and df2[col_cos_tt].notna().any():
        tt = df2[col_cos_tt].astype(str).str.lower()
        def _breast_hits(s):
            m = re.search(r"breast\((\d+)\)", s)
            return int(m.group(1)) if m else 0
        bh = tt.map(_breast_hits)
        w_breast = np.select([bh >= 1000, bh >= 200], [1.0, 0.5], default=0.0)

    patho_score = w_impact + w_sift + w_poly + w_cadd + w_breast
    keep_predpath = (patho_score >= 2.0)

    # Final Step-2 keepers (OR logic)
    df3 = df2[ keep_lof | keep_oncokb | keep_hotspot | keep_predpath ].copy()

    # --- STEP 4A: extra strict gate for MB5 & MB6 with rescue for OncoKB/hotspot ---
    if col_sample is not None and STRICT_SAMPLES:
        # Recompute VAF/Depth on df3; use rounded VAF for thresholding here too
        if col_vaf and df3[col_vaf].notna().any():
            vaf3 = to_num(df3[col_vaf])
        else:
            denom3 = (df3[col_tref].fillna(0) + df3[col_talt].fillna(0)).replace(0, pd.NA)
            vaf3 = df3[col_talt] / denom3
        vaf3_cmp = round_to_step_half_up_series(vaf3, 0.005)

        if col_tdepth and df3[col_tdepth].notna().any():
            depth3 = to_num(df3[col_tdepth])
        else:
            depth3 = (df3[col_tref].fillna(0) + df3[col_talt].fillna(0))

        samp = df3[col_sample].astype(str)

        # Rescue signals
        is_onco = pd.Series(False, index=df3.index)
        if col_onco:
            is_onco = df3[col_onco].astype(str).str.strip().str.lower().isin(
                {"oncogenic", "likely oncogenic"}
            )

        is_hotspot = pd.Series(False, index=df3.index)
        if col_hot:
            is_hotspot = (df3[col_hot].astype(str).str.strip().eq("+")) & df3[col_class].eq("Missense_Mutation")

        before_n = len(df3)
        ok_masks = []
        for s, (min_alt, min_depth, min_vaf) in STRICT_SAMPLES.items():
            sel = samp.eq(s)

            # Base strict gate for that sample
            base_ok = (
                (df3[col_talt] >= min_alt) &
                (depth3 >= min_depth) &
                (vaf3_cmp.fillna(0) >= min_vaf)
            )

            # Rescue: OncoKB oncogenic/likely OR hotspot (+ & missense) → bypass gate
            rescue_ok = (is_onco | is_hotspot)

            # Keep if not that sample, or if base_ok, or if rescued
            ok_masks.append(~sel | base_ok | rescue_ok)

        if ok_masks:
            df3 = df3[np.logical_and.reduce(ok_masks)]

        print(f"[strict+rescue] filtered {before_n - len(df3)} rows with MB5/MB6 strict rule (OncoKB/hotspot rescue)")

    # --- STEP 4b: conflict-based pruning (applies if both FILTER and ML columns exist) ---
    if (col_filter is not None) and (col_ml is not None):
        s_filter = df3[col_filter].astype(str).str.strip().str.upper()
        filter_is_pass = s_filter.isin({"PASS", ".", ""})
        ml_bool = coerce_bool(df3[col_ml])

        # Recompute VAF and depth for df3 (raw values are fine here)
        if col_vaf and df3[col_vaf].notna().any():
            vaf3b = to_num(df3[col_vaf])
        else:
            denom3b = (df3[col_tref].fillna(0) + df3[col_talt].fillna(0)).replace(0, pd.NA)
            vaf3b = df3[col_talt] / denom3b

        if col_tdepth and df3[col_tdepth].notna().any():
            depth3b = to_num(df3[col_tdepth])
        else:
            depth3b = (df3[col_tref].fillna(0) + df3[col_talt].fillna(0))

        # define "conflict"
        conflict = (
            (ml_bool.fillna(False) & ~filter_is_pass) |
            ((ml_bool == False).fillna(False) & filter_is_pass)
        )
        conflict = conflict.fillna(False)

        # Somatic evidence (NO COSMIC-any here)
        som_onc = False
        if col_onco:
            onc2 = df3[col_onco].astype(str).str.lower().str.strip()
            som_onc = onc2.isin({"oncogenic", "likely oncogenic"})
        som_hot = False
        if col_hot:
            hot2 = df3[col_hot].astype(str).str.strip()
            som_hot = (hot2 == "+") & df3[col_class].eq("Missense_Mutation")
        som_status = False
        if col_mutsom:
            ms = df3[col_mutsom].astype(str).str.strip().str.lower()
            som_status = ms.str.contains(r"\bconfirmed\s+somatic\s+variant\b")

        # PureCN posterior
        col_post = "PureCN_POSTERIOR_SOMATIC" if "PureCN_POSTERIOR_SOMATIC" in df3.columns else None
        posterior = to_num(df3[col_post]) if col_post else pd.Series(pd.NA, index=df3.index)
        post_somatic = posterior.fillna(0) >= 0.90
        post_strong  = posterior.fillna(0) >= 0.98

        # COSMIC breast evidence — ONLY for rescue (no "any breast" in is_somatic)
        def _breast_strong_func(s):
            counts, top_tissue, maxcnt = parse_cosmic_tissues(s)
            if not counts:
                return False
            breast = counts.get("breast", 0)
            if breast >= 1000:
                return True
            if top_tissue is None:
                return (breast == maxcnt) and (maxcnt > 0)
            return (top_tissue == "breast") or (maxcnt > 0 and breast >= 0.8 * maxcnt)

        som_breast_strong = pd.Series(False, index=df3.index)
        col_cos_tt_conf = pick(df3, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")
        if col_cos_tt_conf:
            tt = df3[col_cos_tt_conf].astype(str).str.lower()
            som_breast_strong = tt.apply(_breast_strong_func)

        # --- combine somatic evidence (removed COSMIC-any) ---
        truly_somatic = (som_onc | som_hot | som_status | post_somatic)

        # Support logic
        base_support = (vaf3b.fillna(0) >= 0.15)
        post_rescue  = post_strong
        depth3b_safe = depth3b.fillna(0)
        breast_rescue = (som_breast_strong.astype(bool) &
                         ((vaf3b.fillna(0) >= 0.10) | (depth3b_safe >= 40)))
        support_ok = base_support | post_rescue | breast_rescue

        # New: strong-support override clears conflicts alone
        strong_support_override = (
            (vaf3b.fillna(0) >= CONFLICT_STRONG_VAF) &
            (depth3b_safe >= CONFLICT_STRONG_DEPTH)
        )

        # Identify "germline, no-conflict" rows to exclude from main and (optionally) export
        ml_is_false     = (ml_bool == False).fillna(False)
        not_pass_filter = ~filter_is_pass
        germline_no_conflict = (ml_is_false & not_pass_filter)
        if col_mutsom:
            ms2 = df3[col_mutsom].astype(str).str.strip().str.lower()
            germline_no_conflict |= ms2.str.contains(r"\bgermline\b", na=False)

        df_germline_nc = df3[germline_no_conflict].copy()

        keep_conflict = (
            (~conflict) |
            (conflict & ((support_ok & truly_somatic) | post_rescue | strong_support_override))
        )
        keep_conflict = keep_conflict & (~germline_no_conflict)

        # --- Conflict reporting (per-sample) ---
        if col_sample is not None:
            samp2 = df3[col_sample].astype(str)
            conflict_idx = conflict[conflict].index
            kept_idx     = keep_conflict[keep_conflict].index
            kept_conflict_idx   = conflict_idx.intersection(kept_idx)
            dropped_conflict_idx= conflict_idx.difference(kept_idx)

            total_conflicts = len(conflict_idx)
            kept_conflicts  = len(kept_conflict_idx)
            dropped_conflicts = len(dropped_conflict_idx)

            print(f"[conflict] total conflicts: {total_conflicts}  kept: {kept_conflicts}  dropped: {dropped_conflicts}")

            if total_conflicts > 0:
                per_sample_conf = samp2.loc[conflict_idx].value_counts()
                print("[conflict] conflicts per sample:")
                for s, n in per_sample_conf.items():
                    print(f"  {s}: {n}")

                if kept_conflicts > 0:
                    per_sample_kept = samp2.loc[kept_conflict_idx].value_counts()
                    print("[conflict] kept conflicts per sample:")
                    for s, n in per_sample_kept.items():
                        print(f"  {s}: {n}")

                if dropped_conflicts > 0:
                    per_sample_dropped = samp2.loc[dropped_conflict_idx].value_counts()
                    print("[conflict] dropped conflicts per sample:")
                    for s, n in per_sample_dropped.items():
                        print(f"  {s}: {n}")

        before_n = len(df3)
        df3 = df3[keep_conflict].copy()
        print(f"[conflict] dropped {before_n - len(df3)} rows after strict conflict rule (post germline exclusion)")

    # --- STEP 5: de-duplicate on available keys ---
    key_cols = [c for c in [
        col_sample, col_gene, col_chr, col_start, col_end, col_class, col_vtype, col_ref, col_t1, col_t2
    ] if c is not None and c in df3.columns]
    if key_cols:
        before_n = len(df3)
        df3 = df3.drop_duplicates(subset=key_cols)
        print(f"[dedup] removed {before_n - len(df3)} duplicate rows (subset={','.join(key_cols)})")

    # --- SAVE main output ---
    df3.to_csv(args.outfile, sep="\t", index=False)

    # Optional: write germline-no-conflict variants to their own file
    try:
        if (col_filter is not None) and (col_ml is not None) and ('df_germline_nc' in locals()) and (args.outfile_germline):
            if 'key_cols' in locals() and key_cols:
                df_germline_nc = df_germline_nc.drop_duplicates(subset=key_cols)
            df_germline_nc.to_csv(args.outfile_germline, sep="\t", index=False)
            print(f"[germline] wrote {len(df_germline_nc)} germline (no-conflict) rows → {args.outfile_germline}")
    except Exception as e:
        print(f"[germline] failed to write outfile_germline: {e}")

    # =========================
    # Extra breast-focused outputs
    # =========================
    if col_cos_tt is None:
        print("[breast] COSMIC_tissue_types_affected/TUMOR_TYPE_SUMMARY column not found — skipping breast outputs.")
    else:
        cosmic_series = df3[col_cos_tt].astype(str)

        # breast_top: dominant (≥70% of max tissue count)
        def is_breast_dominant(s, frac=0.7):
            counts, _, _ = parse_cosmic_tissues(s)
            if not counts or "breast" not in counts:
                return False
            breast = counts["breast"]
            maxval = max(counts.values())
            return breast >= frac * maxval

        mask_breast_top = cosmic_series.apply(is_breast_dominant)
        df_breast_top = df3[mask_breast_top].copy()
        df_breast_top.to_csv(args.outfile_breast_top, sep="\t", index=False)
        print(f"[breast] breast_top kept {len(df_breast_top)} rows → {args.outfile_breast_top}")

    # --- brief summary ---
    print(f"[summary] start={len(df)}  step2={len(df2)}  step3={len(df3)}")
    print(f"wrote main → {args.outfile}")

if __name__ == "__main__":
    main()
