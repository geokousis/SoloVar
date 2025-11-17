#!/usr/bin/env python3
"""
Two-stage variant filtering pipeline.

Note: thresholds / parameters kept identical to your original script.
This rewrite only improves structure + readability.
"""

import argparse
import re

import numpy as np
import pandas as pd


# =========================
# Config
# =========================

# Extra strict gate for noisy samples; rescue (OncoKB/hotspot) bypasses this.
STRICT_SAMPLES = {
    "MB5": (10, 40, 0.05),  # (min_alt, min_depth, min_vaf)
    "MB6": (10, 40, 0.05),
}

# Conflict override: allow strong base support alone to clear conflicts
CONFLICT_STRONG_VAF = 0.30
CONFLICT_STRONG_DEPTH = 50


# =========================
# Helpers
# =========================

def pick(df: pd.DataFrame, *candidates):
    """Return the first existing column name from candidates, or None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def to_num(s):
    """Convenience wrapper for numeric conversion with NaNs."""
    return pd.to_numeric(s, errors="coerce")


def coerce_bool(series: pd.Series) -> pd.Series:
    """
    Coerce a column with strings / numbers to nullable boolean (True/False/NA).
    """
    s = series.astype(str).str.strip().str.lower()
    true_vals = {"true", "t", "1", "yes", "y"}
    false_vals = {"false", "f", "0", "no", "n"}

    out = pd.Series(pd.NA, index=series.index, dtype="boolean")
    out[s.isin(true_vals)] = True
    out[s.isin(false_vals)] = False
    return out


# parse strings like "biliary_tract(2)|breast(11)|kidney(92)"
_pattern_token = re.compile(r"\s*([^()|]+)\((\d+)\)\s*")


def parse_cosmic_tissues(s: str):
    """
    Parse COSMIC_tissue_types_affected / TUMOR_TYPE_SUMMARY style strings.

    Returns
    -------
    counts : dict[tissue -> count]
    max_tissue : str | None      # tissue with highest count (if unique)
    max_count  : int | None
    """
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


def compute_vaf_and_depth(
    df: pd.DataFrame,
    col_tref: str,
    col_talt: str,
    col_tdepth: str | None,
    col_vaf: str | None,
    round_step: float | None = None,
):
    """
    Compute VAF and depth using the same logic you used everywhere in the
    original script. Optionally returns a rounded VAF for thresholding.
    """
    # ensure numeric
    df[col_tref] = to_num(df[col_tref])
    df[col_talt] = to_num(df[col_talt])

    # depth
    if col_tdepth and df[col_tdepth].notna().any():
        depth = to_num(df[col_tdepth])
    else:
        depth = df[col_tref].fillna(0) + df[col_talt].fillna(0)

    # VAF
    if col_vaf and df[col_vaf].notna().any():
        vaf = to_num(df[col_vaf])
    else:
        denom = (df[col_tref].fillna(0) + df[col_talt].fillna(0)).replace(0, pd.NA)
        vaf = df[col_talt] / denom

    if round_step is not None:
        vaf_cmp = round_to_step_half_up_series(vaf, round_step)
    else:
        vaf_cmp = None

    return vaf, depth, vaf_cmp


def compute_pathogenicity_score(
    df2: pd.DataFrame,
    col_imp: str | None,
    col_sift: str | None,
    col_poly: str | None,
    col_cadd: str | None,
    col_cos_tt: str | None,
) -> pd.Series:
    """Exact copy of the original patho-score logic."""

    def _num_in_parens_or_self(x):
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

    # IMPACT: HIGH=+2, MODERATE=+1
    if col_imp:
        imp = df2[col_imp].astype(str).str.upper()
        w_impact = np.where(
            imp.eq("HIGH"),
            2.0,
            np.where(imp.eq("MODERATE"), 1.0, 0.0),
        )
    else:
        w_impact = np.zeros(len(df2), dtype=float)

    # SIFT
    if col_sift:
        sift_txt = df2[col_sift].astype(str)
        sift_sc = pd.to_numeric(
            sift_txt.map(_num_in_parens_or_self),
            errors="coerce",
        )
        # only *deleterious_low_confidence* gets +0.5
        sift_low = sift_txt.str.contains(
            r"deleterious_low_confidence",
            case=False,
            na=False,
        )
        sift_del = sift_txt.str.contains(
            r"\bdeleterious\b",
            case=False,
            na=False,
        ) & ~sift_low
        sift_strong = sift_del | (sift_sc.fillna(1.0) <= 0.05)
        w_sift = (sift_strong.astype(float) * 1.0) + (
            sift_low.astype(float) * 0.5
        )
    else:
        w_sift = np.zeros(len(df2), dtype=float)

    # PolyPhen
    if col_poly:
        poly_txt = df2[col_poly].astype(str)
        poly_sc = pd.to_numeric(
            poly_txt.map(_num_in_parens_or_self),
            errors="coerce",
        )
        poly_prob = poly_txt.str.contains(
            r"\bprobably_damaging\b",
            case=False,
            na=False,
        ) | (poly_sc.fillna(0.0) >= 0.85)
        poly_poss = poly_txt.str.contains(
            r"\bpossibly_damaging\b",
            case=False,
            na=False,
        ) | poly_sc.between(0.15, 0.85, inclusive="both").fillna(False)

        w_poly = (
            poly_prob.astype(float) * 1.0
            + ((~poly_prob & poly_poss).astype(float) * 0.5)
        )
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
        w_breast = np.select(
            [bh >= 1000, bh >= 200],
            [1.0, 0.5],
            default=0.0,
        )

    patho_score = (
        w_impact.astype(float)
        + w_sift.astype(float)
        + w_poly.astype(float)
        + w_cadd.astype(float)
        + w_breast.astype(float)
    )
    return pd.Series(patho_score, index=df2.index)


def resolve_conflicts(
    df3: pd.DataFrame,
    col_filter: str,
    col_ml: str,
    col_mutsom: str | None,
    col_post_som: str | None,
    col_tref: str,
    col_talt: str,
    col_tdepth: str | None,
    col_vaf: str | None,
    col_onco: str | None,
    col_hot: str | None,
):
    """
    Conflict resolution between FILTER and ML somatic, germline removal,
    somatic rescue, etc. Logic identical to original script, just grouped here.
    """
    # basic filter / ML
    filter_str = df3[col_filter].astype(str).str.upper().str.strip()
    filter_is_pass = filter_str.isin({"PASS", ".", ""})

    ml_bool = coerce_bool(df3[col_ml])

    # recompute raw VAF/depth for df3
    if col_vaf and df3[col_vaf].notna().any():
        vaf3b = to_num(df3[col_vaf])
    else:
        denom3b = (df3[col_tref].fillna(0) + df3[col_talt].fillna(0)).replace(0, pd.NA)
        vaf3b = df3[col_talt] / denom3b

    if col_tdepth and df3[col_tdepth].notna().any():
        depth3b = to_num(df3[col_tdepth])
    else:
        depth3b = df3[col_tref].fillna(0) + df3[col_talt].fillna(0)

    # define conflict
    conflict = (
        (ml_bool.fillna(False) & ~filter_is_pass)
        | ((ml_bool == False).fillna(False) & filter_is_pass)
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
        col_class = pick(df3, "Variant_Classification", "Variant_Classification_func")
        if col_class:
            cls2 = df3[col_class].astype(str)
            som_hot = (hot2 == "+") & cls2.eq("Missense_Mutation")

    som_status = False
    if col_mutsom:
        ms2 = df3[col_mutsom].astype(str).str.lower()
        som_status = ms2.str.contains("confirmed somatic variant", na=False)

    post_somatic = False
    post_strong = False
    if col_post_som and col_post_som in df3.columns:
        post = pd.to_numeric(df3[col_post_som], errors="coerce")
        post_somatic = (post >= 0.90).fillna(False)
        post_strong = (post >= 0.98).fillna(False)

    # breast-strong based on COSMIC (separate from generic "any COSMIC")
    som_breast_strong = pd.Series(False, index=df3.index)
    col_cos_tt_conf = pick(df3, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")
    if col_cos_tt_conf:
        tt = df3[col_cos_tt_conf].astype(str).str.lower()

        def _breast_strong_func(s):
            counts, top_tissue, maxcnt = parse_cosmic_tissues(s)
            if not counts:
                return False
            # If there's an unambiguous top tissue and it's breast
            if top_tissue == "breast" and maxcnt:
                return True
            # otherwise consider breast "strong" if it reaches 80% of max
            breast = counts.get("breast", 0)
            if not maxcnt:
                return False
            return (top_tissue == "breast") or (
                maxcnt > 0 and breast >= 0.8 * maxcnt
            )

        som_breast_strong = tt.apply(_breast_strong_func)

    # combine somatic evidence (removed COSMIC-any)
    truly_somatic = (som_onc | som_hot | som_status | post_somatic)

    # support logic
    base_support = (vaf3b.fillna(0) >= 0.15)
    post_rescue = post_strong
    depth3b_safe = depth3b.fillna(0)

    breast_rescue = (
        som_breast_strong.astype(bool)
        & ((vaf3b.fillna(0) >= 0.10) | (depth3b_safe >= 40))
    )

    support_ok = base_support | post_rescue | breast_rescue

    # strong-support override clears conflicts alone
    strong_support_override = (
        (vaf3b.fillna(0) >= CONFLICT_STRONG_VAF)
        & (depth3b_safe >= CONFLICT_STRONG_DEPTH)
    )

    # Identify "germline, no-conflict" rows to exclude from main and (optionally) export
    ml_is_false = (ml_bool == False).fillna(False)
    not_pass_filter = ~filter_is_pass
    germline_no_conflict = (ml_is_false & not_pass_filter)

    if col_mutsom:
        ms = df3[col_mutsom].astype(str).str.lower()
        germline_text = ms.str.contains("germline", na=False)
        germline_no_conflict = germline_no_conflict | germline_text

    germline_no_conflict = germline_no_conflict.fillna(False)

    # final keep mask
    keep_conflict = (~conflict)
    keep_conflict |= (
        conflict
        & (
            (support_ok & truly_somatic)
            | post_rescue
            | strong_support_override
        )
    )
    keep_conflict &= ~germline_no_conflict
    keep_conflict = keep_conflict.fillna(False)

    return keep_conflict, germline_no_conflict


# =========================
# Main
# =========================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument(
        "--outfile_breast_top",
        required=True,
        help="Variants where breast is dominant among tissues",
    )
    ap.add_argument(
        "--outfile_germline",
        required=False,
        help="Variants indicated as germline with no conflict (separate output)",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.infile, sep="\t", low_memory=False)

    # --- column handles (support _func or base names) ---
    col_sample = pick(df, "Tumor_Sample_Barcode", "Tumor_Sample_Barcode_func")
    col_gene = pick(df, "Hugo_Symbol", "Hugo_Symbol_func")
    col_class = pick(df, "Variant_Classification", "Variant_Classification_func")
    col_vtype = pick(df, "Variant_Type")
    col_ref = pick(df, "Reference_Allele")

    col_t1 = pick(df, "Tumor_Seq_Allele1", "Tumor_Seq_Allele1_func")
    col_t2 = pick(df, "Tumor_Seq_Allele2", "Tumor_Seq_Allele2_func")

    col_chr = pick(df, "Chromosome")
    col_start = pick(df, "Start_Position")
    col_end = pick(df, "End_Position")

    col_tref = pick(df, "t_ref_count", "t_ref_count_func")
    col_talt = pick(df, "t_alt_count", "t_alt_count_func")
    col_tdepth = pick(df, "t_depth", "t_depth_func")
    col_vaf = pick(df, "VAF", "VAF_func")

    col_onco = pick(df, "ONCOGENIC")
    col_cos_tt = pick(df, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")
    col_hot = pick(df, "hotspot")

    # for conflict step
    col_mutsom = pick(df, "MUTATION_SOMATIC_STATUS")
    col_filter = pick(df, "FILTER")
    col_ml = pick(df, "PureCN_ML_SOMATIC", "ML.SOMATIC", "ml.somatic")
    col_post_som = pick(df, "PureCN_POSTERIOR_SOMATIC", "POSTERIOR.SOMATIC")

    essentials = [col_sample, col_gene, col_class, col_tref, col_talt]
    if any(c is None for c in essentials):
        missing = [
            n
            for (n, c) in zip(
                ["sample", "gene", "class", "t_ref", "t_alt"],
                essentials,
            )
            if c is None
        ]
        raise SystemExit(
            f"Missing required columns for second filter: {missing}"
        )

    # =========================
    # STEP 1 – Global QC
    # =========================
    vaf, depth, vaf_cmp = compute_vaf_and_depth(
        df,
        col_tref=col_tref,
        col_talt=col_talt,
        col_tdepth=col_tdepth,
        col_vaf=col_vaf,
        round_step=0.005,
    )

    # original QC logic:
    valid1 = (((vaf_cmp >= 0.02) | (df[col_talt] >= 10)) & (depth >= 20))
    valid2 = (vaf_cmp >= 0.10)
    df1 = df[valid1 | valid2].copy()
    print(f"[step1] kept {len(df1)}/{len(df)} rows after basic QC")

    # =========================
    # STEP 2 – Drop uninteresting classes + evidence keepers
    # (2a: bio-irrelevant removal, 2b: LOF/onco/hotspot/patho)
    # =========================
    # 2a) remove silent / UTR / intron / RNA etc.
    cls = df1[col_class].astype(str)
    bad = (
        cls.str.contains(r"\bSilent\b", case=False)
        | cls.str.contains(r"\bRNA\b", case=False)
        | cls.str.contains(r"\bIGR\b", case=False)
        | cls.str.contains(r"\bIntron\b", case=False)
        | cls.str.contains(r"5.?UTR", case=False)
        | cls.str.contains(r"3.?UTR", case=False)
    )
    df2 = df1[~bad].copy()

    # 2b) evidence keepers
    col_class2 = pick(df2, "Variant_Classification", "Variant_Classification_func")
    col_onco2 = pick(df2, "ONCOGENIC")
    col_hot2 = pick(df2, "hotspot")
    col_imp = pick(df2, "IMPACT")
    col_sift = pick(df2, "SIFT")
    col_poly = pick(df2, "PolyPhen")
    col_cadd = pick(df2, "CADD")
    col_cos_tt2 = pick(df2, "COSMIC_tissue_types_affected", "TUMOR_TYPE_SUMMARY")

    # 1) LOF keeper
    lof_types = {
        "Splice_Site",
        "Nonsense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
    }
    keep_lof = df2[col_class2].isin(lof_types)

    # 2) OncoKB keeper
    keep_oncokb = False
    if col_onco2:
        onc = df2[col_onco2].astype(str).str.lower().str.strip()
        keep_oncokb = onc.isin({"oncogenic", "likely oncogenic"})

    # 3) hotspot missense
    keep_hotspot = False
    if col_hot2 and col_class2:
        hot = df2[col_hot2].astype(str).str.strip()
        cls2 = df2[col_class2].astype(str)
        keep_hotspot = (hot == "+") & cls2.eq("Missense_Mutation")

    # 4) patho score
    patho_score = compute_pathogenicity_score(
        df2,
        col_imp=col_imp,
        col_sift=col_sift,
        col_poly=col_poly,
        col_cadd=col_cadd,
        col_cos_tt=col_cos_tt2,
    )
    keep_predpath = (patho_score >= 2.0)

    # Final Step-2 keepers (OR logic)
    df3 = df2[keep_lof | keep_oncokb | keep_hotspot | keep_predpath].copy()
    print(f"[step2] evidence keepers: {len(df3)}/{len(df2)}")

    # =========================
    # STEP 3 – Extra strict gate for problematic samples (high contamination)
    # =========================
    if col_sample is not None and STRICT_SAMPLES:
        # Recompute VAF/Depth on df3; rounded VAF for thresholding here too
        vaf3, depth3, vaf3_cmp = compute_vaf_and_depth(
            df3,
            col_tref=col_tref,
            col_talt=col_talt,
            col_tdepth=col_tdepth,
            col_vaf=col_vaf,
            round_step=0.005,
        )

        samp = df3[col_sample].astype(str)

        # Rescue signals
        is_onco = pd.Series(False, index=df3.index)
        if col_onco:
            onc3 = df3[col_onco].astype(str).str.lower().str.strip()
            is_onco = onc3.isin({"oncogenic", "likely oncogenic"})

        is_hot = pd.Series(False, index=df3.index)
        if col_hot and col_class:
            hot3 = df3[col_hot].astype(str).str.strip()
            cls3 = df3[col_class].astype(str)
            is_hot = (hot3 == "+") & cls3.eq("Missense_Mutation")

        rescue = (is_onco | is_hot)

        keep_mask = pd.Series(True, index=df3.index)
        for sample_name, (min_alt, min_depth, min_vaf) in STRICT_SAMPLES.items():
            mask_samp = samp.eq(sample_name)
            if not mask_samp.any():
                continue

            strict_keep = (
                (df3[col_talt] >= min_alt)
                & (depth3 >= min_depth)
                & (vaf3_cmp >= min_vaf)
            )
            mask_rescue = rescue & mask_samp
            keep_samp = mask_samp & (strict_keep | mask_rescue)
            keep_mask &= (~mask_samp | keep_samp)

        before = len(df3)
        df3 = df3[keep_mask].copy()
        after = len(df3)
        print(f"[step3] STRICT_SAMPLES removed {before - after} rows (now {after})")

    # =========================
    # STEP 4 – Conflict resolution & germline handling
    # =========================
    df_germline_nc = None
    if col_filter is not None and col_ml is not None:
        keep_conflict, germline_no_conflict = resolve_conflicts(
            df3=df3,
            col_filter=col_filter,
            col_ml=col_ml,
            col_mutsom=col_mutsom,
            col_post_som=col_post_som,
            col_tref=col_tref,
            col_talt=col_talt,
            col_tdepth=col_tdepth,
            col_vaf=col_vaf,
            col_onco=col_onco,
            col_hot=col_hot,
        )

        # export germline no-conflict separately if requested
        df_germline_nc = df3[germline_no_conflict].copy()
        df3 = df3[keep_conflict].copy()
        print(
            f"[step4] resolved conflicts: kept {len(df3)}, "
            f"germline_no_conflict={len(df_germline_nc)}"
        )

        if args.outfile_germline:
            df_germline_nc.to_csv(
                args.outfile_germline,
                sep="\t",
                index=False,
            )
            print(
                f"[germline] wrote {len(df_germline_nc)} germline/no-conflict "
                f"variants → {args.outfile_germline}"
            )

    # =========================
    # STEP 5 – De-duplicate
    # =========================
    key_cols = []
    for c in [
        col_sample,
        col_gene,
        col_chr,
        col_start,
        col_end,
        col_class,
        col_vtype,
        col_ref,
        col_t1,
        col_t2,
    ]:
        if c:
            key_cols.append(c)

    if key_cols:
        before = len(df3)
        df3 = df3.drop_duplicates(subset=key_cols)
        print(
            f"[step5] de-duplicated on {len(key_cols)} cols: "
            f"{before} → {len(df3)}"
        )

    # --- write main ---
    df3.to_csv(args.outfile, sep="\t", index=False)
    print(f"[main] wrote {len(df3)} rows → {args.outfile}")

    # =========================
    # Extra: breast-dominant COSMIC subset
    # =========================
    if col_cos_tt and col_cos_tt in df3.columns:
        cosmic_series = df3[col_cos_tt].astype(str)

        def is_breast_dominant(s: str, frac: float = 0.7) -> bool:
            counts, top_tissue, maxcnt = parse_cosmic_tissues(s)
            if not counts or not maxcnt:
                return False
            breast = counts.get("breast", 0)
            maxval = max(counts.values())
            return breast >= frac * maxval

        mask_breast_top = cosmic_series.apply(is_breast_dominant)
        df_breast_top = df3[mask_breast_top].copy()
        df_breast_top.to_csv(
            args.outfile_breast_top,
            sep="\t",
            index=False,
        )
        print(
            f"[breast] breast_top kept {len(df_breast_top)} rows → "
            f"{args.outfile_breast_top}"
        )

    # --- brief summary ---
    print(f"[summary] start={len(df)}  after_step2={len(df2)}  final={len(df3)}")
    print(f"wrote main → {args.outfile}")


if __name__ == "__main__":
    main()
