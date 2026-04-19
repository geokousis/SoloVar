import pandas as pd
import numpy as np
import re

# === load ===
df = pd.read_csv("1F_slim_pure.tsv", sep="\t", dtype=str, low_memory=False)
df = df.replace(r'^\s*$|^\.$', pd.NA, regex=True)

# coerce numerics
for c in ["t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","gnomADe_AF",
          "PureCN_POSTERIOR_SOMATIC","VAF","COSMIC_total_alterations_in_gene"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# compute VAF if missing
if "VAF" not in df.columns or df["VAF"].isna().all():
    denom = (df.get("t_ref_count",0).fillna(0) + df.get("t_alt_count",0).fillna(0)).replace(0, np.nan)
    df["VAF"] = df.get("t_alt_count",0) / denom

# helpers
def has_text(col, pattern):
    if col not in df.columns: return pd.Series(False, index=df.index)
    s = df[col].astype(str).str.lower()
    return s.str.contains(pattern, na=False)

def equals_any(col, values):
    if col not in df.columns: return pd.Series(False, index=df.index)
    s = df[col].astype(str).str.strip().str.upper()
    return s.isin(set(v.upper() for v in values))

def is_true_like(col):
    if col not in df.columns: return pd.Series(False, index=df.index)
    s = df[col].astype(str).str.strip().str.lower()
    return s.isin({"true","t","1","yes","y"})

# --- 1) QC gate ---
is_indel = df["Variant_Type"].eq("INS") | df["Variant_Type"].eq("DEL")

support_snp = (((df["VAF"] >= 0.05) | (df["t_alt_count"] >= 8)) & (df["t_depth"] >= 40)) | (df["VAF"] >= 0.10)
support_indel = (((df["VAF"] >= 0.03) | (df["t_alt_count"] >= 5)) & (df["t_depth"] >= 30)) | (df["VAF"] >= 0.08)
support_ok = np.where(is_indel.fillna(False), support_indel, support_snp)

filter_pass = equals_any("FILTER", {"PASS", ".", ""})
posterior_high = (df.get("PureCN_POSTERIOR_SOMATIC", pd.Series(0, index=df.index)) >= 0.90)
qc_ok = (support_ok) & (filter_pass | posterior_high)

# matched-normal leak
if "n_alt_count" in df.columns and "n_ref_count" in df.columns:
    n_denom = (df["n_alt_count"].fillna(0) + df["n_ref_count"].fillna(0)).replace(0, np.nan)
    n_vaf = df["n_alt_count"] / n_denom
    normal_leak = (n_vaf.fillna(0) >= 0.02)
else:
    normal_leak = pd.Series(False, index=df.index)

# --- 2) Population filter ---
pop_common = (df.get("gnomADe_AF", pd.Series(0, index=df.index)).fillna(0) >= 0.01)

# exemptions for pop filter
oncokb_present = df.get("HIGHEST_LEVEL", pd.Series("", index=df.index)).notna()
hotspot_plus = df.get("hotspot", pd.Series("", index=df.index)).astype(str).str.strip().eq("+") & df["Variant_Classification"].eq("Missense_Mutation")
som_status_confirmed = has_text("MUTATION_SOMATIC_STATUS", r"\bconfirmed\s+somatic\s+variant\b")
pop_exempt = oncokb_present | hotspot_plus | som_status_confirmed

# --- 3) Conflict rule ---
ml_true  = is_true_like("PureCN_ML_SOMATIC")
conflict = (filter_pass & (~ml_true.fillna(False))) | ((~filter_pass) & ml_true.fillna(False))

som_evidence = ml_true | posterior_high | som_status_confirmed

conflict_keep = (~conflict) | (conflict & ((df["VAF"] >= 0.10) | (df["t_depth"] >= 50) | posterior_high) & som_evidence)

# --- 4) Driver scoring ---
lof_types = {"Splice_Site","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Frame_Shift_Ins","Frame_Shift"}
is_lof = df["Variant_Classification"].isin(lof_types)

impact = df.get("IMPACT", pd.Series("", index=df.index)).astype(str).str.upper()
impact_score = np.select([impact.eq("HIGH"), impact.eq("MODERATE")], [2,1], default=0)

siftpoly = df.get("SIFTPolyPhen", pd.Series("", index=df.index)).astype(str).str.lower()
pred_score = (siftpoly.str.contains("deleterious", na=False) | siftpoly.str.contains("probably_damaging", na=False)).astype(int)

oncogenic = has_text("ONCOGENIC", r"^oncogenic$|^likely\s+oncogenic$")
oncokb_score = oncokb_present.astype(int) * 3
hotspot_score = (hotspot_plus.astype(int) * 3)
lof_score = (is_lof.astype(int) * 3)

me_text = (df.get("MUTATION_EFFECT", pd.Series("", index=df.index)).astype(str).str.lower() + " " +
           df.get("MUTATION_EFFECT_DESCRIPTION", pd.Series("", index=df.index)).astype(str).str.lower())
effect_score = me_text.str.contains(r"loss of function|gain of function|activating|inactivating", na=False).astype(int)

cgc_present = df.get("CGC_Mutation_Type", pd.Series("", index=df.index)).notna().astype(int)
cgc_breast = has_text("CGC_Tumor_Types_Somatic", r"breast").astype(int)

cosmic_gene_alt = (df.get("COSMIC_total_alterations_in_gene", pd.Series(0, index=df.index)).fillna(0) >= 100).astype(int)

driver_score = (
    oncokb_score + hotspot_score + lof_score + impact_score + pred_score + effect_score + cgc_present + cgc_breast + cosmic_gene_alt
)

driver_tier = np.select(
    [
        (driver_score >= 6) | (oncokb_present & (hotspot_plus | is_lof)),
        (driver_score >= 4)
    ],
    ["Driver (HC)","Likely driver"],
    default="Passenger/Unknown"
)

# --- 5) Final keep & reasons ---
keep = qc_ok & (~normal_leak | (oncokb_present | hotspot_plus | posterior_high)) & (~(pop_common & ~pop_exempt)) & conflict_keep

reasons = []
for i in df.index:
    r = []
    if not qc_ok[i]: r.append("QC_fail")
    if normal_leak[i] and not (oncokb_present[i] or hotspot_plus[i] or posterior_high[i]): r.append("Normal_leak")
    if pop_common[i] and not pop_exempt[i]: r.append("Common_pop")
    if not conflict_keep[i]: r.append("Conflict_drop")
    reasons.append(",".join(r) if r else "OK")
df["Filter_reason"] = reasons
df["Driver_score"] = driver_score
df["Driver_tier"]  = driver_tier
df["Keep"] = keep

# --- outputs ---
df_keep = df[keep].copy()
df_drop = df[~keep].copy()

df_keep.to_csv("filtered_keep.tsv", sep="\t", index=False)
df_drop.to_csv("filtered_drop.tsv", sep="\t", index=False)
print(f"kept {len(df_keep)}/{len(df)} rows → filtered_keep.tsv")
print("top reasons (drops):", pd.Series([r for r in reasons if r!="OK"]).value_counts().head(10).to_dict())
