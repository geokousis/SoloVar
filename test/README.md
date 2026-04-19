# SoloVar CDH1 Smoke Test

Tests the MAF post-processing pipeline on 3 synthetic CDH1 variants.  
No BAMs, no mapping, no pysam — runs in under 10 seconds.

## What it tests

| Step | Script | What it checks |
|------|--------|----------------|
| 1 | `slim_maf.py` | Column coalescing, VAF computation, COSMIC ID extraction |
| 2 | `filter_basic.py` | gnomADe_AF gate (≤0.001), FILTER/ML conflict annotation |
| 3 | `filter_biological.py` | OncoKB / hotspot / pathogenicity-based biological filter |

## Run

```bash
cd SoloVar/test/
bash run_test.sh
```

## Expected output

```
Input variants:                3
After slim_maf:                3
After filter_basic:            2
After filter_biological:       2
Breast-dominant:               2
```

- **Variant 1** (chr16:68777022 C>T): Oncogenic, LEVEL_1, hotspot (+), HIGH IMPACT, confirmed somatic → **kept**
- **Variant 2** (chr16:68773838 G>A): Likely Oncogenic, MODERATE IMPACT, deleterious SIFT → **kept**
- **Variant 3** (chr16:68775640 T>A): gnomADe_AF=0.002 (common variant), byFrequency dbSNP flag → **dropped at filter_basic**

## Test data files

| File | Description |
|------|-------------|
| `data/cdh1_region.fa` | 10 kb CDH1 reference sequence (chr16:68771576-68781575, hg38) |
| `data/cdh1_test.vcf` | 3 raw somatic VCF records for CDH1 |
| `data/cdh1_test.maf` | Pre-annotated MAF (simulates vcf_annotation_pipeline output) |
| `data/cdh1.bed` | CDH1 capture region BED (for use with other pipeline stages) |

## What is NOT tested here

- `filter_af.py` — requires `pysam` and bgzipped+indexed VCFs (run via `vcf_annotation_pipeline.sh`)
- `merge_purecn.py` — requires PureCN output directory
- `top_mutated_genes.py` — requires a multi-sample MAF
- CNV scripts — require CNVkit `.call.cns` files
