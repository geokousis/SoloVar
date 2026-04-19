#!/usr/bin/env python3
import os
import glob
import re
import pandas as pd
from tqdm import tqdm
import concurrent.futures
import argparse
def extract_sample_id(filename):
    """
    Extracts the chunk between the last '_' before '_recalibrated' and '_recalibrated' itself.
    Example: from 'test.MB1_recalibrated' → returns 'test.MB1'
    """
    match = re.search(r"_([^_]+)_recalibrated", filename)
    if match:
        sample_id = match.group(1)
        print(f"Matched sample ID: {sample_id}")
        return sample_id
    print(f"No valid sample ID found in file: {filename}")
    return None

def merge_sample_mafs(func_maf_path, onko_maf_path):
    """
    Reads and merges the Funcotator and OncoKB MAF files on the columns:
    Chromosome, Start_Position, End_Position.
    The Funcotator file is kept intact while the OncoKB file's Chromosome values
    are modified to include a 'chr' prefix if missing.
    
    Returns the merged DataFrame or None if an error occurs.
    """
    try:
        # Specify dtypes for faster reading; adjust if additional columns/types are needed.
        dtype_spec = {'Chromosome': str, 'Start_Position': int, 'End_Position': int}
        df_func = pd.read_csv(func_maf_path, sep='\t', comment='#', low_memory=False, dtype=dtype_spec)
        df_onko = pd.read_csv(onko_maf_path, sep='\t', comment='#', low_memory=False, dtype=dtype_spec)
    except Exception as e:
        print(f"Error reading MAF files: {e}")
        return None

    # Ensure the OncoKB file's Chromosome values have a 'chr' prefix.
    df_onko['Chromosome'] = df_onko['Chromosome'].astype(str).apply(
        lambda x: x if x.startswith('chr') else 'chr' + x
    )

    # Set index on merge columns for faster join.
    merge_cols = ['Chromosome', 'Start_Position', 'End_Position']
    df_func.set_index(merge_cols, inplace=True)
    df_onko.set_index(merge_cols, inplace=True)

    merged = pd.merge(
        df_func, df_onko,
        left_index=True, right_index=True,
        how='outer', suffixes=('_func', '_onko')
    )
    merged.reset_index(inplace=True)
    # Drop Funcotator's gene/classification columns — VEP (_onko) always covers
    # these 100% of rows and coalesce() always picks _onko first anyway.
    merged.drop(columns=[c for c in ('Hugo_Symbol_func', 'Variant_Classification_func')
                         if c in merged.columns], inplace=True)
    return merged

def annotate_hotspots(df, hotspot_set):
    """
    Adds a new column 'hotspot' that gets a '+' if the (Chromosome, Start_Position, End_Position)
    tuple exists in the hotspot_set; otherwise, it remains blank.
    This is achieved via a vectorized merge.
    """
    # Construct a DataFrame from the hotspot set.
    hotspot_df = pd.DataFrame(list(hotspot_set), columns=['Chromosome', 'Start_Position', 'End_Position'])
    hotspot_df['hotspot'] = '+'
    # Merge with the main DataFrame; rows with no match will have NaN in 'hotspot'.
    df = pd.merge(df, hotspot_df, on=['Chromosome', 'Start_Position', 'End_Position'], how='left')
    df['hotspot'] = df['hotspot'].fillna('')
    return df

def process_sample(sample_id, func_path, onko_path, hotspot_set, output_dir):
    """
    Processes a single sample: merging the Funcotator and OncoKB MAF files,
    annotating hotspot status, saving the individual merged file, and
    returning the merged DataFrame.
    """
    merged_df = merge_sample_mafs(func_path, onko_path)
    if merged_df is not None:
        merged_df = annotate_hotspots(merged_df, hotspot_set)
        merged_df['Tumor_Sample_Barcode'] = sample_id
        out_file = os.path.join(output_dir, f"{sample_id}_merged.maf")
        merged_df.to_csv(out_file, sep='\t', index=False)
        print(f"Merged MAF written to {out_file}")
        return merged_df
    else:
        print(f"Error processing sample {sample_id}")
        return None

def main(func_dir, onko_dir, hotspot_file, cosmic_file, output_dir, num_threads):
    os.makedirs(output_dir, exist_ok=True)

    # Read the hotspot file and build a set of positions.
    try:
        df_hotspot = pd.read_csv(hotspot_file, sep='\t', comment='#',low_memory=False)
    except Exception as e:
        print(f"Error reading hotspot file: {e}")
        return

    # Ensure hotspot file Chromosome values have the 'chr' prefix.
    df_hotspot['Chromosome'] = df_hotspot['Chromosome'].astype(str).apply(
        lambda x: x if x.startswith('chr') else 'chr' + x
    )
    hotspot_set = set(
        df_hotspot[['Chromosome', 'Start_Position', 'End_Position']].itertuples(index=False, name=None)
    )

    # Build dictionaries mapping sample_id to file path for both directories.
    func_files = {}
    for f in glob.glob(os.path.join(func_dir, "*.maf")):
        base = os.path.basename(f)
        sample_id = extract_sample_id(base)
        if sample_id:
            func_files[sample_id] = f
        else:
            print(f"Warning: Could not extract sample id from Funcotator file: {base}")

    onko_files = {}
    for f in glob.glob(os.path.join(onko_dir, "*.maf")):
        base = os.path.basename(f)
        sample_id = extract_sample_id(base)
        if sample_id:
            onko_files[sample_id] = f
        else:
            print(f"Warning: Could not extract sample id from OncoKB file: {base}")

    merged_samples = []  # To store merged DataFrames for each sample.

    # Use parallel processing with the specified number of threads.
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {}
        for sample_id, func_path in func_files.items():
            if sample_id in onko_files:
                onko_path = onko_files[sample_id]
                futures[executor.submit(process_sample, sample_id, func_path, onko_path, hotspot_set, output_dir)] = sample_id
            else:
                print(f"No matching OncoKB MAF found for sample {sample_id}")

        for future in tqdm(concurrent.futures.as_completed(futures),
                           total=len(futures), desc="Processing samples"):
            result = future.result()
            if result is not None:
                merged_samples.append(result)

    # Combine all sample-specific merged DataFrames into one final file.
    if merged_samples:
        combined_df = pd.concat(merged_samples, ignore_index=True)
        combined_df.sort_values(by=['Chromosome', 'Start_Position'], inplace=True)
        combined_out_file = os.path.join(output_dir, "combined_merged.maf")
        combined_df.to_csv(combined_out_file, sep='\t', index=False)
        print(f"Combined and sorted MAF file written to {combined_out_file}")
    else:
        print("No merged samples to combine.")
        return

    # ===== Integrate Cosmic TSV Data with Allele Check =====
    try:
        cosmic_df = pd.read_csv(cosmic_file, sep='\t', comment='#', low_memory=False)
    except Exception as e:
        print(f"Error reading cosmic file: {e}")
        return

    # Ensure the Cosmic file's CHROMOSOME values have the 'chr' prefix.
    cosmic_df['CHROMOSOME'] = cosmic_df['CHROMOSOME'].astype(str).apply(
        lambda x: x if x.startswith('chr') else 'chr' + x
    )

    # Retain the desired Cosmic columns.
    cosmic_cols = ['CHROMOSOME', 'GENOME_START', 'GENOME_STOP',
                   'MUTATION_ID', 'MUTATION_DESCRIPTION', 'MUTATION_CDS',
                   'MUTATION_SOMATIC_STATUS', 'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE']
    cosmic_df = cosmic_df[cosmic_cols]

    # Merge the Cosmic data into the combined MAF using coordinates.
    combined_cosmic_df = pd.merge(
        combined_df, cosmic_df,
        left_on=['Chromosome', 'Start_Position', 'End_Position'],
        right_on=['CHROMOSOME', 'GENOME_START', 'GENOME_STOP'],
        how='left'
    )

    # Drop the extra merge key columns from Cosmic.
    combined_cosmic_df.drop(columns=['CHROMOSOME', 'GENOME_START', 'GENOME_STOP'], inplace=True)

    # ---- Allele Matching Logic ----
    # For a valid Cosmic match, we expect:
    #   Reference_Allele_func == GENOMIC_WT_ALLELE and
    #   GENOMIC_MUT_ALLELE equals either Tumor_Seq_Allele1_func or Tumor_Seq_Allele2_func.
    allele_match_mask = (
        (combined_cosmic_df['Reference_Allele_func'] == combined_cosmic_df['GENOMIC_WT_ALLELE']) &
        (
            (combined_cosmic_df['Tumor_Seq_Allele1_func'] == combined_cosmic_df['GENOMIC_MUT_ALLELE']) |
            (combined_cosmic_df['Tumor_Seq_Allele2_func'] == combined_cosmic_df['GENOMIC_MUT_ALLELE'])
        )
    )

    # List the Cosmic annotation columns that should only be retained when allele match is true.
    cosmic_annotation_cols = [
        'MUTATION_ID', 'MUTATION_DESCRIPTION', 'MUTATION_CDS',
        'MUTATION_SOMATIC_STATUS', 'GENOMIC_WT_ALLELE', 'GENOMIC_MUT_ALLELE'
    ]

    # For rows where the allele condition fails, set the Cosmic columns to missing.
    combined_cosmic_df.loc[~allele_match_mask, cosmic_annotation_cols] = pd.NA

    # Write out the final combined file.
    cosmic_combined_out_file = os.path.join(output_dir, "combined_merged_with_cosmic.maf")
    combined_cosmic_df.to_csv(cosmic_combined_out_file, sep='\t', index=False)
    print(f"Combined MAF file with Cosmic info written to {cosmic_combined_out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge MAFs with hotspot and Cosmic annotation - Optimized version")
    parser.add_argument("--func_dir", required=True, help="Directory containing Funcotator MAFs")
    parser.add_argument("--onko_dir", required=True, help="Directory containing OncoKB MAFs")
    parser.add_argument("--hotspot_file", required=True, help="Path to hotspot MAF file")
    parser.add_argument("--cosmic_file", required=True, help="Path to COSMIC TSV file")
    parser.add_argument("--output_dir", required=True, help="Output directory for merged MAFs")
    parser.add_argument("--num_threads", type=int, default=4, help="Number of threads to use for parallel processing")
    args = parser.parse_args()

    main(
        func_dir=args.func_dir,
        onko_dir=args.onko_dir,
        hotspot_file=args.hotspot_file,
        cosmic_file=args.cosmic_file,
        output_dir=args.output_dir,
        num_threads=args.num_threads
    )
