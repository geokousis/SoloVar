#!/usr/bin/env python3
import sys
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd


def calculate_percentile(data, value):
    size = len(data)
    if size == 0:
        return 0.0
    sorted_data = sorted(data)
    count_below = sum(1 for x in sorted_data if x < value)
    return (count_below / size) * 100


def parse_combined_log(file_path):
    """
    Parse a combined flagstat log file where each sample block is introduced
    by a header line of the form:  === LOG_<samplename>.txt ===

    samtools flagstat lines are parsed by content (keyword match) so the
    function is robust to both the 13-line (samtools < 1.13) and the
    16-line (samtools >= 1.13) output formats.
    """
    with open(file_path, "r") as fh:
        raw = fh.read()

    # Split into per-sample blocks on the === header lines
    blocks = re.split(r"^=== .+ ===$", raw, flags=re.MULTILINE)
    headers = re.findall(r"^=== (.+) ===$", raw, flags=re.MULTILINE)

    sample_names = []
    mapping_numbers = []
    properly_paired_numbers = []
    properly_paired_percentages = []

    for header, block in zip(headers, blocks[1:]):  # blocks[0] is text before first header (empty)
        sample_name = os.path.basename(header.strip().replace("LOG_", "").replace(".txt", ""))
        lines = [l for l in block.splitlines() if l.strip()]

        mapped_count = None
        paired_count = None
        paired_pct = None

        for line in lines:
            # "mapped" line — prefer "primary mapped" if present, fall back to "mapped"
            if re.search(r"\bprimary mapped\b", line):
                m = re.match(r"(\d+)\s*\+\s*\d+\s+primary mapped", line)
                if m:
                    mapped_count = int(m.group(1))
            elif mapped_count is None and re.search(r"\bmapped\b", line) and "supplementary" not in line:
                m = re.match(r"(\d+)\s*\+\s*\d+\s+mapped", line)
                if m:
                    mapped_count = int(m.group(1))

            # "properly paired" line
            if re.search(r"\bproperly paired\b", line):
                m = re.match(r"(\d+)\s*\+\s*\d+\s+properly paired\s+\(([0-9.]+)%", line)
                if m:
                    paired_count = int(m.group(1))
                    paired_pct = float(m.group(2))

        if mapped_count is None or paired_count is None:
            print(f"Warning: could not parse flagstat block for sample '{sample_name}'; skipping.")
            continue

        sample_names.append(sample_name)
        mapping_numbers.append(mapped_count)
        properly_paired_numbers.append(paired_count)
        properly_paired_percentages.append(paired_pct)
        print(f"Processed sample: {sample_name}  mapped={mapped_count:,}  paired={paired_count:,}  paired%={paired_pct:.2f}")

    return sample_names, mapping_numbers, properly_paired_numbers, properly_paired_percentages


def process_file_and_plot(file_path, output_dir):
    sample_names, mapping_numbers, properly_paired_numbers, properly_paired_percentages = \
        parse_combined_log(file_path)

    if not sample_names:
        print("Error: no samples were successfully parsed from the log file.")
        sys.exit(1)

    data = pd.DataFrame({
        "Sample": sample_names,
        "Mapping Number": mapping_numbers,
        "Properly Paired Number": properly_paired_numbers,
        "Properly Paired Percentage": properly_paired_percentages,
    })

    # Save statistics to CSV
    csv_file = os.path.join(output_dir, "samples_stats.txt")
    data[["Sample", "Mapping Number", "Properly Paired Number", "Properly Paired Percentage"]].to_csv(
        csv_file, index=False
    )

    # Convert raw counts to millions for plotting
    data["Mapping Number (Millions)"] = data["Mapping Number"] / 1e6
    data["Properly Paired Number (Millions)"] = data["Properly Paired Number"] / 1e6
    data["Color"] = data["Mapping Number"].apply(lambda x: "orange" if x <= 1_000_000 else "blue")

    THRESHOLD = 1_000_000
    filtered_data = data[data["Mapping Number"] > THRESHOLD]

    threshold_percentile = calculate_percentile(mapping_numbers, THRESHOLD)
    percentile_10th = calculate_percentile(mapping_numbers, 0.10)
    print(f"10th Percentile of Mapping Number: {percentile_10th:.2f}%")
    print(f"Percentile for threshold {THRESHOLD:,}: {threshold_percentile:.2f}%")

    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["font.family"] = "sans-serif"

    fig, axs = plt.subplots(3, 1, figsize=(20, 28), gridspec_kw={"height_ratios": [1.5, 3, 1]})

    # --- Subplot 1: Histogram of mapping numbers ---
    max_millions = data["Mapping Number (Millions)"].max()
    bins = range(0, int(max_millions) + 2)
    counts, edges, patches = axs[0].hist(data["Mapping Number (Millions)"], bins=bins, alpha=0.5)
    for patch, edge in zip(patches, edges):
        patch.set_facecolor("orange" if edge < 1 else "blue")
    axs[0].axvline(1, color="red", linestyle="dashed", linewidth=1)
    axs[0].set_title("Histogram of Primary Mapped Reads")
    axs[0].set_xlabel("Primary Mapped Reads (Millions)")
    axs[0].set_ylabel("Counts")
    axs[0].set_xticks(range(0, int(max_millions) + 2))
    ylim_top = axs[0].get_ylim()[1]
    axs[0].annotate(
        f"{threshold_percentile:.2f} percentile",
        xy=(1, ylim_top * 0.9),
        xytext=(1.2, ylim_top * 0.9),
        arrowprops=dict(facecolor="black", shrink=0.05),
        fontsize=12,
        color="black",
    )

    # --- Subplot 2: Bar plot mapped + properly paired per sample ---
    bar_colors = data["Color"].tolist()
    x_pos = range(len(data))

    axs[1].bar(x_pos, data["Mapping Number (Millions)"], color=bar_colors, alpha=0.8, label="_nolegend_")
    axs[1].bar(x_pos, data["Properly Paired Number (Millions)"], color="green", alpha=0.5, label="_nolegend_")
    axs[1].set_xticks(list(x_pos))
    axs[1].set_xticklabels(data["Sample"].tolist(), rotation=90, ha="center", fontsize=12)
    axs[1].set_title("Primary Mapped and Properly Paired Reads per Sample")
    axs[1].set_xlabel("Samples")
    axs[1].set_ylabel("Reads (Millions)")

    legend_handles = [
        mpatches.Patch(color="blue", label="Mapped Reads > 1 Million"),
        mpatches.Patch(color="orange", label="Mapped Reads \u2264 1 Million"),
        mpatches.Patch(color="green", alpha=0.5, label="Properly Paired Reads"),
    ]
    axs[1].legend(handles=legend_handles, loc="upper right")

    # --- Subplot 3: Histogram of properly paired percentages ---
    if not filtered_data.empty:
        sns.histplot(
            filtered_data["Properly Paired Percentage"],
            bins=30,
            ax=axs[2],
            color="green",
            kde=False,
        )
    else:
        axs[2].text(0.5, 0.5, "No samples above 1M reads", ha="center", va="center",
                    transform=axs[2].transAxes)
    axs[2].set_title("Histogram of Properly Paired Reads Percentages (Samples > 1 Million)")
    axs[2].set_xlabel("Properly Paired Reads Percentage")
    axs[2].set_ylabel("Counts")

    plt.tight_layout()
    plot_file = os.path.join(output_dir, "combined_mapping_plot_1.png")
    plt.savefig(plot_file, dpi=150)
    print(f"Plot saved to {plot_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: qc_pipeline.py <combined_log_file> <output_directory>")
        sys.exit(1)
    process_file_and_plot(sys.argv[1], sys.argv[2])
