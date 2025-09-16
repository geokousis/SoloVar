#!/usr/bin/env python3                                                                                                                                                                                                                                                                                                                                                                                        
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def calculate_percentile(data, value):
    size = len(data)
    sorted_data = sorted(data)
    count_below = len([x for x in sorted_data if x < value])
    return (count_below / size) * 100

def process_file_and_plot(file_path, output_dir):
    # Read the combined log file                                                                                                                                                                                                                                                                                                                                                                              
    with open(file_path, 'r') as file:
        lines = file.readlines()

    sample_names = []
    mapping_numbers = []
    properly_paired_numbers = []
    properly_paired_percentages = []

    # Process each block of 19 lines.                                                                                                                                                                                                                                                                                                                                                                         
    for i in range(0, len(lines), 19):
        try:
            sample_name = os.path.basename(lines[i+1].strip().replace("LOG_", "").replace(".txt", ""))
            sample_names.append(sample_name)
            mapping_info = lines[i+10].strip()
            mapped_number = int(mapping_info.split()[0].split('+')[0])
            mapping_numbers.append(mapped_number)
            properly_paired_info = lines[i+14].strip()
            properly_paired_number = int(properly_paired_info.split()[0].split('+')[0])
            properly_paired_percentage = float(properly_paired_info.split('(')[1].split('%')[0])
            properly_paired_numbers.append(properly_paired_number)
            properly_paired_percentages.append(properly_paired_percentage)
            print(f"Processed sample: {sample_name}")
        except (IndexError, ValueError) as e:
            print(f"Error processing block starting at line {i}: {e}")

    # Create a DataFrame from the processed data.                                                                                                                                                                                                                                                                                                                                                             
    data = pd.DataFrame({
        'Sample': sample_names,
        'Mapping Number': mapping_numbers,
        'Properly Paired Number': properly_paired_numbers,
        'Properly Paired Percentage': properly_paired_percentages,
        'Color': ['orange' if x > 1000000 else 'blue' for x in mapping_numbers]
    })

    # Convert raw counts to millions for plotting.                                                                                                                                                                                                                                                                                                                                                            
    data['Mapping Number (Millions)'] = data['Mapping Number'] / 1e6
    data['Properly Paired Number (Millions)'] = data['Properly Paired Number'] / 1e6
    filtered_data = data[data['Mapping Number'] > 1000000]

    # Save statistics to CSV.                                                                                                                                                                                                                                                                                                                                                                                 
    csv_file = f"{output_dir}/samples_stats.txt"
    data[['Sample', 'Mapping Number', 'Properly Paired Number', 'Properly Paired Percentage']].to_csv(csv_file, index=False)

    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["font.family"] = "sans-serif"

    # Create a figure with three subplots.                                                                                                                                                                                                                                                                                                                                                                    
    fig, axs = plt.subplots(3, 1, figsize=(20,28), gridspec_kw={'height_ratios': [1.5, 3, 1]})

    # Plot histogram of mapping numbers.                                                                                                                                                                                                                                                                                                                                                                      
    bins = range(0, int(data['Mapping Number (Millions)'].max())+2)
    counts, edges, patches = axs[0].hist(data['Mapping Number (Millions)'], bins=bins, alpha=0.5)
    for patch, edge in zip(patches, edges):
        patch.set_facecolor('orange' if edge < 1 else 'blue')
    axs[0].axvline(1, color='red', linestyle='dashed', linewidth=1)
    axs[0].set_title('Histogram of Primary Mapped Reads')
    axs[0].set_xlabel('Primary Mapped Reads (Millions)')
    axs[0].set_ylabel('Counts')
    max_mapping_value = int(data['Mapping Number (Millions)'].max()) + 1
    axs[0].set_xticks(range(0, max_mapping_value+1))
    axs[0].set_yticks(range(0, int(axs[0].get_ylim()[1]) + 1))
    threshold = 1000000
    threshold_percentile = calculate_percentile(mapping_numbers, threshold)
    axs[0].annotate(f'{threshold_percentile:.2f} percentile',
                    xy=(1, axs[0].get_ylim()[1] * 0.9),
                    xytext=(1.2, axs[0].get_ylim()[1] * 0.9),
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    fontsize=12, color='black')

    # Bar plot of mapped and properly paired reads.                                                                                                                                                                                                                                                                                                                                                           
    sns.barplot(x='Sample', y='Mapping Number (Millions)', hue='Color', dodge=False,
                data=data, palette={'orange': 'orange', 'blue': 'blue'}, ax=axs[1])
    for index, row in data.iterrows():
        axs[1].bar(row['Sample'], row['Properly Paired Number (Millions)'], color='green', alpha=0.5)
    axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=90, ha="center", fontsize=12)
    axs[1].set_title('Primary Mapped and Properly Paired Reads per Sample')
    axs[1].set_xlabel('Samples')
    axs[1].set_ylabel('Reads (Millions)')

    # Fix for missing legend handles                                                                                                                                                                                                                                                                                                                                                                          
    handles, labels = axs[1].get_legend_handles_labels()

    # Ensure there are at least two handles before modifying them                                                                                                                                                                                                                                                                                                                                             
    if len(handles) >= 2:
        handles[0] = plt.Rectangle((0, 0), 1, 1, color='orange')
        handles[1] = plt.Rectangle((0, 0), 1, 1, color='blue')

    # Always add the 'Properly Paired Reads' handle                                                                                                                                                                                                                                                                                                                                                           
    handles.append(plt.Rectangle((0, 0), 1, 1, color='green', alpha=0.5))

    # Ensure labels match handles                                                                                                                                                                                                                                                                                                                                                                             
    labels = ['Mapped Reads > 1 Million', 'Mapped Reads < 1 Million', 'Properly Paired Reads']

    # Set the legend only if handles exist                                                                                                                                                                                                                                                                                                                                                                    
    if handles:
        axs[1].legend(handles=handles, labels=labels, loc='upper right')

    # Histogram of properly paired read percentages.                                                                                                                                                                                                                                                                                                                                                          
    sns.histplot(filtered_data['Properly Paired Percentage'], bins=30, ax=axs[2], color='green', kde=False)
    axs[2].set_title('Histogram of Properly Paired Reads Percentages (Samples > 1 Million)')
    axs[2].set_xlabel('Properly Paired Reads Percentage')
    axs[2].set_ylabel('Counts')

    plt.tight_layout()
    plot_file = f"{output_dir}/combined_mapping_plot_1.png"
    plt.savefig(plot_file)

    percentile_10th = calculate_percentile(mapping_numbers, 0.10)
    print(f"10th Percentile of Mapping Number: {percentile_10th}")
    print(f"Percentile for the threshold value {threshold}: {threshold_percentile:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: qc_pipeline.py <combined_log_file> <output_directory>")
        sys.exit(1)
    combined_log = sys.argv[1]
    output_dir = sys.argv[2]
    process_file_and_plot(combined_log, output_dir)
