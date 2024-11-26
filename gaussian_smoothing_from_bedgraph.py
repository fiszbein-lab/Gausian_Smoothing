import argparse
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d

def expand_bedgraph(df):
    """Expands variable-length bedGraph intervals to 1-nt resolution for intervals longer than 500 bp."""
    expanded_data = []
    for _, row in df.iterrows():
        chrom = str(row['chrom'])  # Ensure chromosome is treated as a string
        start = int(row['start'])  # Ensure start is an integer
        end = int(row['end'])      # Ensure end is an integer
        value = row['value']
        
        # Only expand if interval is longer than 500 bp
        if (end - start) > 500:
            for pos in range(start, end):
                expanded_data.append([chrom, pos, pos + 1, value])
        else:
            # If interval <= 500 bp, keep as-is
            expanded_data.append([chrom, start, end, value])
    
    expanded_df = pd.DataFrame(expanded_data, columns=['chrom', 'start', 'end', 'value'])
    return expanded_df

def apply_gaussian_smoothing(expanded_df, std_dev):
    """Apply Gaussian smoothing to the entire data, including zero regions."""
    values = expanded_df['value'].to_numpy()
    
    # Apply Gaussian smoothing across the entire series of values
    smoothed_values = gaussian_filter1d(values, sigma=std_dev)
    
    # Assign smoothed values back to the DataFrame
    expanded_df['smoothed_value'] = smoothed_values
    return expanded_df

def collapse_bedgraph(df):
    """Collapse consecutive rows with the same smoothed value, rounded to two decimal places."""
    collapsed_data = []
    prev_chrom = df.iloc[0]['chrom']
    prev_start = df.iloc[0]['start']
    prev_end = df.iloc[0]['end']
    prev_value = round(df.iloc[0]['smoothed_value'], 2)  # Round the initial value

    for _, row in df.iloc[1:].iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        value = round(row['smoothed_value'], 2)  # Round the current value
        
        if chrom == prev_chrom and value == prev_value:
            # Extend the previous interval if the rounded value is the same
            prev_end = end
        else:
            # Append the previous interval to collapsed_data
            collapsed_data.append([prev_chrom, prev_start, prev_end, prev_value])
            
            # Start a new interval
            prev_chrom = chrom
            prev_start = start
            prev_end = end
            prev_value = value

    # Append the last interval
    collapsed_data.append([prev_chrom, prev_start, prev_end, prev_value])

    # Convert to DataFrame
    collapsed_df = pd.DataFrame(collapsed_data, columns=['chrom', 'start', 'end', 'value'])
    return collapsed_df

def main(input_file, window_size, std_dev, output_file, chrom_to_work):
    # Load bedGraph data, filter only for the chrom of interest
    df = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'], dtype={'chrom': str, 'start': int, 'end': int, 'value': float})
    df_chr_to_work = df[df['chrom'] == chrom_to_work].reset_index(drop=True)
    
    # Process: expand, smooth, and collapse
    expanded_df = expand_bedgraph(df_chr_to_work)
    smoothed_df = apply_gaussian_smoothing(expanded_df, std_dev)
    collapsed_df = collapse_bedgraph(smoothed_df)
    
    # Write the final collapsed data to output file
    collapsed_df.to_csv(output_file, sep='\t', header=False, index=False, float_format='%.2f')
    print(f"Smoothing complete for {chrom_to_work}. Smoothed and collapsed bedGraph saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimized Gaussian smoothing for bedGraph files.")
    parser.add_argument("input_file", help="Path to the input bedGraph file.")
    parser.add_argument("window_size", type=int, help="Window size for smoothing.")
    parser.add_argument("std_dev", type=float, help="Standard deviation for Gaussian smoothing.")
    parser.add_argument("output_file", help="Path to the output smoothed bedGraph file.")
    parser.add_argument("chrom", help="Chromosome to smooth")
    args = parser.parse_args()
    
    main(args.input_file, args.window_size, args.std_dev, args.output_file, args.chrom)
