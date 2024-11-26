# Gaussian Smoothing for BedGraph Files

This repository contains a Python script for performing Gaussian smoothing on genomic data stored in BedGraph format. The script is optimized for noise reduction and signal enhancement in genomic datasets, such as ChIP-seq or CUT&RUN data, by applying Gaussian smoothing to user-defined chromosome intervals.

## Features
- **Chromosome-Specific Smoothing**: Process data for a specific chromosome to enhance flexibility and efficiency.
- **Gaussian Smoothing**: Apply smoothing with a user-defined window size and standard deviation for noise reduction.
- **BedGraph Expansion**: Expand variable-length intervals to 1-nt resolution for accurate smoothing.
- **Value Collapsing**: Collapse consecutive intervals with identical smoothed values (rounded to two decimals) to reduce file size.

## Usage

### Input Requirements
The script expects a BedGraph file with the following columns:
- `chrom`: Chromosome name.
- `start`: Start coordinate of the interval.
- `end`: End coordinate of the interval.
- `value`: Signal intensity for the interval.

### Command-Line Arguments
```bash
python gaussian_smoothing_for_bg_specific_chrom.py <input_file> <window_size> <std_dev> <output_file> <chrom>
```

###Parameters
 - input_file: Path to the input BedGraph file.
 - window_size: Window size for smoothing (e.g., 50 nt).
 - std_dev: Standard deviation for Gaussian smoothing (e.g., 5 nt).
 - output_file: Path to the output BedGraph file with smoothed values.
 - chrom: Chromosome to process (e.g., chr7).

###Implementation Details
1. Expansion: Expands intervals longer than 500 bp to 1-nt resolution for precise smoothing.
2. Smoothing: Applies a Gaussian filter to the signal intensity across the expanded intervals.
3. Collapsing: Merges consecutive intervals with the same smoothed value, rounded to two decimals, for compact output.

###Citation
If you use this script in your research, please cite:
Kim et al., ??, 2025.
