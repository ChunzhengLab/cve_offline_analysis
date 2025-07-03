#!/usr/bin/env python3
"""
Feeddown comparison plotting script
This script compares observables before and after feeddown correction
Extracted from plot_stepbystep.py feeddown comparison functionality
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from pathlib import Path
import argparse

# Set up modern plotting style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif'],
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.figsize': (10, 8),
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 2.0,
    'errorbar.capsize': 3,
})

# Define color schemes
COLORS = {
    'original_SS': '#1f77b4',     # blue - original SS
    'original_OS': '#ff7f0e',     # orange - original OS  
    'original_Del': '#2ca02c',    # green - original Del
    'feeddown_SS': '#d62728',     # red - feeddown SS
    'feeddown_OS': '#9467bd',     # purple - feeddown OS
    'feeddown_Del': '#8c564b',    # brown - feeddown Del
}

def load_data(file_path):
    """
    Load data from CSV file
    
    Args:
        file_path: Path to CSV file
        
    Returns:
        DataFrame or None if file not found
    """
    try:
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}")
            return None
            
        df = pd.read_csv(file_path)
        
        # Remove resolution column if present
        if 'resolution' in df.columns:
            df = df.drop(columns=['resolution'])
            
        # Filter centrality range (0-60)
        df = df[df['centrality'] <= 60]
        
        return df
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return None

def prepare_output_directory(directory):
    """Create output directory if it doesn't exist"""
    Path(directory).mkdir(parents=True, exist_ok=True)
    print(f"Output directory prepared: {directory}")

def plot_feeddown_comparison(original_data_path, feeddown_data_path, 
                           diff_type=None, diff_bin=None, output_dir="./plots"):
    """
    Plot feeddown comparison for both delta and gamma observables
    
    Args:
        original_data_path: Path to original data
        feeddown_data_path: Path to feeddown corrected data
        diff_type: Differential type (e.g. 'DEta', 'SPt', 'Intg')
        diff_bin: Differential bin value
        output_dir: Output directory
    """
    # Load data
    original_data = load_data(original_data_path)
    feeddown_data = load_data(feeddown_data_path)
    
    if original_data is None or feeddown_data is None:
        print("Skipping plot: missing data files")
        return
        
    # Prepare figure with 2 subplots (delta and gamma)
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"Lambda-Proton Correlation: Feeddown Comparison", fontsize=18)
    
    # Set subplot titles
    axes[0].set_title("Delta Observable", fontsize=16)
    axes[1].set_title("Gamma Observable", fontsize=16)
    
    # Set axis labels and ranges
    for ax in axes:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_xlim(0, 60)
        ax.grid(True, alpha=0.3)
    
    axes[0].set_ylabel('Delta Value', fontsize=14)
    axes[1].set_ylabel('Gamma Value', fontsize=14)
    
    # Apply filters
    filter_cond = True
    if diff_type:
        filter_cond = (original_data['diff_type'] == diff_type) & (feeddown_data['diff_type'] == diff_type)
    if diff_bin is not None:
        filter_cond = filter_cond & (original_data['diff_bin'] == diff_bin) & (feeddown_data['diff_bin'] == diff_bin)
    
    # Filter data for each pair type
    pair_types = ['SS', 'OS', 'Del']
    filtered_data = {}
    
    for pair_type in pair_types:
        filtered_data[f'original_{pair_type}'] = original_data[filter_cond & (original_data['pair_type'] == pair_type)]
        filtered_data[f'feeddown_{pair_type}'] = feeddown_data[filter_cond & (feeddown_data['pair_type'] == pair_type)]
    
    # Plot for both delta and gamma
    for i, value_type in enumerate(['delta', 'gamma']):
        ax = axes[i]
        
        # Plot each pair type
        for pair_type in pair_types:
            orig_data = filtered_data[f'original_{pair_type}']
            feed_data = filtered_data[f'feeddown_{pair_type}']
            
            # Create proper labels
            pair_label = "OS-SS" if pair_type == "Del" else pair_type
            
            # Plot original data
            ax.errorbar(
                orig_data['centrality'], orig_data[value_type], orig_data[f"{value_type}_err"],
                marker='o', linestyle='-', color=COLORS[f'original_{pair_type}'], 
                label=f"Before feeddown removal {pair_label}", alpha=0.8, markerfacecolor='none'
            )
            
            # Plot feeddown corrected data
            ax.errorbar(
                feed_data['centrality'], feed_data[value_type], feed_data[f"{value_type}_err"],
                marker='s', linestyle='--', color=COLORS[f'feeddown_{pair_type}'], 
                label=f"After feeddown removal {pair_label}", alpha=0.8
            )
        
        ax.legend(loc='best', frameon=True, framealpha=0.9)
    
    # Create output directory
    prepare_output_directory(output_dir)
    
    # Save plot
    filename = f"feeddown_comparison"
    if diff_type:
        filename += f"_{diff_type}"
    if diff_bin is not None:
        filename += f"_{diff_bin}"
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir}/{filename}.pdf")

def plot_all_diff_types_comparison(original_data_path, feeddown_data_path, output_dir="./plots"):
    """
    Plot feeddown comparison for all differential types
    
    Args:
        original_data_path: Path to original data
        feeddown_data_path: Path to feeddown corrected data
        output_dir: Output directory
    """
    # Load data
    original_data = load_data(original_data_path)
    feeddown_data = load_data(feeddown_data_path)
    
    if original_data is None or feeddown_data is None:
        print("Skipping plot: missing data files")
        return
    
    # Get all available differential type and bin combinations
    diff_combinations = []
    
    # Get all combinations from original data
    for diff_type in original_data['diff_type'].unique():
        for diff_bin in original_data[original_data['diff_type'] == diff_type]['diff_bin'].unique():
            diff_combinations.append((diff_type, diff_bin))
    
    # Generate plot for each combination
    for diff_type, diff_bin in diff_combinations:
        print(f"Generating comparison plot for {diff_type} {diff_bin}...")
        plot_feeddown_comparison(
            original_data_path, feeddown_data_path, 
            diff_type, diff_bin, output_dir
        )

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Feeddown comparison plotting")
    parser.add_argument("--original", type=str, required=True, help="Path to original data file")
    parser.add_argument("--feeddown", type=str, required=True, help="Path to feeddown corrected data file")
    parser.add_argument("--diff_type", type=str, help="Differential type (e.g. DEta, SPt, Intg)")
    parser.add_argument("--diff_bin", type=float, help="Differential bin value")
    parser.add_argument("--output_dir", type=str, default="./plots", help="Output directory")
    parser.add_argument("--all_types", action='store_true', help="Generate plots for all differential types")
    
    args = parser.parse_args()
    
    if args.all_types:
        # Generate comparison plots for all differential types
        print("Generating comparison plots for all differential types...")
        plot_all_diff_types_comparison(
            args.original, args.feeddown, args.output_dir
        )
    else:
        # Generate single comparison plot
        print("Generating feeddown comparison plot...")
        plot_feeddown_comparison(
            args.original, args.feeddown, 
            args.diff_type, args.diff_bin, args.output_dir
        )
    
    print("Plotting completed!")

if __name__ == "__main__":
    main()