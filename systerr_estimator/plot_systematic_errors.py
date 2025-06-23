#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from pathlib import Path
import matplotlib.backends.backend_pdf

def get_default_data(input_file):
    """
    Get default data from the corresponding original file
    
    Parameters:
    -----------
    input_file : str
        Path to the barlow input file
        
    Returns:
    --------
    pandas.DataFrame
        Default data
    """
    # Derive the original file path from the barlow file
    input_path = Path(input_file)
    particle_type = input_path.stem.replace('barlow_finalise_sys_', '')
    
    # Try to find the original file
    original_file = input_path.parent / f"finalise_sys_{particle_type}.csv"
    if not original_file.exists():
        print(f"Error: Original file {original_file} not found")
        
        # Try to look in parent directory as fallback
        original_file = input_path.parent.parent / f"finalise_sys_{particle_type}.csv"
        if not original_file.exists():
            print(f"Error: Could not find original file in parent directory either")
            return None
    
    print(f"Reading default data from: {original_file}")
    try:
        # Read the original file
        original_df = pd.read_csv(original_file)
        
        # Check if 'source' column exists
        if 'source' not in original_df.columns:
            print(f"Error: Original file {original_file} doesn't have a 'source' column")
            return None
        
        # Extract default data
        default_df = original_df[original_df['source'] == 'default']
        
        if len(default_df) == 0:
            print(f"Error: No default data found in {original_file}")
            return None
            
        print(f"Found {len(default_df)} default data points")
        return default_df
    
    except Exception as e:
        print(f"Error reading original file {original_file}: {e}")
        return None

def create_plots_for_source(df, source_name, default_df, output_dir, particle_type, max_cent=60):
    """
    Create delta and gamma plots for a specific source
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Data for the specific source
    source_name : str
        Name of the source
    default_df : pandas.DataFrame
        Default data for comparison
    output_dir : str
        Output directory for plots
    particle_type : str
        Type of particle (Lambda, Proton, Hadron)
    max_cent : int
        Maximum centrality to plot
    """
    # Filter data by max_cent
    df_filtered = df[df['centrality'] <= max_cent].copy()
    default_filtered = default_df[default_df['centrality'] <= max_cent].copy()
    
    # Get unique diff combinations
    has_diff_type = 'diff_type' in df.columns and not df['diff_type'].isna().all() and not (df['diff_type'] == '').all()
    has_diff_bin = 'diff_bin' in df.columns and not df['diff_bin'].isna().all() and not (df['diff_bin'] == '').all()
    
    print(f"  has_diff_type: {has_diff_type}, has_diff_bin: {has_diff_bin}")
    
    if has_diff_type and has_diff_bin:
        # For Proton/Hadron data, handle both diff_type and diff_bin
        diff_combinations = df_filtered[['diff_type', 'diff_bin']].drop_duplicates().values.tolist()
        diff_combinations = [(str(dt) if pd.notna(dt) else '', str(db) if pd.notna(db) else '') for dt, db in diff_combinations]
        # Remove combinations where both are empty
        diff_combinations = [(dt, db) for dt, db in diff_combinations if dt or db]
    elif has_diff_type:
        diff_combinations = [(str(dt) if pd.notna(dt) else '', '') for dt in df_filtered['diff_type'].unique()]
        diff_combinations = [(dt, db) for dt, db in diff_combinations if dt]
    elif has_diff_bin:
        diff_combinations = [('', str(db) if pd.notna(db) else '') for db in df_filtered['diff_bin'].unique()]
        diff_combinations = [(dt, db) for dt, db in diff_combinations if db]
    else:
        # For Lambda data without meaningful diff columns
        diff_combinations = [('', '')]
    
    print(f"  Found diff combinations: {diff_combinations}")
    
    # Process each diff combination
    for diff_type, diff_bin in diff_combinations:
        # Filter data for this diff combination
        if has_diff_type and has_diff_bin:
            source_data = df_filtered[(df_filtered['diff_type'].astype(str) == diff_type) & 
                                    (df_filtered['diff_bin'].astype(str) == diff_bin)]
            default_data = default_filtered[(default_filtered['diff_type'].astype(str) == diff_type) & 
                                          (default_filtered['diff_bin'].astype(str) == diff_bin)]
            diff_label = f"[{diff_type}]_[{diff_bin}]" if diff_type and diff_bin else ""
        elif has_diff_type:
            source_data = df_filtered[df_filtered['diff_type'].astype(str) == diff_type]
            default_data = default_filtered[default_filtered['diff_type'].astype(str) == diff_type]
            diff_label = f"[{diff_type}]_[]" if diff_type else ""
        elif has_diff_bin:
            source_data = df_filtered[df_filtered['diff_bin'].astype(str) == diff_bin]
            default_data = default_filtered[default_filtered['diff_bin'].astype(str) == diff_bin]
            diff_label = f"[]_[{diff_bin}]" if diff_bin else ""
        else:
            source_data = df_filtered
            default_data = default_filtered
            diff_label = "[]_[]"
        
        print(f"    Processing diff combination: {diff_type}, {diff_bin}")
        print(f"    Source data: {len(source_data)} rows, Default data: {len(default_data)} rows")
        
        if len(source_data) == 0:
            print(f"    No source data for combination {diff_type}, {diff_bin}")
            continue
            
        # Create delta plot
        create_single_plot(source_data, default_data, source_name, 'delta', 
                          diff_type, diff_bin, output_dir, particle_type, diff_label, max_cent)
        
        # Create gamma plot
        create_single_plot(source_data, default_data, source_name, 'gamma', 
                          diff_type, diff_bin, output_dir, particle_type, diff_label, max_cent)

def create_single_plot(source_data, default_data, source_name, value_type, 
                      diff_type, diff_bin, output_dir, particle_type, diff_label, max_cent):
    """
    Create a single plot (delta or gamma) with two subplots
    """
    # Setup figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    # Define colors and markers for pair types
    pair_colors = {'OS': 'red', 'SS': 'blue', 'Del': 'green'}
    pair_markers = {'OS': 'o', 'SS': 's', 'Del': '^'}
    
    # Offset for x-axis positioning
    offset = 1.0
    
    # Get centralities
    centralities = sorted(source_data['centrality'].unique())
    centralities = [c for c in centralities if c <= max_cent]
    
    # First subplot: Values vs centrality
    for pair_type in ['OS', 'SS', 'Del']:
        # Default data
        default_pair = default_data[default_data['pair_type'] == pair_type]
        if len(default_pair) > 0:
            x_default = [c - offset for c in default_pair['centrality']]
            y_default = default_pair[value_type]
            ax1.errorbar(x_default, y_default, yerr=default_pair[f'{value_type}_err'],
                        fmt=pair_markers[pair_type], color=pair_colors[pair_type], 
                        alpha=0.7, label=f'Default {pair_type}', markersize=6)
        
        # Source data
        source_pair = source_data[source_data['pair_type'] == pair_type]
        if len(source_pair) > 0:
            x_source = [c + offset for c in source_pair['centrality']]
            y_source = source_pair[value_type]
            ax1.errorbar(x_source, y_source, yerr=source_pair[f'{value_type}_err'],
                        fmt=pair_markers[pair_type], color=pair_colors[pair_type], 
                        alpha=1.0, label=f'{source_name} {pair_type}', markersize=6,
                        markerfacecolor='none', markeredgewidth=2)
    
    ax1.set_xlabel('Centrality (%)')
    ax1.set_ylabel(f'{value_type.capitalize()}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, max_cent + 5)
    
    # Title for first subplot
    title1 = f'{source_name} - {value_type.capitalize()} vs Centrality'
    if diff_type or diff_bin:
        title1 += f' (diff_type: {diff_type}, diff_bin: {diff_bin})'
    ax1.set_title(title1)
    
    # Second subplot: Barlow ratio vs centrality (histogram style)
    barlow_col = f'{value_type}_barlow_ratio'
    
    # Create histogram-style plot for each pair type
    for pair_type in ['OS', 'SS', 'Del']:
        source_pair = source_data[source_data['pair_type'] == pair_type]
        if len(source_pair) > 0:
            x_values = source_pair['centrality']
            y_values = source_pair[barlow_col]
            
            # Only plot non-NaN values
            valid_mask = ~pd.isna(y_values)
            if valid_mask.sum() > 0:
                x_valid = x_values[valid_mask]
                y_valid = y_values[valid_mask]
                
                # Draw horizontal lines for each bin (histogram style)
                for x, y in zip(x_valid, y_valid):
                    # Draw horizontal line from x-2.5 to x+2.5
                    ax2.plot([x-2.5, x+2.5], [y, y], color=pair_colors[pair_type], 
                            linewidth=3, alpha=0.8, label=f'{pair_type}' if x == x_valid.iloc[0] else "")
                    # Add vertical lines at the edges for clarity
                    ax2.plot([x-2.5, x-2.5], [0, y], color=pair_colors[pair_type], 
                            linewidth=1, alpha=0.5)
                    ax2.plot([x+2.5, x+2.5], [0, y], color=pair_colors[pair_type], 
                            linewidth=1, alpha=0.5)
    
    # Add reference line at y=1
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.7, label='Reference (y=1)')
    
    ax2.set_xlabel('Centrality (%)')
    ax2.set_ylabel(f'{value_type.capitalize()} Barlow Ratio')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, max_cent + 5)
    
    # Title for second subplot
    title2 = f'{source_name} - {value_type.capitalize()} Barlow Ratio vs Centrality'
    if diff_type or diff_bin:
        title2 += f' (diff_type: {diff_type}, diff_bin: {diff_bin})'
    ax2.set_title(title2)
    
    plt.tight_layout()
    
    # Create particle-specific output directory
    particle_output_dir = os.path.join(output_dir, particle_type)
    os.makedirs(particle_output_dir, exist_ok=True)
    
    # Save plot with improved filename format
    filename = f"{source_name}_{value_type}_{diff_label}.pdf"
    filepath = os.path.join(particle_output_dir, filename)
    plt.savefig(filepath, format='pdf', bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"    Saved plot: {filepath}")

def main():
    parser = argparse.ArgumentParser(description='Generate systematic error plots for each source')
    parser.add_argument('--input', type=str, required=True,
                      help='Input CSV file (e.g., barlow_finalise_sys_Proton.csv)')
    parser.add_argument('--output-dir', type=str, default='./plots',
                      help='Output directory for plots (default: ./plots)')
    parser.add_argument('--max-cent', type=int, default=60,
                      help='Maximum centrality to plot (default: 60)')
    parser.add_argument('--sources', type=str, nargs='*', default=None,
                      help='Specific sources to plot (default: all sources)')
    parser.add_argument('--original', type=str, default=None,
                      help='Path to original file with default data (optional)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read data
    print(f"Reading data from {args.input}")
    df = pd.read_csv(args.input)
    
    # Check required columns
    required_cols = ['source', 'centrality', 'pair_type', 'delta', 'delta_err', 
                    'gamma', 'gamma_err', 'delta_barlow_ratio', 'gamma_barlow_ratio']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {', '.join(missing_cols)}")
        return
    
    # Separate default and source data
    default_df = df[df['source'] == 'default'].copy() if 'source' in df.columns else pd.DataFrame()
    source_df = df[df['source'] != 'default'].copy() if 'source' in df.columns else df.copy()
    
    # If no default data found in the barlow file, try to get it from original file
    if len(default_df) == 0:
        print("No default data found in the barlow file, trying to get from original file...")
        if args.original:
            print(f"Using specified original file: {args.original}")
            try:
                original_df = pd.read_csv(args.original)
                default_df = original_df[original_df['source'] == 'default'].copy()
                if len(default_df) == 0:
                    print(f"Error: No default data found in specified original file {args.original}")
                    return
                else:
                    print(f"Found {len(default_df)} default data points from specified file")
            except Exception as e:
                print(f"Error reading specified original file {args.original}: {e}")
                return
        else:
            default_df = get_default_data(args.input)
            if default_df is None or len(default_df) == 0:
                print("Error: Could not obtain default data")
                return
    
    if len(source_df) == 0:
        print("Warning: No source data found in the file")
        return
    
    # Get unique sources
    sources = source_df['source'].unique()
    
    # Filter sources if specified
    if args.sources:
        sources = [s for s in sources if s in args.sources]
        if not sources:
            print("Error: None of the specified sources found in the data")
            return
    
    print(f"Found {len(sources)} sources to plot")
    print(f"Default data contains {len(default_df)} rows")
    
    # Determine particle type from input filename
    input_path = Path(args.input)
    particle_type = input_path.stem.replace('barlow_finalise_sys_', '')
    print(f"Detected particle type: {particle_type}")
    
    # Create plots for each source
    for source in sources:
        print(f"\nProcessing source: {source}")
        source_data = source_df[source_df['source'] == source]
        
        if len(source_data) == 0:
            print(f"Warning: No data found for source {source}")
            continue
            
        try:
            create_plots_for_source(source_data, source, default_df, 
                                  args.output_dir, particle_type, args.max_cent)
        except Exception as e:
            print(f"Error processing source {source}: {e}")
            continue
    
    print(f"\nAll plots saved to {args.output_dir}")

if __name__ == "__main__":
    main()