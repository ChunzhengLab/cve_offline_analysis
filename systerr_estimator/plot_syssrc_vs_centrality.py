#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path

def centre_of(label):
    """Convert centrality label to numerical center value"""
    if isinstance(label, (int, float)):
        return float(label)
    m = re.match(r"(\d*\.?\d+)\D+(\d*\.?\d+)", str(label))
    if m:
        a, b = map(float, m.groups())
        return 0.5 * (a + b)
    try:
        return float(label)
    except Exception:
        return np.nan

def get_pair_label(particle_name):
    """Get pair label based on particle type"""
    if particle_name == 'Proton':
        return 'Lambda-Proton Pair'
    elif particle_name == 'Lambda':
        return 'Lambda-Lambda Pair'
    elif particle_name == 'Hadron':
        return 'Lambda-Hadron'
    else:
        return f'{particle_name} Pair'

def get_diff_label(diff_type, diff_bin):
    """Get differential label based on diff_type and diff_bin"""
    if diff_type == 'Intg':
        return 'Intg'
    elif diff_type == 'SPt':
        if diff_bin == '0.5':
            return '1.0 < pT < 3.0 GeV/c'
        elif diff_bin == '1.5':
            return '3.0 < pT < 5.0 GeV/c'
        elif diff_bin == '2.5':
            return '5.0 < pT < 8.0 GeV/c'
        else:
            return f'SPt {diff_bin}'
    elif diff_type == 'DEta':
        if diff_bin == '0.5':
            return '|Δη| < 0.6'
        elif diff_bin == '1.5':
            return '|Δη| > 0.6'
        else:
            return f'DEta {diff_bin}'
    else:
        return f'{diff_type} {diff_bin}' if diff_bin else str(diff_type)

def plot_syssrc_vs_centrality(input_file, output_dir=None, 
                             diff_type=None, diff_bin=None, 
                             min_points=2, figsize=(20, 6)):
    """
    Plot |var-def| vs centrality for systematic error sources
    Creates separate plots for SS, OS, Del with delta on left and gamma on right
    
    Parameters:
    -----------
    input_file : str
        Path to syssrc_barlow_finalise_sys_<particle>.csv file
    output_dir : str
        Output directory for plots
    diff_type : str
        Filter by specific diff_type (optional)  
    diff_bin : str
        Filter by specific diff_bin (optional)
    min_points : int
        Minimum number of points required to draw a line for a source
    figsize : tuple
        Figure size for subplots
    """
    
    # Read the data
    try:
        df = pd.read_csv(input_file)
        print(f"Successfully read {len(df)} rows from {input_file}")
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return
    
    # Check required columns
    required_cols = ['centrality', 'source', 'delta_diff', 'gamma_diff', 
                    'is_delta_syssrc', 'is_gamma_syssrc', 'pair_type']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        return
    
    # Apply filters
    filtered_df = df.copy()
    
    if diff_type is not None:
        filtered_df = filtered_df[filtered_df['diff_type'] == diff_type]
        print(f"Filtered by diff_type={diff_type}: {len(filtered_df)} rows")
        
    if diff_bin is not None:
        filtered_df = filtered_df[filtered_df['diff_bin'] == diff_bin]
        print(f"Filtered by diff_bin={diff_bin}: {len(filtered_df)} rows")
    
    if len(filtered_df) == 0:
        print("No data remaining after filtering")
        return
    
    # Convert centrality to numerical values and filter centrality >= 60%
    filtered_df['cent_numeric'] = filtered_df['centrality'].apply(centre_of)
    filtered_df = filtered_df.dropna(subset=['cent_numeric'])
    
    # Remove data points with centrality >= 60%
    filtered_df = filtered_df[filtered_df['cent_numeric'] < 60]
    print(f"After filtering centrality < 60%: {len(filtered_df)} rows")
    
    # Determine particle name for titles
    particle_name = Path(input_file).stem.replace('syssrc_barlow_finalise_sys_', '')
    
    # Get labels for titles
    pair_label = get_pair_label(particle_name)
    diff_label = get_diff_label(diff_type, diff_bin) if diff_type else "All"
    
    # Set output directory
    if output_dir is None:
        output_dir = Path(input_file).parent / "plots_syssrc_vs_cent"
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Get all unique combinations of diff_type and diff_bin
    if diff_type is None and diff_bin is None:
        # Get all available combinations
        combinations = []
        for _, row in filtered_df[['diff_type', 'diff_bin']].drop_duplicates().iterrows():
            dt = row['diff_type'] if pd.notna(row['diff_type']) else None
            db = row['diff_bin'] if pd.notna(row['diff_bin']) else None
            combinations.append((dt, db))
        
        if not combinations:
            combinations = [(None, None)]
    else:
        combinations = [(diff_type, diff_bin)]
    
    print(f"Will create plots for combinations: {combinations}")
    
    # Create plots for each combination
    for diff_type_val, diff_bin_val in combinations:
        # Filter data for this combination
        combo_df = filtered_df.copy()
        if diff_type_val is not None:
            combo_df = combo_df[combo_df['diff_type'] == diff_type_val]
        if diff_bin_val is not None:
            combo_df = combo_df[combo_df['diff_bin'] == diff_bin_val]
        
        if len(combo_df) == 0:
            continue
            
        combo_label = get_diff_label(diff_type_val, diff_bin_val)
        print(f"\nProcessing combination: {combo_label}")
        
        # Get available pair types for this combination
        available_pair_types = ['SS', 'OS', 'Del']
        actual_pair_types = [pt for pt in available_pair_types if pt in combo_df['pair_type'].values]
        
        if not actual_pair_types:
            print(f"No SS, OS, or Del pair types found for {combo_label}")
            continue
        
        # Create plots for each pair type
        for pair_type in actual_pair_types:
            pair_df = combo_df[combo_df['pair_type'] == pair_type].copy()
            
            if len(pair_df) == 0:
                continue
                
            print(f"  Processing pair_type: {pair_type}")
            
            # Create subplot figure with delta on left, gamma on right
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
            
            # Function to plot one variable
            def plot_one_variable(ax, var_name, var_full_name, syssrc_col, diff_col, color_map):
                var_df = pair_df[pair_df[syssrc_col] == True].copy()
                
                if len(var_df) == 0:
                    ax.text(0.5, 0.5, f'No {var_name} systematic\nerror sources', 
                           transform=ax.transAxes, ha='center', va='center', fontsize=12)
                    ax.set_xlabel('Centrality (%)')
                    ax.set_ylabel('|var - def|')
                    
                    # Add annotation even when no data
                    info_text = f"{pair_label}\n{pair_type}\n{combo_label}\n{var_full_name}"
                    ax.text(0.05, 0.95, info_text, transform=ax.transAxes, 
                           verticalalignment='top', fontsize=10,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))
                    return 0, [], {}
                
                # Calculate absolute differences
                var_df['abs_diff'] = np.abs(var_df[diff_col])
                
                # Get unique sources
                sources = var_df['source'].unique()
                
                plotted_sources = []
                source_colors = {}  # Store the color mapping
                
                for source in sources:
                    source_data = var_df[var_df['source'] == source].copy()
                    source_data = source_data.sort_values('cent_numeric')
                    
                    # Get color from global mapping
                    source_color = color_map[source]
                    source_colors[source] = source_color
                    
                    # Plot horizontal lines for each centrality bin
                    for _, row in source_data.iterrows():
                        cent = row['cent_numeric']
                        value = row['abs_diff']
                        
                        # Draw horizontal line representing the centrality bin width
                        # Assume centrality bins are about 10% wide
                        bin_width = 5  # Half-width of the horizontal line
                        ax.plot([cent - bin_width, cent + bin_width], [value, value], 
                               color=source_color, linewidth=3, alpha=0.8)
                    
                    # Add source to legend (no connecting lines)
                    if len(source_data) > 0:
                        ax.plot([], [], color=source_color, label=source, linewidth=3)
                        plotted_sources.append(source)
                
                # Calculate total systematic error (quadrature sum)
                centralities = var_df['cent_numeric'].unique()
                centralities.sort()
                
                total_syst_err = []
                for cent in centralities:
                    cent_data = var_df[var_df['cent_numeric'] == cent]
                    # Sum of squares of absolute differences
                    sum_sq = (cent_data['abs_diff'] ** 2).sum()
                    total_syst_err.append(np.sqrt(sum_sq))
                
                # Plot total systematic error as black horizontal lines (no connecting line)
                if len(total_syst_err) > 0:
                    bin_width = 5  # Same as individual sources
                    
                    # Draw horizontal lines for each centrality bin
                    for cent, total_err in zip(centralities, total_syst_err):
                        ax.plot([cent - bin_width, cent + bin_width], [total_err, total_err], 
                               'k-', linewidth=4, alpha=0.9)
                    
                    # Add to legend (no connecting line)
                    ax.plot([], [], 'k-', linewidth=4, label='Total Systematic Error')
                
                ax.set_xlabel('Centrality (%)')
                ax.set_ylabel('|var - def|')
                ax.grid(True, alpha=0.3)
                
                # Add information annotation
                info_text = f"{pair_label}\n{pair_type}\n{combo_label}\n{var_full_name}"
                ax.text(0.05, 0.95, info_text, transform=ax.transAxes, 
                       verticalalignment='top', fontsize=10,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))
                
                return len(plotted_sources), plotted_sources, source_colors
            
            # Get all unique sources first to ensure consistent coloring
            delta_df = pair_df[pair_df['is_delta_syssrc'] == True]
            gamma_df = pair_df[pair_df['is_gamma_syssrc'] == True]
            all_unique_sources = list(set(delta_df['source'].unique().tolist() + gamma_df['source'].unique().tolist()))
            global_colors = plt.cm.tab20(np.linspace(0, 1, len(all_unique_sources)))
            global_color_map = {source: global_colors[i] for i, source in enumerate(all_unique_sources)}
            
            # Plot delta (left) and gamma (right) with consistent colors
            delta_count, delta_sources, delta_colors = plot_one_variable(ax1, 'δ', 'Delta', 'is_delta_syssrc', 'delta_diff', global_color_map)
            gamma_count, gamma_sources, gamma_colors = plot_one_variable(ax2, 'γ', 'Gamma', 'is_gamma_syssrc', 'gamma_diff', global_color_map)
            
            # Create shared legend on the right side
            all_sources = list(set(delta_sources + gamma_sources))
            if all_sources:
                # Create legend with all sources using global color mapping
                legend_elements = []
                for source in all_sources:
                    legend_elements.append(plt.Line2D([0], [0], color=global_color_map[source], lw=3, label=source))
                
                # Add total systematic error line
                legend_elements.append(plt.Line2D([0], [0], color='black', lw=4, label='Total Systematic Error'))
                
                # Place legend outside the plot area on the right
                fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.83, 0.5), 
                          fontsize=8, ncol=1)
            
            # Set main title (simplified since info is now in each subplot)
            main_title = f'{pair_label} - {pair_type} - {combo_label}'
            fig.suptitle(main_title, fontsize=14, y=0.95)
            
            # Adjust layout to make room for right legend
            plt.tight_layout()
            plt.subplots_adjust(top=0.85, right=0.82)
            
            # Save plot
            filter_str = []
            if diff_type_val: filter_str.append(str(diff_type_val))
            if diff_bin_val: filter_str.append(str(diff_bin_val))
            filter_suffix = "_".join(filter_str) if filter_str else "all"
            
            filename = f"syssrc_vs_centrality_{particle_name}_{pair_type}_{filter_suffix}.pdf"
            filepath = output_dir / filename
            plt.savefig(filepath, bbox_inches='tight', dpi=300)
            plt.close()
            
            print(f"    Saved plot: {filepath}")
            print(f"    Delta sources: {delta_count}, Gamma sources: {gamma_count}")
            print(f"    Total {pair_type} data points: {len(pair_df)}")

def main():
    parser = argparse.ArgumentParser(description='Plot systematic error sources vs centrality')
    parser.add_argument('--input', required=True,
                       help='Input CSV file (syssrc_barlow_finalise_sys_<particle>.csv)')
    parser.add_argument('--output-dir', 
                       help='Output directory for plots (default: plots_syssrc_vs_cent/)')
    parser.add_argument('--diff-type',
                       help='Filter by specific diff_type')
    parser.add_argument('--diff-bin',
                       help='Filter by specific diff_bin')
    parser.add_argument('--min-points', type=int, default=2,
                       help='Minimum points to draw line for a source (default: 2)')
    parser.add_argument('--figsize', nargs=2, type=float, default=[20, 6],
                       help='Figure size in inches (default: 20 6)')
    parser.add_argument('--list-filters', action='store_true',
                       help='List available filter values and exit')
    
    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        return
    
    # List available filter values if requested
    if args.list_filters:
        try:
            df = pd.read_csv(args.input)
            print("Available filter values:")
            for col in ['pair_type', 'diff_type', 'diff_bin']:
                if col in df.columns:
                    unique_vals = df[col].dropna().unique()
                    print(f"  {col}: {list(unique_vals)}")
                else:
                    print(f"  {col}: column not found")
            return
        except Exception as e:
            print(f"Error reading file for filter listing: {e}")
            return
    
    # Create the plots
    plot_syssrc_vs_centrality(
        input_file=args.input,
        output_dir=args.output_dir,
        diff_type=args.diff_type,
        diff_bin=args.diff_bin,
        min_points=args.min_points,
        figsize=tuple(args.figsize)
    )

if __name__ == "__main__":
    main()