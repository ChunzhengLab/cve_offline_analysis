#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import sys
from pathlib import Path


def calculate_contribution_fractions(df, group_by_cols=None):
    """
    Calculate contribution fractions for systematic error sources
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with systematic source determinations
    group_by_cols : list
        Columns to group by when calculating contribution fractions
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe with added contribution fraction columns
    """
    # Create a copy of the input dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Initialize new columns
    result_df['delta_contribution_fraction'] = 0.0
    result_df['gamma_contribution_fraction'] = 0.0
    
    # If no group_by columns specified, use all relevant columns
    if group_by_cols is None:
        group_by_cols = ['centrality', 'diff_type', 'diff_bin', 'pair_type']
    
    # Ensure all group_by columns exist in the dataframe
    group_by_cols = [col for col in group_by_cols if col in result_df.columns]
    
    # Group by specified columns
    groups = result_df.groupby(group_by_cols)
    
    # Process each group
    for group_name, group_data in groups:
        # Only consider rows marked as systematic sources
        delta_sources = group_data[group_data['is_delta_syssrc'] == True]
        gamma_sources = group_data[group_data['is_gamma_syssrc'] == True]
        
        # Calculate total squared sum for delta and gamma
        if len(delta_sources) > 0:
            delta_squared_sum = np.sum(delta_sources['delta_diff'] ** 2)
        else:
            delta_squared_sum = 0
            
        if len(gamma_sources) > 0:
            gamma_squared_sum = np.sum(gamma_sources['gamma_diff'] ** 2)
        else:
            gamma_squared_sum = 0
        
        # Calculate contribution fractions for each source in the group
        for idx, row in group_data.iterrows():
            # Delta contribution
            if row['is_delta_syssrc'] and delta_squared_sum > 0:
                delta_contrib = (row['delta_diff'] ** 2) / delta_squared_sum
                result_df.loc[idx, 'delta_contribution_fraction'] = delta_contrib
            
            # Gamma contribution
            if row['is_gamma_syssrc'] and gamma_squared_sum > 0:
                gamma_contrib = (row['gamma_diff'] ** 2) / gamma_squared_sum
                result_df.loc[idx, 'gamma_contribution_fraction'] = gamma_contrib
    
    return result_df


def save_formatted_output(result_df, output_file, verbose=False):
    """
    Save contribution fractions in a formatted text file
    
    Parameters:
    -----------
    result_df : pandas.DataFrame
        Dataframe with contribution fraction results
    output_file : str
        Path to the output text file
    verbose : bool
        Whether to print detailed information
    """
    # Verify contribution fractions sum to approximately 1
    def verify_fractions(data, column, group_cols):
        groups = data[data[column] > 0].groupby(group_cols)
        for name, group in groups:
            sum_val = group[column].sum()
            if abs(sum_val - 1.0) > 0.01:  # Allow for small floating-point errors
                print(f"Warning: {column} sum = {sum_val:.4f} for group {name} (should be 1.0)")
                
    # Verify delta and gamma contribution fractions
    group_cols = ['centrality', 'pair_type']
    if 'diff_type' in result_df.columns:
        group_cols.append('diff_type')
    if 'diff_bin' in result_df.columns:
        group_cols.append('diff_bin')
        
    if verbose:
        verify_fractions(result_df, 'delta_contribution_fraction', group_cols)
        verify_fractions(result_df, 'gamma_contribution_fraction', group_cols)
    # Get unique diff_type and diff_bin combinations
    diff_combinations = []
    
    # Check if diff columns exist
    has_diff_type = 'diff_type' in result_df.columns and not result_df['diff_type'].isna().all() and not (result_df['diff_type'] == '').all()
    has_diff_bin = 'diff_bin' in result_df.columns and not result_df['diff_bin'].isna().all() and not (result_df['diff_bin'] == '').all()
    
    if has_diff_type and has_diff_bin:
        # Get unique combinations of diff_type and diff_bin
        diff_combinations = result_df[['diff_type', 'diff_bin']].drop_duplicates().values.tolist()
    elif has_diff_type:
        # Only diff_type is used
        diff_combinations = [[dt, ''] for dt in result_df['diff_type'].unique() if dt != '']
    elif has_diff_bin:
        # Only diff_bin is used
        diff_combinations = [['', db] for db in result_df['diff_bin'].unique() if db != '']
    else:
        # No diff columns used
        diff_combinations = [['', '']]
    
    # Map centrality values to readable ranges
    centrality_map = {
        5.0: '0-10%',
        15.0: '10-20%',
        25.0: '20-30%',
        35.0: '30-40%',
        45.0: '40-50%',
        55.0: '50-60%',
        65.0: '60-70%',
        75.0: '70-80%',
        85.0: '80-90%'
    }
    
    # Get unique pair_types
    pair_types = sorted(result_df['pair_type'].unique())
    
    # Open output file
    with open(output_file, 'w') as f:
        # Process each diff combination
        for diff_type, diff_bin in diff_combinations:
            # Write header for this diff combination
            f.write(f"[{diff_type}], [{diff_bin}]\n")
            
            # Filter data for this diff combination
            if has_diff_type and has_diff_bin:
                group_data = result_df[(result_df['diff_type'] == diff_type) & (result_df['diff_bin'] == diff_bin)]
            elif has_diff_type:
                group_data = result_df[result_df['diff_type'] == diff_type]
            elif has_diff_bin:
                group_data = result_df[result_df['diff_bin'] == diff_bin]
            else:
                group_data = result_df
            
            # Get unique sources and centralities
            sources = sorted(group_data['source'].unique())
            centralities = sorted(group_data['centrality'].unique())
            
            if 'default' in sources:
                sources.remove('default')
            
            # Process each pair_type separately
            for pair_type in pair_types:
                # Filter data for this pair_type
                pair_data = group_data[group_data['pair_type'] == pair_type]
                
                if len(pair_data) == 0:
                    continue
                
                f.write(f"\n=== Pair Type: {pair_type} ===\n")
                
                # Write column headers
                f.write("Centrality,")
                f.write(",".join([f" {source}" for source in sources]))
                f.write("\n")
                
                # Write delta contribution fractions
                f.write("--- Delta contributions ---\n")
                for cent in centralities:
                    cent_data = pair_data[pair_data['centrality'] == cent]
                    if len(cent_data) == 0:
                        continue
                        
                    cent_str = centrality_map.get(cent, str(cent))
                    f.write(f"{cent_str},")
                    
                    # Write fraction for each source
                    fractions = []
                    for source in sources:
                        source_data = cent_data[cent_data['source'] == source]
                        if len(source_data) > 0 and source_data['is_delta_syssrc'].any():
                            frac = source_data['delta_contribution_fraction'].values[0]
                            fractions.append(f" {frac:.4f}")
                        else:
                            fractions.append(" 0.0000")
                    
                    f.write(",".join(fractions))
                    f.write("\n")
                
                # Write gamma contribution fractions
                f.write("\n--- Gamma contributions ---\n")
                for cent in centralities:
                    cent_data = pair_data[pair_data['centrality'] == cent]
                    if len(cent_data) == 0:
                        continue
                        
                    cent_str = centrality_map.get(cent, str(cent))
                    f.write(f"{cent_str},")
                    
                    # Write fraction for each source
                    fractions = []
                    for source in sources:
                        source_data = cent_data[cent_data['source'] == source]
                        if len(source_data) > 0 and source_data['is_gamma_syssrc'].any():
                            frac = source_data['gamma_contribution_fraction'].values[0]
                            fractions.append(f" {frac:.4f}")
                        else:
                            fractions.append(" 0.0000")
                    
                    f.write(",".join(fractions))
                    f.write("\n")
            
            # Add a blank line between diff combinations
            f.write("\n\n")
    
    # Verify and report total contribution sums
    if verbose:
        # Re-verify after writing to ensure consistency
        delta_sources = result_df[result_df['is_delta_syssrc'] == True]
        gamma_sources = result_df[result_df['is_gamma_syssrc'] == True]
        
        print("Contribution fraction verification:")
        for cols in group_cols:
            delta_totals = delta_sources.groupby(cols)['delta_contribution_fraction'].sum()
            gamma_totals = gamma_sources.groupby(cols)['gamma_contribution_fraction'].sum()
            
            if not delta_totals.empty:
                print(f"  Delta contribution sums by {cols}: Min={delta_totals.min():.4f}, Max={delta_totals.max():.4f}")
            if not gamma_totals.empty:
                print(f"  Gamma contribution sums by {cols}: Min={gamma_totals.min():.4f}, Max={gamma_totals.max():.4f}")
    
    print(f"Successfully created formatted output file: {output_file}")

def process_file(input_file, output_file=None, formatted_output=None, verbose=False):
    """
    Process a single file and calculate contribution fractions
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file (from sys_src_determiner.py)
    output_file : str, optional
        Path to the output CSV file. If None, will be derived from input_file
    formatted_output : str, optional
        Path to the formatted text output file. If None, will be derived from input_file
    verbose : bool
        Whether to print detailed information
        
    Returns:
    --------
    str
        Path to the output file
    """
    try:
        print(f"Reading input file: {input_file}")
        # Read the input file
        df = pd.read_csv(input_file)
        
        # Verify the file has the required columns
        required_cols = [
            'source', 'centrality', 'pair_type', 
            'delta_diff', 'gamma_diff',
            'is_delta_syssrc', 'is_gamma_syssrc'
        ]
        
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: Input file {input_file} is missing required columns: {', '.join(missing_cols)}")
            sys.exit(1)
        
        # Skip default data if present
        if 'source' in df.columns and 'default' in df['source'].values:
            print("Removing default data from analysis")
            df = df[df['source'] != 'default']
        
        if len(df) == 0:
            print(f"Warning: No data found in {input_file} after filtering")
            return None
        
        # Determine group by columns based on available data
        group_by_cols = ['centrality', 'pair_type']
        if 'diff_type' in df.columns and not df['diff_type'].isna().all() and not (df['diff_type'] == '').all():
            group_by_cols.append('diff_type')
        if 'diff_bin' in df.columns and not df['diff_bin'].isna().all() and not (df['diff_bin'] == '').all():
            group_by_cols.append('diff_bin')
            
        print(f"Grouping by columns: {group_by_cols}")
        
        # Calculate contribution fractions
        result_df = calculate_contribution_fractions(df, group_by_cols)
        
        # Determine output file paths if not provided
        input_path = Path(input_file)
        if output_file is None:
            output_file = str(input_path.parent / f"contrib_{input_path.name}")
        if formatted_output is None:
            formatted_output = str(input_path.parent / f"contrib_formatted_{input_path.stem}.txt")
        
        # Define columns for output
        output_columns = [
            'source', 'centrality', 'diff_type', 'diff_bin', 'pair_type', 
            'delta', 'delta_err', 'gamma', 'gamma_err',
            'delta_diff', 'gamma_diff', 
            'delta_relative_diff', 'gamma_relative_diff',
            'delta_barlow_ratio', 'gamma_barlow_ratio',
            'is_delta_over_thr', 'is_gamma_over_thr',
            'is_delta_syssrc', 'is_gamma_syssrc',
            'delta_contribution_fraction', 'gamma_contribution_fraction'
        ]
        
        # Ensure all columns exist
        for col in output_columns:
            if col not in result_df.columns and col in ['diff_type', 'diff_bin']:
                result_df[col] = ''
            elif col not in result_df.columns:
                result_df[col] = np.nan
        
        # Save to CSV
        result_df[output_columns].to_csv(output_file, index=False)
        print(f"Successfully created {output_file}")
        
        # Save formatted output
        save_formatted_output(result_df, formatted_output, verbose)
        
        # Print summary
        delta_syssrc_count = result_df['is_delta_syssrc'].sum()
        gamma_syssrc_count = result_df['is_gamma_syssrc'].sum()
        
        print(f"Summary:")
        print(f"  Total sources analyzed: {len(result_df['source'].unique())}")
        print(f"  Sources identified as delta systematic sources: {delta_syssrc_count}")
        print(f"  Sources identified as gamma systematic sources: {gamma_syssrc_count}")
        
        # Calculate statistics for contribution fractions
        delta_contrib = result_df[result_df['is_delta_syssrc']]['delta_contribution_fraction']
        gamma_contrib = result_df[result_df['is_gamma_syssrc']]['gamma_contribution_fraction']
        
        if len(delta_contrib) > 0:
            print(f"  Delta contribution statistics:")
            print(f"    Min: {delta_contrib.min():.4f}, Max: {delta_contrib.max():.4f}, Mean: {delta_contrib.mean():.4f}")
            
        if len(gamma_contrib) > 0:
            print(f"  Gamma contribution statistics:")
            print(f"    Min: {gamma_contrib.min():.4f}, Max: {gamma_contrib.max():.4f}, Mean: {gamma_contrib.mean():.4f}")
        
        return output_file
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate contribution fractions for systematic error sources.')
    parser.add_argument('--input', type=str, required=True,
                      help='Input CSV file (from sys_src_determiner.py)')
    parser.add_argument('--output', type=str, default=None,
                      help='Output CSV file (default: contrib_<input_filename>)')
    parser.add_argument('--formatted-output', type=str, default=None,
                      help='Formatted text output file (default: contrib_formatted_<input_filename>.txt)')
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)
    
    # Process the file
    process_file(args.input, args.output, args.formatted_output, args.verbose)


if __name__ == "__main__":
    main()