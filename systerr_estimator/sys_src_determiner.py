#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import sys
from pathlib import Path


def determine_systematic_sources(df, min_occurrences=2, centrality_col='centrality', max_cent=60, verbose=False):
    """
    Determine which sources are systematic error sources based on two rules:
    1. Rule 1: For each row, if is_delta_over_thr/is_gamma_over_thr is True, 
       mark that row as a systematic error source.
    2. Rule 2: For each (diff_type, diff_bin) combination and each source, 
       if the total count of is_over_thr = True across all centralities >= min_occurrences,
       mark that source in ALL centralities for that combination as a systematic error source.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with Barlow ratio results
    min_occurrences : int
        Minimum total count of over-threshold values across all centralities
        for a source to be considered systematic for a given diff combination
    centrality_col : str
        Name of the column containing centrality information
    max_cent : int
        Maximum centrality value to consider, points with centrality > max_cent will be excluded
    verbose : bool
        Whether to print detailed information
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe with added is_delta_syssrc and is_gamma_syssrc columns
    """
    # Create a copy of the input dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Initialize new columns
    result_df['is_delta_syssrc'] = False
    result_df['is_gamma_syssrc'] = False
    
    # Get unique sources and centrality bins, filtering out high centrality values
    sources = result_df['source'].unique()
    
    # Filter centralities based on max_cent parameter
    centralities = result_df[result_df[centrality_col] <= max_cent][centrality_col].unique()
    
    if verbose:
        print(f"Found {len(sources)} unique sources and {len(centralities)} centrality bins (≤ {max_cent})")
    
    # Check column existence and whether they contain meaningful data
    has_diff_type = 'diff_type' in df.columns and not df['diff_type'].isna().all() and not (df['diff_type'] == '').all()
    has_diff_bin = 'diff_bin' in df.columns and not df['diff_bin'].isna().all() and not (df['diff_bin'] == '').all()
    
    if verbose:
        print(f"Using diff_type for grouping: {has_diff_type}")
        print(f"Using diff_bin for grouping: {has_diff_bin}")
    
    # Get unique diff combinations
    if has_diff_type and has_diff_bin:
        diff_combinations = result_df[['diff_type', 'diff_bin']].drop_duplicates().values.tolist()
    else:
        diff_combinations = [('', '')]  # For Lambda data without diff columns
    
    # Process delta and gamma separately
    # First process delta
    for source in sources:
        source_df = result_df[result_df['source'] == source]
        
        # Apply Rule 1 for delta - mark individual rows
        for centrality in centralities:
            cent_df = source_df[source_df[centrality_col] == centrality]
            
            # Mark individual rows that are over threshold
            for idx, row in cent_df.iterrows():
                if row['is_delta_over_thr'] == True:
                    result_df.loc[idx, 'is_delta_syssrc'] = True
                    
                    if verbose:
                        if has_diff_type and has_diff_bin:
                            print(f"Delta Rule 1: Marked source={source}, cent={centrality}, diff_type={row['diff_type']}, diff_bin={row['diff_bin']}, pair_type={row['pair_type']}")
                        else:
                            print(f"Delta Rule 1: Marked source={source}, cent={centrality}, pair_type={row['pair_type']}")
        
        # Apply Rule 2 for delta - for each diff combination separately
        for diff_combo in diff_combinations:
            diff_type, diff_bin = diff_combo
            
            # Filter data for this source and diff combination
            if has_diff_type and has_diff_bin:
                combo_df = source_df[(source_df['diff_type'] == diff_type) & (source_df['diff_bin'] == diff_bin)]
            else:
                combo_df = source_df
            
            # Count over-threshold occurrences across all centralities for this diff combination
            over_thr_count = combo_df['is_delta_over_thr'].sum()
            
            if over_thr_count >= min_occurrences:
                if verbose:
                    if has_diff_type and has_diff_bin:
                        print(f"Delta Rule 2: Source {source} in diff_type={diff_type}, diff_bin={diff_bin} has {over_thr_count} over-threshold occurrences (≥ {min_occurrences})")
                    else:
                        print(f"Delta Rule 2: Source {source} has {over_thr_count} over-threshold occurrences (≥ {min_occurrences})")
                    print(f"  Marking ALL bins for this source in this diff combination as delta systematic source")
                
                # Mark ALL rows for this source in this diff combination
                if has_diff_type and has_diff_bin:
                    mask = (result_df['source'] == source) & (result_df['diff_type'] == diff_type) & (result_df['diff_bin'] == diff_bin)
                else:
                    mask = result_df['source'] == source
                result_df.loc[mask, 'is_delta_syssrc'] = True
    
    # Now process gamma
    for source in sources:
        source_df = result_df[result_df['source'] == source]
        
        # Apply Rule 1 for gamma - mark individual rows
        for centrality in centralities:
            cent_df = source_df[source_df[centrality_col] == centrality]
            
            # Mark individual rows that are over threshold
            for idx, row in cent_df.iterrows():
                if row['is_gamma_over_thr'] == True:
                    result_df.loc[idx, 'is_gamma_syssrc'] = True
                    
                    if verbose:
                        if has_diff_type and has_diff_bin:
                            print(f"Gamma Rule 1: Marked source={source}, cent={centrality}, diff_type={row['diff_type']}, diff_bin={row['diff_bin']}, pair_type={row['pair_type']}")
                        else:
                            print(f"Gamma Rule 1: Marked source={source}, cent={centrality}, pair_type={row['pair_type']}")
        
        # Apply Rule 2 for gamma - for each diff combination separately
        for diff_combo in diff_combinations:
            diff_type, diff_bin = diff_combo
            
            # Filter data for this source and diff combination
            if has_diff_type and has_diff_bin:
                combo_df = source_df[(source_df['diff_type'] == diff_type) & (source_df['diff_bin'] == diff_bin)]
            else:
                combo_df = source_df
            
            # Count over-threshold occurrences across all centralities for this diff combination
            over_thr_count = combo_df['is_gamma_over_thr'].sum()
            
            if over_thr_count >= min_occurrences:
                if verbose:
                    if has_diff_type and has_diff_bin:
                        print(f"Gamma Rule 2: Source {source} in diff_type={diff_type}, diff_bin={diff_bin} has {over_thr_count} over-threshold occurrences (≥ {min_occurrences})")
                    else:
                        print(f"Gamma Rule 2: Source {source} has {over_thr_count} over-threshold occurrences (≥ {min_occurrences})")
                    print(f"  Marking ALL bins for this source in this diff combination as gamma systematic source")
                
                # Mark ALL rows for this source in this diff combination
                if has_diff_type and has_diff_bin:
                    mask = (result_df['source'] == source) & (result_df['diff_type'] == diff_type) & (result_df['diff_bin'] == diff_bin)
                else:
                    mask = result_df['source'] == source
                result_df.loc[mask, 'is_gamma_syssrc'] = True
    
    # Count systematic sources and bins
    delta_syssrc_rows = result_df[result_df['is_delta_syssrc'] == True]
    gamma_syssrc_rows = result_df[result_df['is_gamma_syssrc'] == True]
    
    delta_syssrc_sources = delta_syssrc_rows['source'].unique()
    gamma_syssrc_sources = gamma_syssrc_rows['source'].unique()
    
    print(f"Identified {len(delta_syssrc_sources)} delta systematic sources with {len(delta_syssrc_rows)} affected bins")
    print(f"Identified {len(gamma_syssrc_sources)} gamma systematic sources with {len(gamma_syssrc_rows)} affected bins")
    
    return result_df


def process_file(input_file, output_file=None, min_occurrences=2, max_cent=60, verbose=False):
    """
    Process a single file and determine systematic sources
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file (from barlow_ratio_calculator.py)
    output_file : str, optional
        Path to the output CSV file. If None, will be derived from input_file
    min_occurrences : int
        Minimum number of occurrences to consider a source as systematic
    max_cent : int
        Maximum centrality value to consider, points with centrality > max_cent will be excluded
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
            'delta', 'delta_err', 'gamma', 'gamma_err',
            'delta_barlow_ratio', 'gamma_barlow_ratio',
            'is_delta_over_thr', 'is_gamma_over_thr'
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
        
        # Determine systematic sources
        result_df = determine_systematic_sources(df, min_occurrences, 'centrality', max_cent, verbose=verbose)
        
        # Determine output file path if not provided
        if output_file is None:
            input_path = Path(input_file)
            output_file = str(input_path.parent / f"syssrc_{input_path.name}")
        
        # Define columns for output
        output_columns = [
            'source', 'centrality', 'diff_type', 'diff_bin', 'pair_type', 
            'delta', 'delta_err', 'gamma', 'gamma_err',
            'delta_diff', 'gamma_diff', 
            'delta_relative_diff', 'gamma_relative_diff',
            'delta_barlow_ratio', 'gamma_barlow_ratio',
            'is_delta_over_thr', 'is_gamma_over_thr',
            'is_delta_syssrc', 'is_gamma_syssrc'
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
        
        # Print summary
        delta_syssrc_rows = result_df[result_df['is_delta_syssrc'] == True]
        gamma_syssrc_rows = result_df[result_df['is_gamma_syssrc'] == True]
        
        delta_syssrc_sources = len(delta_syssrc_rows['source'].unique())
        gamma_syssrc_sources = len(gamma_syssrc_rows['source'].unique())
        
        total_sources = len(result_df['source'].unique())
        
        print(f"Summary:")
        print(f"  Total sources analyzed: {total_sources}")
        print(f"  Sources identified as delta systematic sources: {delta_syssrc_sources} ({delta_syssrc_sources/total_sources*100:.1f}%)")
        print(f"  Sources identified as gamma systematic sources: {gamma_syssrc_sources} ({gamma_syssrc_sources/total_sources*100:.1f}%)")
        print(f"  Total bins marked as delta systematic: {len(delta_syssrc_rows)}")
        print(f"  Total bins marked as gamma systematic: {len(gamma_syssrc_rows)}")
        
        return output_file
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return None


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Determine systematic error sources from Barlow ratio calculations.')
    parser.add_argument('--input', type=str, required=True,
                      help='Input CSV file (from barlow_ratio_calculator.py)')
    parser.add_argument('--output', type=str, default=None,
                      help='Output CSV file (default: syssrc_<input_filename>)')
    parser.add_argument('--min-occurrences', type=int, default=3,
                      help='Minimum total count of over-threshold values across all centralities to mark a source as systematic (default: 3)')
    parser.add_argument('--max-cent', type=int, default=60,
                      help='Maximum centrality value to consider, points with centrality > max_cent will be excluded (default: 60)')
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)
    
    # Process the file
    process_file(args.input, args.output, args.min_occurrences, args.max_cent, args.verbose)


if __name__ == "__main__":
    main()