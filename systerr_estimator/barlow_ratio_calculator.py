#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
import sys
from pathlib import Path


def calculate_barlow_ratios(df, default_df, threshold=1.0, input_file=None):
    """
    Calculate Barlow ratios for each row in the dataframe by comparing with default values
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with source data
    default_df : pandas.DataFrame
        Default dataframe to compare against
    threshold : float
        Threshold value for determining if a Barlow ratio is significant
    input_file : str, optional
        Path to the input file, used to determine special handling for diff columns
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe with added Barlow ratio columns
    """
    # Create a copy of the input dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Ensure all required columns exist in default_df
    required_cols = ['centrality', 'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']
    for col in required_cols:
        if col not in default_df.columns:
            print(f"Error: Default data is missing required column: {col}")
            sys.exit(1)
    
    # Initialize new columns
    result_df['delta_diff'] = np.nan
    result_df['gamma_diff'] = np.nan
    result_df['delta_relative_diff'] = np.nan
    result_df['gamma_relative_diff'] = np.nan
    result_df['delta_barlow_ratio'] = np.nan
    result_df['gamma_barlow_ratio'] = np.nan
    result_df['is_delta_over_thr'] = False
    result_df['is_gamma_over_thr'] = False
    
    # Determine folder name from input file
    folder_name = 'unknown'
    if input_file:
        folder_name = Path(input_file).stem.replace('finalise_sys_', '')
        print(f"Processing data from folder: {folder_name}")
    
    # Iterate through each row in the dataframe
    for idx, row in result_df.iterrows():
        # Find matching default row based on centrality, diff_type, diff_bin, and pair_type
        try:
            match_conditions = True
            for col in ['centrality', 'pair_type']:
                match_conditions = match_conditions & (default_df[col] == row[col])
            
            # Handle diff columns special case
            # For Lambda data or when diff columns are empty, don't use them in matching
            
            if folder_name == 'Lambda':
                # For Lambda, ignore diff columns in matching
                pass
            elif ('diff_type' in row and (row['diff_type'] == '' or pd.isna(row['diff_type'])) and 
                  'diff_bin' in row and (row['diff_bin'] == '' or pd.isna(row['diff_bin']))):
                # Both diff columns are empty or NaN, don't use in matching
                pass
            else:
                # Use diff columns in matching if they exist and have values
                if 'diff_type' in default_df.columns and 'diff_type' in row:
                    if not pd.isna(row['diff_type']) and row['diff_type'] != '':
                        match_conditions = match_conditions & (default_df['diff_type'] == row['diff_type'])
                    else:
                        # If source has empty diff_type, match against empty diff_type in default
                        match_conditions = match_conditions & ((default_df['diff_type'] == '') | pd.isna(default_df['diff_type']))
                
                if 'diff_bin' in default_df.columns and 'diff_bin' in row:
                    if not pd.isna(row['diff_bin']) and row['diff_bin'] != '':
                        match_conditions = match_conditions & (default_df['diff_bin'] == row['diff_bin'])
                    else:
                        # If source has empty diff_bin, match against empty diff_bin in default
                        match_conditions = match_conditions & ((default_df['diff_bin'] == '') | pd.isna(default_df['diff_bin']))
                
            default_row = default_df[match_conditions]
            
            if len(default_row) == 0:
                print(f"Warning: No matching default row found for {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}")
                print(f"Available default data contains centralities: {default_df['centrality'].unique()}")
                print(f"Available default data contains pair_types: {default_df['pair_type'].unique()}")
            
                # Print more diagnostic info to understand the matching issue
                if 'diff_type' in row and 'diff_bin' in row:
                    print(f"Source diff columns: diff_type='{row['diff_type']}' (is NaN: {pd.isna(row['diff_type'])}), diff_bin='{row['diff_bin']}' (is NaN: {pd.isna(row['diff_bin'])})")
        
                if 'diff_type' in default_df.columns and 'diff_bin' in default_df.columns:
                    diff_type_vals = default_df['diff_type'].unique()
                    diff_bin_vals = default_df['diff_bin'].unique()
                    print(f"Default diff_type values: {diff_type_vals}")
                    print(f"Default diff_bin values: {diff_bin_vals}")
                    print(f"NaN in default diff_type: {pd.isna(default_df['diff_type']).any()}")
                    print(f"NaN in default diff_bin: {pd.isna(default_df['diff_bin']).any()}")
                
                # Try a simpler match for debugging - just match on centrality and pair_type
                simple_match = default_df[(default_df['centrality'] == row['centrality']) & 
                                         (default_df['pair_type'] == row['pair_type'])]
                if len(simple_match) > 0:
                    print(f"Found {len(simple_match)} matches when ignoring diff columns - this suggests diff columns are causing the issue")
                else:
                    print(f"No matches found even when ignoring diff columns - check centrality and pair_type values carefully")
                continue
        except Exception as e:
            print(f"Error while matching default data: {e}")
            continue
        
        if len(default_row) > 1:
            print(f"Warning: Multiple matching default rows found for {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}. Using the first match.")
            
        default_row = default_row.iloc[0]
        
        # Calculate differences
        try:
            # Convert values to float to ensure proper numerical operations
            source_delta = float(row['delta'])
            source_gamma = float(row['gamma'])
            default_delta = float(default_row['delta'])
            default_gamma = float(default_row['gamma'])
            
            delta_diff = source_delta - default_delta
            gamma_diff = source_gamma - default_gamma
            
            # Print diagnostic info to help debug
            if np.isnan(delta_diff) or np.isnan(gamma_diff):
                print(f"Diagnostic: NaN detected in diff calculation for {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}")
                print(f"Source values: delta={source_delta}, gamma={source_gamma}")
                print(f"Default values: delta={default_delta}, gamma={default_gamma}")
        except Exception as e:
            print(f"Error calculating differences: {e}")
            print(f"Source values: delta={row['delta']} (type: {type(row['delta'])}), gamma={row['gamma']} (type: {type(row['gamma'])})")
            print(f"Default values: delta={default_row['delta']} (type: {type(default_row['delta'])}), gamma={default_row['gamma']} (type: {type(default_row['gamma'])})")
            continue
        
        # Calculate relative differences
        delta_relative_diff = delta_diff / default_row['delta'] if default_row['delta'] != 0 else np.nan
        gamma_relative_diff = gamma_diff / default_row['gamma'] if default_row['gamma'] != 0 else np.nan
        
        # Calculate Barlow ratios - using the modified formula:
        # |diff| / sqrt(|err^2[source] - err^2[default]|)
        # Note: We use absolute value in the denominator to handle cases where source error < default error
        try:
            delta_err_diff_sq = row['delta_err']**2 - default_row['delta_err']**2
            gamma_err_diff_sq = row['gamma_err']**2 - default_row['gamma_err']**2
            
            # For cases where source error is smaller than default error, we use absolute value of the difference
            if delta_err_diff_sq <= 0:
                print(f"Warning: Source error is smaller than default error for delta in {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}")
                # Taking absolute value of the error difference square
                delta_err_diff = np.sqrt(abs(delta_err_diff_sq))
                delta_barlow_ratio = abs(delta_diff) / delta_err_diff
            else:
                delta_err_diff = np.sqrt(delta_err_diff_sq)
                delta_barlow_ratio = abs(delta_diff) / delta_err_diff
                
            if gamma_err_diff_sq <= 0:
                print(f"Warning: Source error is smaller than default error for gamma in {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}")
                # Taking absolute value of the error difference square
                gamma_err_diff = np.sqrt(abs(gamma_err_diff_sq))
                gamma_barlow_ratio = abs(gamma_diff) / gamma_err_diff
            else:
                gamma_err_diff = np.sqrt(gamma_err_diff_sq)
                gamma_barlow_ratio = abs(gamma_diff) / gamma_err_diff
        except Exception as e:
            print(f"Error calculating Barlow ratios for {row['source']}, centrality={row['centrality']}, pair_type={row['pair_type']}: {e}")
            delta_barlow_ratio = np.nan
            gamma_barlow_ratio = np.nan
        
        # Store calculated values
        result_df.at[idx, 'delta_diff'] = delta_diff
        result_df.at[idx, 'gamma_diff'] = gamma_diff
        result_df.at[idx, 'delta_relative_diff'] = delta_relative_diff
        result_df.at[idx, 'gamma_relative_diff'] = gamma_relative_diff
        result_df.at[idx, 'delta_barlow_ratio'] = delta_barlow_ratio
        result_df.at[idx, 'is_delta_over_thr'] = False if np.isnan(delta_barlow_ratio) else (delta_barlow_ratio > threshold)
        result_df.at[idx, 'gamma_barlow_ratio'] = gamma_barlow_ratio
        result_df.at[idx, 'is_gamma_over_thr'] = False if np.isnan(gamma_barlow_ratio) else (gamma_barlow_ratio > threshold)
    
    return result_df


def extract_default_data(input_file, strip_diff_columns=False):
    """
    Extract default data from the input file
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file
    strip_diff_columns : bool
        Whether to set diff_type and diff_bin columns to empty strings
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe containing only the default data
    """
    """
    Extract default data from the input file
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file
        
    Returns:
    --------
    pandas.DataFrame
        Dataframe containing only the default data
    """
    print(f"Extracting default data from {input_file}")
    try:
        df = pd.read_csv(input_file)
        print(f"Successfully read {len(df)} rows from input file")
        
        # Check if 'source' column exists
        if 'source' not in df.columns:
            print(f"Error: Input file {input_file} doesn't have a 'source' column")
            print(f"Available columns: {df.columns.tolist()}")
            sys.exit(1)
        
        # Find rows with source='default'
        default_df = df[df['source'] == 'default']
        print(f"Found {len(default_df)} rows with source='default'")
        
        if len(default_df) == 0:
            print("No default data found in input file. Searching for separate default file...")
            # Try to look for default file directly if not found in the combined file
            input_path = Path(input_file)
            folder_name = input_path.stem.replace('finalise_sys_', '')
            print(f"Detected folder name: {folder_name}")
            
            # Try both possible default file patterns
            possible_default_files = []
            if folder_name in ['Lambda', 'Hadron']:
                possible_default_files.append(input_path.parent / f"finalise_default.csv")
            elif folder_name == 'Proton':
                possible_default_files.append(input_path.parent / f"finalise_feeddown_dispose_default.csv")
            else:
                # Try both patterns if folder type can't be determined
                possible_default_files = [
                    input_path.parent / f"finalise_default.csv",
                    input_path.parent / f"finalise_feeddown_dispose_default.csv"
                ]
                print(f"Could not determine folder type from '{folder_name}', will try multiple default file patterns")
            
            default_file = None
            for possible_file in possible_default_files:
                if possible_file.exists():
                    default_file = possible_file
                    print(f"Found separate default file: {default_file}")
                    break
            
            if default_file is not None:
                try:
                    default_df = pd.read_csv(default_file)
                    print(f"Successfully read {len(default_df)} rows from default file")
                    # Add source column if it doesn't exist
                    if 'source' not in default_df.columns:
                        default_df['source'] = 'default'
                        print("Added 'source' column to default data")
                    
                    # Check for required columns
                    missing_cols = []
                    for col in ['centrality', 'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']:
                        if col not in default_df.columns:
                            missing_cols.append(col)
                        
                    if missing_cols:
                        print(f"Error: Default file {default_file} is missing required columns: {', '.join(missing_cols)}")
                        print(f"Available columns: {default_df.columns.tolist()}")
                        sys.exit(1)
                            
                    # Add diff columns if they don't exist
                    if 'diff_type' not in default_df.columns:
                        default_df['diff_type'] = ''
                        print("Added 'diff_type' column to default data")
                    if 'diff_bin' not in default_df.columns:
                        default_df['diff_bin'] = ''
                        print("Added 'diff_bin' column to default data")
                            
                    print(f"Successfully loaded default data from separate file with {len(default_df)} rows")
                except Exception as e:
                    print(f"Error reading default file {default_file}: {e}")
                    sys.exit(1)
            else:
                print(f"Error: No default data found in {input_file} and no suitable default file found in {input_path.parent}")
                print(f"Tried: {[str(f) for f in possible_default_files]}")
                sys.exit(1)
        
        # Print summary of default data
        print(f"Default data summary:")
        print(f"  Centrality values: {default_df['centrality'].unique()}")
        print(f"  Pair types: {default_df['pair_type'].unique()}")
        
        # For Lambda, consider removing diff columns from default data to simplify matching
        folder_name = Path(input_file).stem.replace('finalise_sys_', '')
        if folder_name == 'Lambda' or strip_diff_columns:
            if 'diff_type' in default_df.columns:
                default_df.loc[:, 'diff_type'] = default_df['diff_type'].fillna('')
                print("Set all 'diff_type' values to empty strings for Lambda data")
            if 'diff_bin' in default_df.columns:
                default_df.loc[:, 'diff_bin'] = default_df['diff_bin'].fillna('')
                print("Set all 'diff_bin' values to empty strings for Lambda data")
            
        return default_df
        
    except Exception as e:
        print(f"Unexpected error extracting default data: {e}")
        sys.exit(1)


def process_file(input_file, output_file=None, threshold=1.0, verbose=False):
    """
    Process a single file and calculate Barlow ratios
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file
    output_file : str, optional
        Path to the output CSV file. If None, will be derived from input_file
    threshold : float
        Threshold value for determining if a Barlow ratio is significant
    verbose : bool
        Whether to print detailed information during processing
        
    Returns:
    --------
    str
        Path to the output file
    """
    try:
        print(f"Reading input file: {input_file}")
        # Read the input file
        try:
            df = pd.read_csv(input_file)
            print(f"Successfully read {len(df)} rows from input file")
            if verbose:
                print(f"First few rows of input data:")
                print(df.head())
        except Exception as e:
            print(f"Error reading input file {input_file}: {e}")
            sys.exit(1)
        
        # Verify the file has the required columns
        required_cols = ['source', 'centrality', 'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: Input file {input_file} is missing required columns: {', '.join(missing_cols)}")
            sys.exit(1)
                
        # Convert NaN values in diff columns to empty strings for consistency
        if 'diff_type' in df.columns:
            df['diff_type'] = df['diff_type'].fillna('')
        if 'diff_bin' in df.columns:
            df['diff_bin'] = df['diff_bin'].fillna('')
            
        # Extract default data (either from this file or from a separate default file)
        folder_name = Path(input_file).stem.replace('finalise_sys_', '')
        strip_diff = (folder_name == 'Lambda')
        default_df = extract_default_data(input_file, strip_diff_columns=strip_diff)
        print(f"Found default data with {len(default_df)} rows")
        
        if verbose:
            print("Default data preview:")
            print(default_df.head())
        
        # Remove default data from the main dataframe to avoid comparing default with itself
        df = df[df['source'] != 'default']
        
        if len(df) == 0:
            print(f"Warning: No non-default data found in {input_file}")
            return None
        
        # Calculate Barlow ratios (using absolute value for error differences when needed)
        print(f"Calculating Barlow ratios for file: {input_file}")
        result_df = calculate_barlow_ratios(df, default_df, threshold, input_file)
        
        # Check if the calculation was successful
        if result_df is None:
            print("Error: Barlow ratio calculation failed")
            return None
            
        # Verify results to ensure calculations were performed
        empty_calcs = True
        for col in ['delta_diff', 'gamma_diff', 'delta_barlow_ratio', 'gamma_barlow_ratio']:
            if col in result_df.columns and not result_df[col].isnull().all():
                empty_calcs = False
                break
                
        if empty_calcs:
            print("Warning: All calculated values are empty/NaN. This suggests a problem with matching default data.")
            if verbose:
                print("Result preview (with empty calculations):")
                print(result_df.head())
        
        # Determine output file path if not provided
        if output_file is None:
            input_path = Path(input_file)
            output_file = str(input_path.parent / f"barlow_{input_path.name}")
        
        # Define columns for output
        output_columns = [
            'source', 'centrality', 'diff_type', 'diff_bin', 'pair_type', 
            'delta', 'delta_err', 'gamma', 'gamma_err',
            'delta_diff', 'gamma_diff', 
            'delta_relative_diff', 'gamma_relative_diff',
            'delta_barlow_ratio', 'gamma_barlow_ratio',
            'is_delta_over_thr', 'is_gamma_over_thr'
        ]
        
        # Ensure all columns exist
        for col in output_columns:
            if col not in result_df.columns:
                result_df[col] = np.nan
        
        # Save to CSV
        result_df[output_columns].to_csv(output_file, index=False)
        print(f"Successfully created {output_file}")
        
        # Print summary
        delta_over = result_df['is_delta_over_thr'].sum()
        gamma_over = result_df['is_gamma_over_thr'].sum()
        total_rows = len(result_df)
        
        print(f"Summary:")
        print(f"  Total rows processed: {total_rows}")
        print(f"  Rows with delta Barlow ratio > {threshold}: {delta_over} ({delta_over/total_rows*100:.1f}%)")
        print(f"  Rows with gamma Barlow ratio > {threshold}: {gamma_over} ({gamma_over/total_rows*100:.1f}%)")
        
        return output_file
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return None


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate Barlow ratios for systematic error estimation.')
    parser.add_argument('--input', type=str, required=True,
                      help='Input CSV file (from syssrc_data_collector.py)')
    parser.add_argument('--output', type=str, default=None,
                      help='Output CSV file (default: barlow_<input_filename>)')
    parser.add_argument('--threshold', type=float, default=1.0,
                      help='Threshold for Barlow ratio significance (default: 1.0). Ratios above this are considered significant.')
    parser.add_argument('--ignore-diff', action='store_true',
                      help='Ignore diff_type and diff_bin columns in matching (useful for Lambda data)')
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose output')
    parser.add_argument('--debug', action='store_true',
                      help='Enable debug mode with detailed diagnostic information')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)
    
    # Show pandas version for debugging
    if args.debug:
        print(f"Using pandas version: {pd.__version__}")
        print(f"Using numpy version: {np.__version__}")
        
    # Process the file
    folder_name = Path(args.input).stem.replace('finalise_sys_', '')
    strip_diff = (folder_name == 'Lambda') or args.ignore_diff
    
    if strip_diff:
        print("Note: Ignoring diff_type and diff_bin columns in matching")
    
    process_file(args.input, args.output, args.threshold, args.verbose or args.debug)


if __name__ == "__main__":
    main()