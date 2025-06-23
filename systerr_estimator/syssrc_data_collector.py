#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import re
import sys
from pathlib import Path


def check_folder_name(folder_path):
    """Check if the folder name is valid (Lambda, Proton, or Hadron)"""
    folder_name = os.path.basename(folder_path)
    if folder_name not in ["Lambda", "Proton", "Hadron"]:
        print(f"Warning: Folder name '{folder_name}' is not one of the expected names (Lambda, Proton, Hadron)")
    return folder_name


def get_default_file(folder_path):
    """Get the path of the required default file in the folder based on folder name"""
    folder_name = os.path.basename(folder_path)
    
    if folder_name == "Proton":
        default_file = os.path.join(folder_path, "finalise_feeddown_dispose_default.csv")
        file_name = "finalise_feeddown_dispose_default.csv"
    else:  # Lambda or Hadron
        default_file = os.path.join(folder_path, "finalise_default.csv")
        file_name = "finalise_default.csv"
        
    if not os.path.exists(default_file):
        print(f"Error: Required file '{file_name}' not found in {folder_path}")
        sys.exit(1)
    return default_file


def get_task_files(folder_path, folder_name, exclude_default=False):
    """Get all task files according to the folder name pattern"""
    files = []
    
    if folder_name in ["Lambda", "Hadron"]:
        pattern = re.compile(r"finalise_(.+)\.csv$")
    elif folder_name == "Proton":
        pattern = re.compile(r"finalise_feeddown_dispose_(.+)\.csv$")
    else:
        return files
    
    for file in os.listdir(folder_path):
        match = pattern.match(file)
        if match:
            task_name = match.group(1)
            
            # Skip default files if requested
            if exclude_default and task_name == "default":
                continue
            
            full_path = os.path.join(folder_path, file)
            
            # Verify the file is not empty
            try:
                if os.path.getsize(full_path) == 0:
                    print(f"Warning: Skipping empty file {file}")
                    continue
                    
                # Quick check that it's a valid CSV
                with open(full_path, 'r') as f:
                    first_line = f.readline().strip()
                    if not first_line or ',' not in first_line:
                        print(f"Warning: File {file} does not appear to be a valid CSV, skipping")
                        continue
                        
                files.append((full_path, task_name))
                
            except Exception as e:
                print(f"Error checking file {file}: {e}")
                continue
    
    return files


def process_file(file_path, task_name, folder_name):
    """Process a single file and return a dataframe with standardized format"""
    try:
        df = pd.read_csv(file_path)
        
        # Verify required columns exist
        required_cols = ['centrality', 'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error in {file_path}: Missing required columns: {', '.join(missing_cols)}")
            return None
        
        # Check which format we have
        has_diff_cols = all(col in df.columns for col in ['diff_type', 'diff_bin'])
        
        # Create a copy of the original dataframe to avoid modification warnings
        result_df = df.copy()
        
        # Add source column (derived from task_name)
        result_df['source'] = task_name
        
        # If diff columns don't exist, add them as empty
        if not has_diff_cols:
            result_df['diff_type'] = ''
            result_df['diff_bin'] = ''
        
        # Ensure all required columns are present
        for col in ['source', 'centrality', 'diff_type', 'diff_bin', 
                    'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']:
            if col not in result_df.columns:
                print(f"Warning: Adding missing column '{col}' to result")
                result_df[col] = ''
        
        return result_df
        
    except pd.errors.EmptyDataError:
        print(f"Error: {file_path} is empty")
        return None
    except pd.errors.ParserError:
        print(f"Error: Cannot parse {file_path} as a CSV file")
        return None
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process CSV files for systematic error estimation.')
    parser.add_argument('--inputdir', type=str, default='../csv_data_point/Proton',
                        help='Input directory containing CSV files (default: ../csv_data_point/Proton)')
    parser.add_argument('--outputdir', type=str, default='./',
                        help='Output directory for processed files (default: ./)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.exists(args.inputdir):
        print(f"Error: Input directory {args.inputdir} does not exist")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outputdir, exist_ok=True)
    
    # Check folder name
    folder_name = check_folder_name(args.inputdir)
    
    # Get the default file path
    default_file_path = get_default_file(args.inputdir)
    
    # Get task files (excluding default since we'll handle it separately)
    task_files = get_task_files(args.inputdir, folder_name, exclude_default=True)
    
    if not task_files:
        print(f"Warning: No non-default task files found in {args.inputdir}")
        print("Only the default file will be processed")
    else:
        print(f"Found {len(task_files)} task files in {args.inputdir}")
    
    # Process all files
    all_results = []
    processed_count = 0
    failed_count = 0
    
    # First process the default file to ensure it's included
    print(f"Processing default file: {default_file_path}")
    default_result_df = process_file(default_file_path, "default", folder_name)
    if default_result_df is not None:
        all_results.append(default_result_df)
        processed_count += 1
        print("Successfully processed default file")
    else:
        failed_count += 1
        print("Error: Failed to process the default file, which is required")
        sys.exit(1)
    
    # Then process all other files
    for file_path, task_name in task_files:
        print(f"Processing {file_path} (task: {task_name})")
        result_df = process_file(file_path, task_name, folder_name)
        if result_df is not None:
            all_results.append(result_df)
            processed_count += 1
        else:
            failed_count += 1
            
    if failed_count > 0:
        print(f"Warning: {failed_count} files failed to process")
    
    print(f"Successfully processed {processed_count} files")
    
    # Combine all results
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Ensure columns are in the correct order
        column_order = ['source', 'centrality', 'diff_type', 'diff_bin', 
                         'pair_type', 'delta', 'delta_err', 'gamma', 'gamma_err']
        
        # Make sure all required columns exist
        for col in column_order:
            if col not in combined_df.columns:
                combined_df[col] = ''
                
        combined_df = combined_df[column_order]
        
        # Convert empty values to empty strings for consistency
        combined_df = combined_df.fillna('')
        
        # Create output filename
        output_file = os.path.join(args.outputdir, f"finalise_sys_{folder_name}.csv")
        
        # Save to CSV
        try:
            combined_df.to_csv(output_file, index=False)
            print(f"Successfully created {output_file} with {len(combined_df)} rows of data")
            
            if args.verbose:
                # Show summary information
                print("\nSummary of processed data:")
                print(f"Total sources: {combined_df['source'].nunique()}")
                print(f"Centrality bins: {combined_df['centrality'].nunique()}")
                print(f"Pair types: {combined_df['pair_type'].nunique()}")
                if combined_df['diff_type'].any():
                    print(f"Differential types: {combined_df['diff_type'].nunique()}")
        except Exception as e:
            print(f"Error saving output file: {e}")
            sys.exit(1)
    else:
        print("No data to process")


if __name__ == "__main__":
    main()