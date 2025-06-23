#!/usr/bin/env python3
"""
Dataset merger for combining LHC18q and LHC18r results.
Merges CSV files using weighted average method.

Author: Dataset Analysis Team
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path
import sys
import argparse

def weighted_merge(value1, err1, value2, err2):
    """
    Merge two measurements using weighted average.

    Returns:
        merged_value, merged_error
    """
    if pd.isna(value1) or pd.isna(err1) or pd.isna(value2) or pd.isna(err2):
        return np.nan, np.nan

    if err1 <= 0 or err2 <= 0:
        return np.nan, np.nan

    weight1 = 1 / (err1 ** 2)
    weight2 = 1 / (err2 ** 2)

    merged_value = (value1 * weight1 + value2 * weight2) / (weight1 + weight2)
    merged_error = np.sqrt(1 / (weight1 + weight2))

    return merged_value, merged_error

def get_task_names(particle_type):
    """
    Get all task names from finalise_*.csv files for a given particle type.
    """
    # Use different directory for Lambda particles
    base_dir = "../obv_fit_lambda2" if particle_type == "Lambda" else "../obv_fit"
    period = "LHC18q"
    pattern = f"{base_dir}/{period}/{particle_type}/finalise_*.csv"
    files = glob.glob(pattern)
    task_names = []

    if not files:
        print(f"    Warning: No finalise_*.csv files found in {base_dir}/{period}/{particle_type}/")
        return task_names

    for file_path in files:
        filename = os.path.basename(file_path)
        # Extract task name from finalise_[task].csv
        task_name = filename.replace('finalise_', '').replace('.csv', '')
        task_names.append(task_name)

    return task_names

def load_and_standardize_data(file_path, particle_type):
    """
    Load CSV data and standardize column names based on particle type.
    """
    try:
        df = pd.read_csv(file_path)
        print(f"      Loaded: {file_path} ({len(df)} rows)")
    except FileNotFoundError:
        print(f"      Warning: File not found: {file_path}")
        return None
    except Exception as e:
        print(f"      Error loading {file_path}: {e}")
        return None

    if df.empty:
        print(f"      Warning: Empty file: {file_path}")
        return None

    # Remove resolution column if present (especially for Lambda data)
    if "resolution" in df.columns:
        df = df.drop(columns=["resolution"])
        print(f"      Removed resolution column")
        
    # Verify required columns are present
    required_cols = ["centrality", "pair_type", "delta", "delta_err", "gamma", "gamma_err"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"      Error: Missing required columns: {missing}")
        return None
        
    # Ensure centrality is always numeric for consistent handling
    df["centrality"] = pd.to_numeric(df["centrality"])
    
    return df

def merge_datasets_for_particle(particle_type):
    """
    Merge datasets for a specific particle type.
    """
    task_names = get_task_names(particle_type)

    if not task_names:
        print(f"  No task files found for {particle_type}")
        return

    print(f"  Found {len(task_names)} tasks: {task_names}")

    # Create output directory
    output_dir = f"./Merged/{particle_type}"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"  Output directory: {output_dir}")

    successful_merges = 0

    for task in task_names:
        print(f"    Processing task: {task}")

        # File paths - use different directory for Lambda particles
        base_dir = "../obv_fit_lambda2" if particle_type == "Lambda" else "../obv_fit"
        file1_path = f"{base_dir}/LHC18q/{particle_type}/finalise_{task}.csv"
        file2_path = f"{base_dir}/LHC18r/{particle_type}/finalise_{task}.csv"
        
        print(f"      Looking for files in {base_dir} directory")

        # Load data
        df1 = load_and_standardize_data(file1_path, particle_type)
        df2 = load_and_standardize_data(file2_path, particle_type)

        if df1 is None or df2 is None:
            print(f"      Skipping {task} due to missing or invalid files")
            continue

        # Merge the datasets
        merged_df = merge_two_dataframes(df1, df2, particle_type)

        if merged_df is not None:
            # Save merged result
            output_path = f"{output_dir}/finalise_{task}.csv"
            try:
                merged_df.to_csv(output_path, index=False)
                print(f"      ✓ Saved: {output_path} ({len(merged_df)} rows)")
                successful_merges += 1
            except Exception as e:
                print(f"      Error saving {output_path}: {e}")
        else:
            print(f"      Failed to merge {task}")

    print(f"  Completed: {successful_merges}/{len(task_names)} tasks merged successfully")

def merge_two_dataframes(df1, df2, particle_type):
    """
    Merge two dataframes with the same structure.
    """
    # Determine grouping columns based on particle type
    if particle_type == "Lambda":
        group_cols = ["centrality", "pair_type"]
    else:
        group_cols = ["centrality", "diff_type", "diff_bin", "pair_type"]

    # Check if required columns exist
    missing_cols_df1 = [col for col in group_cols if col not in df1.columns]
    missing_cols_df2 = [col for col in group_cols if col not in df2.columns]

    if missing_cols_df1 or missing_cols_df2:
        print(f"        Error: Missing columns. DF1: {missing_cols_df1}, DF2: {missing_cols_df2}")
        return None

    # Filter out Del rows for merging (will be calculated later)
    df1_no_del = df1[df1["pair_type"] != "Del"].copy()
    df2_no_del = df2[df2["pair_type"] != "Del"].copy()

    print(f"        Merging {len(df1_no_del)} rows from LHC18q with {len(df2_no_del)} rows from LHC18r")

    # Merge on the grouping columns
    merged = pd.merge(df1_no_del, df2_no_del, on=group_cols, suffixes=('_1', '_2'))

    if merged.empty:
        print("        Warning: No matching rows found for merging")
        return None
        
    # Ensure centrality is always a numeric type
    merged["centrality"] = pd.to_numeric(merged["centrality"])

    # Apply weighted merging for delta and gamma
    print(f"        Performing weighted merge of measurements")
    merged['delta'], merged['delta_err'] = zip(*merged.apply(
        lambda row: weighted_merge(row['delta_1'], row['delta_err_1'],
                                 row['delta_2'], row['delta_err_2']), axis=1))

    merged['gamma'], merged['gamma_err'] = zip(*merged.apply(
        lambda row: weighted_merge(row['gamma_1'], row['gamma_err_1'],
                                 row['gamma_2'], row['gamma_err_2']), axis=1))

    # Keep only the final columns
    final_cols = group_cols + ['delta', 'delta_err', 'gamma', 'gamma_err']
    merged = merged[final_cols]

    # Calculate Del (OS - SS) for each group
    del_rows = []
    print(f"        Calculating Del (OS-SS) rows")

    if particle_type == "Lambda":
        # Simpler approach for Lambda data
        for centrality, group_df in merged.groupby("centrality"):
            os_row = group_df[group_df["pair_type"] == "OS"]
            ss_row = group_df[group_df["pair_type"] == "SS"]

            if len(os_row) == 1 and len(ss_row) == 1:
                os_row = os_row.iloc[0]
                ss_row = ss_row.iloc[0]

                # Calculate Del = OS - SS
                del_delta = os_row["delta"] - ss_row["delta"]
                del_delta_err = np.sqrt(os_row["delta_err"]**2 + ss_row["delta_err"]**2)

                del_gamma = os_row["gamma"] - ss_row["gamma"]
                del_gamma_err = np.sqrt(os_row["gamma_err"]**2 + ss_row["gamma_err"]**2)

                # Create Del row
                del_row = {
                    "centrality": centrality,  # Use the scalar centrality value directly
                    "pair_type": "Del",
                    "delta": del_delta,
                    "delta_err": del_delta_err,
                    "gamma": del_gamma,
                    "gamma_err": del_gamma_err
                }
                del_rows.append(del_row)
    else:
        # Original approach for other particle types
        group_by_cols = ["centrality", "diff_type", "diff_bin"]
        for group_vals, group_df in merged.groupby(group_by_cols):
            os_row = group_df[group_df["pair_type"] == "OS"]
            ss_row = group_df[group_df["pair_type"] == "SS"]

            if len(os_row) == 1 and len(ss_row) == 1:
                os_row = os_row.iloc[0]
                ss_row = ss_row.iloc[0]

                # Calculate Del = OS - SS
                del_delta = os_row["delta"] - ss_row["delta"]
                del_delta_err = np.sqrt(os_row["delta_err"]**2 + ss_row["delta_err"]**2)

                del_gamma = os_row["gamma"] - ss_row["gamma"]
                del_gamma_err = np.sqrt(os_row["gamma_err"]**2 + ss_row["gamma_err"]**2)

                # Create Del row
                del_row = {
                    "centrality": group_vals[0],
                    "diff_type": group_vals[1],
                    "diff_bin": group_vals[2],
                    "pair_type": "Del",
                    "delta": del_delta,
                    "delta_err": del_delta_err,
                    "gamma": del_gamma,
                    "gamma_err": del_gamma_err
                }
                del_rows.append(del_row)

    # Add Del rows to the merged dataframe
    if del_rows:
        del_df = pd.DataFrame(del_rows)
        # Ensure centrality is numeric in del_df too
        del_df["centrality"] = pd.to_numeric(del_df["centrality"])
        merged = pd.concat([merged, del_df], ignore_index=True)

    # Sort the dataframe for consistent output
    if particle_type == "Lambda":
        # Ensure centrality is numeric and define ordering for pair_type
        merged["centrality"] = pd.to_numeric(merged["centrality"])
        # Create a custom ordering for pair_type
        pair_type_order = {"SS": 0, "OS": 1, "Del": 2}
        merged["pair_type_order"] = merged["pair_type"].map(pair_type_order)
        # Sort by centrality and then by pair_type_order
        merged = merged.sort_values(["centrality", "pair_type_order"])
        # Drop the temporary ordering column
        merged = merged.drop(columns=["pair_type_order"]).reset_index(drop=True)
    else:
        merged = merged.sort_values(["centrality", "diff_type", "diff_bin", "pair_type"]).reset_index(drop=True)

    return merged

def main():
    """
    Main function to merge datasets for specified particle types.
    """
    parser = argparse.ArgumentParser(
        description="Merge LHC18q and LHC18r datasets using weighted average method",
        epilog="Example: python dataset_merger.py -p Proton"
    )
    parser.add_argument(
        "--particle", "-p",
        type=str,
        choices=["Proton", "Pion", "Hadron", "Lambda"],
        help="Specify particle type to merge. If not provided, all particles will be processed."
    )


    args = parser.parse_args()

    if args.particle:
        particle_types = [args.particle]
        print(f"Starting dataset merger for {args.particle}...")
    else:
        particle_types = ["Lambda", "Proton", "Pion", "Hadron"]
        print("Starting dataset merger for all particle types...")

    print("=" * 50)

    for particle_type in particle_types:
        print(f"\nProcessing {particle_type}...")
        print(f"Using {'../obv_fit_lambda2' if particle_type == 'Lambda' else '../obv_fit'} as base directory")
        merge_datasets_for_particle(particle_type)

    print("\n" + "=" * 50)
    print("Dataset merging completed!")

    # Check if any merged files were created
    merged_dir = Path("./Merged")
    if merged_dir.exists():
        total_files = sum(len(list(p.glob("*.csv"))) for p in merged_dir.iterdir() if p.is_dir())
        print(f"Total merged files created: {total_files}")
    else:
        print("No merged files were created.")

if __name__ == "__main__":
    main()
