#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data collection script for differential plotting data (Proton only)
Collects all central values, statistical errors, and systematic errors
into a single CSV file for plotting

Usage:
  ./collect_diff.py [--smooth <type>]

Options:
  --smooth <type>  Use smoothed systematic errors: 'exp' for exponential fit,
                  'pol2' for polynomial fit, or omit for default values.
"""

import pandas as pd
from pathlib import Path
import argparse


def collect_data(smooth_type=None):
    """
    Collect differential data for Proton pairs
    
    Args:
        smooth_type (str, optional): Type of smoothing to apply to systematic errors.
                                    Options: 'exp', 'pol2', or None for default.
    """
    base_dir = Path("../csv_data_point")
    
    # Load Proton data (only source needed since we're only processing Proton data)
    proton_default = pd.read_csv(base_dir / "Proton" / "finalise_feeddown_dispose_default.csv")
    proton_combined = pd.read_csv(base_dir / "Proton.csv")  # Contains both central and syst errors
    
    # Filter out Intg data, keeping only differential data
    proton_default = proton_default[proton_default["diff_type"] != "Intg"].copy()
    proton_combined = proton_combined[proton_combined["diff_type"] != "Intg"].copy()
    
    # Filter centrality <= 60%
    proton_default = proton_default[proton_default["centrality"] <= 60].copy()
    proton_combined = proton_combined[proton_combined["centrality"] <= 60].copy()
    
    # Prepare output data
    output_data = []
    
    # Process Proton data (use combined file for both central and systematic errors)
    print("Processing Proton differential data...")
    for _, row in proton_combined.iterrows():
        # Select appropriate systematic error columns based on smooth_type
        delta_syst_col = "delta_syst_err"
        gamma_syst_col = "gamma_syst_err"
        
        if smooth_type == 'exp' and "delta_syst_err_expfit" in row:
            delta_syst_col = "delta_syst_err_expfit"
            gamma_syst_col = "gamma_syst_err_expfit"
        elif smooth_type == 'pol2' and "delta_syst_err_pol2" in row:
            delta_syst_col = "delta_syst_err_pol2"
            gamma_syst_col = "gamma_syst_err_pol2"
        
        output_data.append({
            "diff_type": row["diff_type"],
            "diff_bin": row["diff_bin"],
            "centrality": row["centrality"],
            "pair_type": row["pair_type"],
            "delta": row["delta"],
            "delta_err": row["delta_err"],
            "delta_syst_err": row[delta_syst_col] if delta_syst_col in row else row["delta_syst_err"],
            "gamma": row["gamma"],
            "gamma_err": row["gamma_err"],
            "gamma_syst_err": row[gamma_syst_col] if gamma_syst_col in row else row["gamma_syst_err"]
        })
    
    return pd.DataFrame(output_data)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Collect differential correlator data for plotting (Proton only)")
    parser.add_argument("--smooth", choices=["exp", "pol2"], help="Use smoothed systematic errors: 'exp' for exponential fit, 'pol2' for polynomial fit")
    args = parser.parse_args()
    
    smooth_type = args.smooth
    if smooth_type:
        print(f"Using {smooth_type} smoothing for systematic errors")
    else:
        print("Using default systematic errors (no smoothing)")
    
    print("Collecting differential plotting data for Proton...")
    
    # Collect all data
    df = collect_data(smooth_type)
    
    # Sort by diff_type, diff_bin, centrality and pair_type for easier reading
    df = df.sort_values(["diff_type", "diff_bin", "centrality", "pair_type"]).reset_index(drop=True)
    
    # Save to CSV
    if smooth_type:
        output_path = Path(f"../csv_data_point/diff_{smooth_type}.csv")
    else:
        output_path = Path("../csv_data_point/diff.csv")
    
    df.to_csv(output_path, index=False)
    
    print(f"Data collected and saved to: {output_path}")
    print(f"Total rows: {len(df)}")
    
    # Print summary
    print("\nData summary:")
    for diff_type in df["diff_type"].unique():
        diff_type_data = df[df["diff_type"] == diff_type]
        print(f"  {diff_type}: {len(diff_type_data)} rows")
        for diff_bin in sorted(diff_type_data["diff_bin"].unique()):
            bin_data = diff_type_data[diff_type_data["diff_bin"] == diff_bin]
            print(f"    Bin {diff_bin}: {len(bin_data)} rows")
            for pair_type in bin_data["pair_type"].unique():
                count = len(bin_data[bin_data["pair_type"] == pair_type])
                print(f"      {pair_type}: {count} rows")
    
    # Print sample data
    print("\nSample data (first 5 rows):")
    print(df.head().to_string(index=False))
    
    # Check for any missing systematic errors
    missing_delta_syst = df["delta_syst_err"].isna().sum()
    missing_gamma_syst = df["gamma_syst_err"].isna().sum()
    
    if missing_delta_syst > 0 or missing_gamma_syst > 0:
        print(f"\nWarning: Missing systematic errors:")
        print(f"  Missing delta_syst_err: {missing_delta_syst}")
        print(f"  Missing gamma_syst_err: {missing_gamma_syst}")
    else:
        print("\nAll systematic errors found successfully!")


if __name__ == "__main__":
    main()