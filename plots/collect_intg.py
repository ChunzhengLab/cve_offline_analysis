#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data collection script for Intg plotting data
Collects all central values, statistical errors, and systematic errors
into a single CSV file for plotting

Usage:
  ./collect_intg.py [--smooth <type>]

Options:
  --smooth <type>  Use smoothed systematic errors: 'exp' for exponential fit,
                  'pol2' for polynomial fit, or omit for default values.
"""

import pandas as pd
from pathlib import Path
import argparse



def collect_data(smooth_type=None):
    """
    Collect all data needed for plotting
    
    Args:
        smooth_type (str, optional): Type of smoothing to apply to systematic errors.
                                    Options: 'exp', 'pol2', or None for default.
    """
    base_dir = Path("../csv_data_point")
    
    # Load central values and statistical errors
    lambda_default = pd.read_csv(base_dir / "Lambda" / "finalise_default.csv")
    hadron_default = pd.read_csv(base_dir / "Hadron" / "finalise_default.csv")
    proton_default = pd.read_csv(base_dir / "Proton" / "finalise_feeddown_dispose_default.csv")
    
    # Load systematic errors
    lambda_syst = pd.read_csv(base_dir / "Lambda.csv")
    hadron_syst = pd.read_csv(base_dir / "Hadron.csv")
    proton_combined = pd.read_csv(base_dir / "Proton.csv")  # Contains both central and syst errors
    
    # Filter for Intg data only
    hadron_default = hadron_default[hadron_default["diff_type"] == "Intg"].copy()
    proton_default = proton_default[proton_default["diff_type"] == "Intg"].copy()
    hadron_syst = hadron_syst[hadron_syst["diff_type"] == "Intg"].copy()
    proton_combined = proton_combined[proton_combined["diff_type"] == "Intg"].copy()
    
    # Filter centrality <= 60%
    lambda_default = lambda_default[lambda_default["centrality"] <= 60].copy()
    hadron_default = hadron_default[hadron_default["centrality"] <= 60].copy()
    proton_default = proton_default[proton_default["centrality"] <= 60].copy()
    lambda_syst = lambda_syst[lambda_syst["centrality"] <= 60].copy()
    hadron_syst = hadron_syst[hadron_syst["centrality"] <= 60].copy()
    proton_combined = proton_combined[proton_combined["centrality"] <= 60].copy()
    
    # Prepare output data
    output_data = []
    
    # Process Lambda data
    print("Processing Lambda data...")
    for _, row in lambda_default.iterrows():
        centrality = row["centrality"]
        pair_type = row["pair_type"]
        
        # Find corresponding systematic error
        syst_row = lambda_syst[
            (lambda_syst["centrality"] == centrality) & 
            (lambda_syst["pair_type"] == pair_type)
        ]
        
        if len(syst_row) > 0:
            syst_row = syst_row.iloc[0]
            
            # Select appropriate systematic error columns based on smooth_type
            delta_syst_col = "delta_syst_err"
            gamma_syst_col = "gamma_syst_err"
            
            if smooth_type == 'exp' and f"delta_syst_err_expfit" in syst_row:
                delta_syst_col = "delta_syst_err_expfit"
                gamma_syst_col = "gamma_syst_err_expfit"
            elif smooth_type == 'pol2' and f"delta_syst_err_pol2" in syst_row:
                delta_syst_col = "delta_syst_err_pol2"
                gamma_syst_col = "gamma_syst_err_pol2"
            
            output_data.append({
                "particle": "Lambda",
                "centrality": centrality,
                "pair_type": pair_type,
                "delta": row["delta"],
                "delta_err": row["delta_err"],
                "delta_syst_err": syst_row[delta_syst_col] if delta_syst_col in syst_row else syst_row["delta_syst_err"],
                "gamma": row["gamma"],
                "gamma_err": row["gamma_err"],
                "gamma_syst_err": syst_row[gamma_syst_col] if gamma_syst_col in syst_row else syst_row["gamma_syst_err"]
            })
    
    # Process Proton data (use combined file for both central and systematic errors)
    print("Processing Proton data...")
    for _, row in proton_combined.iterrows():
        # Select appropriate systematic error columns based on smooth_type
        delta_syst_col = "delta_syst_err"
        gamma_syst_col = "gamma_syst_err"
        
        if smooth_type == 'exp' and f"delta_syst_err_expfit" in row:
            delta_syst_col = "delta_syst_err_expfit"
            gamma_syst_col = "gamma_syst_err_expfit"
        elif smooth_type == 'pol2' and f"delta_syst_err_pol2" in row:
            delta_syst_col = "delta_syst_err_pol2"
            gamma_syst_col = "gamma_syst_err_pol2"
        
        output_data.append({
            "particle": "Proton",
            "centrality": row["centrality"],
            "pair_type": row["pair_type"],
            "delta": row["delta"],
            "delta_err": row["delta_err"],
            "delta_syst_err": row[delta_syst_col] if delta_syst_col in row else row["delta_syst_err"],
            "gamma": row["gamma"],
            "gamma_err": row["gamma_err"],
            "gamma_syst_err": row[gamma_syst_col] if gamma_syst_col in row else row["gamma_syst_err"]
        })
    
    # Process Hadron data
    print("Processing Hadron data...")
    for _, row in hadron_default.iterrows():
        centrality = row["centrality"]
        pair_type = row["pair_type"]
        
        # Find corresponding systematic error
        syst_row = hadron_syst[
            (hadron_syst["centrality"] == centrality) & 
            (hadron_syst["pair_type"] == pair_type)
        ]
        
        if len(syst_row) > 0:
            syst_row = syst_row.iloc[0]
            
            # Select appropriate systematic error columns based on smooth_type
            delta_syst_col = "delta_syst_err"
            gamma_syst_col = "gamma_syst_err"
            
            if smooth_type == 'exp' and f"delta_syst_err_expfit" in syst_row:
                delta_syst_col = "delta_syst_err_expfit"
                gamma_syst_col = "gamma_syst_err_expfit"
            elif smooth_type == 'pol2' and f"delta_syst_err_pol2" in syst_row:
                delta_syst_col = "delta_syst_err_pol2"
                gamma_syst_col = "gamma_syst_err_pol2"
            
            output_data.append({
                "particle": "Hadron",
                "centrality": centrality,
                "pair_type": pair_type,
                "delta": row["delta"],
                "delta_err": row["delta_err"],
                "delta_syst_err": syst_row[delta_syst_col] if delta_syst_col in syst_row else syst_row["delta_syst_err"],
                "gamma": row["gamma"],
                "gamma_err": row["gamma_err"],
                "gamma_syst_err": syst_row[gamma_syst_col] if gamma_syst_col in syst_row else syst_row["gamma_syst_err"]
            })
    
    return pd.DataFrame(output_data)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Collect integrated correlator data for plotting")
    parser.add_argument("--smooth", choices=["exp", "pol2"], help="Use smoothed systematic errors: 'exp' for exponential fit, 'pol2' for polynomial fit")
    args = parser.parse_args()
    
    smooth_type = args.smooth
    if smooth_type:
        print(f"Using {smooth_type} smoothing for systematic errors")
    else:
        print("Using default systematic errors (no smoothing)")
    
    print("Collecting Intg plotting data...")
    
    # Collect all data
    df = collect_data(smooth_type)
    
    # Sort by particle, pair_type, and centrality for easier reading
    df = df.sort_values(["particle", "pair_type", "centrality"]).reset_index(drop=True)
    
    # Save to CSV
    if smooth_type:
        output_path = Path(f"../csv_data_point/intg_{smooth_type}.csv")
    else:
        output_path = Path("../csv_data_point/intg.csv")
    
    df.to_csv(output_path, index=False)
    
    print(f"Data collected and saved to: {output_path}")
    print(f"Total rows: {len(df)}")
    
    # Print summary
    print("\nData summary:")
    for particle in df["particle"].unique():
        particle_data = df[df["particle"] == particle]
        print(f"  {particle}: {len(particle_data)} rows")
        for pair_type in particle_data["pair_type"].unique():
            count = len(particle_data[particle_data["pair_type"] == pair_type])
            print(f"    {pair_type}: {count} rows")
    
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