#!/usr/bin/env python3
"""
Convenient script to run feeddown comparison
Uses default paths and settings
"""

import os
import sys
from pathlib import Path

# Add current directory to path
sys.path.append(str(Path(__file__).parent))

from plot_feeddown_comparison import plot_all_diff_types_comparison, plot_feeddown_comparison

def main():
    """Main function: run feeddown comparison plotting"""
    
    # Data path configuration
    base_dir = Path(__file__).parent
    
    # Original data path (from dataset_merger)
    original_data_path = base_dir / "../dataset_merger/Merged/Proton/finalise_default.csv"
    
    # Feeddown corrected data path
    feeddown_data_path = base_dir / "Proton/finalise_feeddown_dispose_default.csv"
    
    # Output directory
    output_dir = base_dir / "plots/feeddown_comparison"
    
    # Check if files exist
    if not original_data_path.exists():
        print(f"Error: Original data file not found: {original_data_path}")
        return
        
    if not feeddown_data_path.exists():
        print(f"Error: Feeddown data file not found: {feeddown_data_path}")
        return
    
    print("Starting feeddown comparison plot generation...")
    print(f"Original data: {original_data_path}")
    print(f"Feeddown data: {feeddown_data_path}")
    print(f"Output directory: {output_dir}")
    
    # Generate comparison plots for all differential types
    print("\nGenerating comparison plots for all differential types...")
    plot_all_diff_types_comparison(
        str(original_data_path), str(feeddown_data_path), 
        str(output_dir)
    )
    
    print("\nAll feeddown comparison plots generated successfully!")
    print(f"Plots saved in: {output_dir}")

if __name__ == "__main__":
    main()