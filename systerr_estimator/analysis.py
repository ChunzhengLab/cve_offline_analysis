#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path

def run_command(cmd, description=None, exit_on_error=True):
    """Run a shell command and print its output"""
    if description:
        print(f"\n=== {description} ===")

    print(f"Running: {cmd}")
    start_time = time.time()

    try:
        result = subprocess.run(cmd, shell=True, check=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(result.stdout)

        elapsed = time.time() - start_time
        print(f"Command completed successfully in {elapsed:.2f} seconds")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print(f"Output: {e.output}")

        if exit_on_error:
            print("Exiting due to error")
            sys.exit(1)
        return False

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Run the complete systematic error analysis pipeline')

    # Required positional argument - particle type
    parser.add_argument('particle_type', type=str, choices=['Lambda', 'Proton', 'Hadron'],
                      help='Particle type to analyze (Lambda, Proton, or Hadron)')

    # Input directory
    parser.add_argument('--inputdir', type=str, default='../csv_data_point',
                      help='Input directory containing CSV files (default: ../csv_data_point)')

    # Output directory
    parser.add_argument('--outputdir', type=str, default='./',
                      help='Output directory for processed files (default: ./)')

    # Parameters for different steps
    parser.add_argument('--threshold', type=float, default=1.0,
                      help='Threshold for Barlow ratio significance (default: 1.0)')
    parser.add_argument('--min-occurrences', type=int, default=2,
                      help='Minimum number of occurrences to consider a source as systematic (default: 2)')
    parser.add_argument('--max-cent', type=int, default=60,
                      help='Maximum centrality value to consider (default: 60)')

    # Output formatting options
    parser.add_argument('--formatted-contrib', type=str, default=None,
                      help='Custom filename for formatted contribution fractions (default: contrib_formatted_<particle_type>.txt)')

    # Control flags
    parser.add_argument('--skip-collector', action='store_true',
                      help='Skip the data collection step')
    parser.add_argument('--skip-barlow', action='store_true',
                      help='Skip the Barlow ratio calculation step')
    parser.add_argument('--skip-syssrc', action='store_true',
                      help='Skip the systematic source determination step')
    parser.add_argument('--skip-contrib', action='store_true',
                      help='Skip the contribution fraction calculation step')
    parser.add_argument('--ignore-diff', action='store_true',
                      help='Ignore diff_type and diff_bin columns in matching (useful for Lambda data)')
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose output in all steps')

    return parser.parse_args()

def main():
    """Main function to run the complete analysis pipeline"""
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.outputdir, exist_ok=True)

    # Input directory for the specific particle type
    particle_input_dir = os.path.join(args.inputdir, args.particle_type)

    # Define file paths
    base_output_dir = args.outputdir

    # Step 1: Collect data from CSV files
    collector_output = os.path.join(base_output_dir, f"finalise_sys_{args.particle_type}.csv")

    if not args.skip_collector:
        collector_cmd = (f"python3 syssrc_data_collector.py "
                        f"--inputdir {particle_input_dir} "
                        f"--outputdir {base_output_dir}")

        run_command(collector_cmd, "Collecting data from CSV files")
    else:
        print("\n=== Skipping data collection step ===")
        # Check if the file exists
        if not os.path.exists(collector_output):
            print(f"Error: {collector_output} does not exist. Cannot skip data collection step.")
            sys.exit(1)

    # Step 2: Calculate Barlow ratios
    barlow_output = os.path.join(base_output_dir, f"barlow_finalise_sys_{args.particle_type}.csv")

    if not args.skip_barlow:
        barlow_cmd = (f"python3 barlow_ratio_calculator.py "
                    f"--input {collector_output} "
                    f"--output {barlow_output} "
                    f"--threshold {args.threshold}")

        if args.ignore_diff:
            barlow_cmd += " --ignore-diff"

        if args.verbose:
            barlow_cmd += " --verbose"

        run_command(barlow_cmd, "Calculating Barlow ratios")
    else:
        print("\n=== Skipping Barlow ratio calculation step ===")
        # Check if the file exists
        if not os.path.exists(barlow_output):
            print(f"Error: {barlow_output} does not exist. Cannot skip Barlow ratio calculation step.")
            sys.exit(1)

    # Step 3: Determine systematic sources
    syssrc_output = os.path.join(base_output_dir, f"syssrc_barlow_finalise_sys_{args.particle_type}.csv")

    if not args.skip_syssrc:
        syssrc_cmd = (f"python3 sys_src_determiner.py "
                    f"--input {barlow_output} "
                    f"--output {syssrc_output} "
                    f"--min-occurrences {args.min_occurrences} "
                    f"--max-cent {args.max_cent}")

        if args.verbose:
            syssrc_cmd += " --verbose"

        run_command(syssrc_cmd, "Determining systematic sources")
    else:
        print("\n=== Skipping systematic source determination step ===")
        # Check if the file exists
        if not os.path.exists(syssrc_output):
            print(f"Error: {syssrc_output} does not exist. Cannot skip systematic source determination step.")
            sys.exit(1)

    # Step 4: Calculate contribution fractions
    contrib_output = os.path.join(base_output_dir, f"contrib_syssrc_barlow_finalise_sys_{args.particle_type}.csv")

    # Set formatted contribution output filename
    if args.formatted_contrib:
        formatted_contrib_output = os.path.join(base_output_dir, args.formatted_contrib)
    else:
        formatted_contrib_output = os.path.join(base_output_dir, f"contrib_formatted_{args.particle_type}.txt")

    if not args.skip_contrib:
        contrib_cmd = (f"python3 contribution_fraction.py "
                    f"--input {syssrc_output} "
                    f"--output {contrib_output} "
                    f"--formatted-output {formatted_contrib_output}")

        if args.verbose:
            contrib_cmd += " --verbose"

        run_command(contrib_cmd, "Calculating contribution fractions")
    else:
        print("\n=== Skipping contribution fraction calculation step ===")

    # Step 5: Generate final output
    final_output = os.path.join(base_output_dir, f"{args.particle_type}.csv")

    finalise_cmd = (f"python3 finalise.py "
                  f"--input {syssrc_output} "
                  f"--output {final_output} "
                  f"--original {collector_output}")

    if args.verbose:
        finalise_cmd += " --verbose"

    run_command(finalise_cmd, "Generating final systematic error values")

    print("\n=== Analysis Complete ===")
    print(f"Final output file: {final_output}")

    # Print a summary of all generated files
    print("\nGenerated files:")
    print(f"1. Data collection: {collector_output}")
    print(f"2. Barlow ratios: {barlow_output}")
    print(f"3. Systematic sources: {syssrc_output}")
    print(f"4. Contribution fractions: {contrib_output}")
    print(f"   Formatted contributions: {formatted_contrib_output}")
    print(f"5. Final results: {final_output}")

if __name__ == "__main__":
    main()
