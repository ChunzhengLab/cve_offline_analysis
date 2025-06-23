#!/usr/bin/env python3
"""
Feeddown dispose script for Proton data correction.
Applies feeddown correction using Lambda secondary data.

Version: 1.0.0
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path
import argparse
import sys


def calculate_fp(centrality):
    """
    Calculate the primary fraction fp from CSV data for given centrality.
    First average 18q and 18r data, then average p and ap.
    """
    # 读取CSV文件
    df = pd.read_csv("../src_fraction_fit/results/integration_results.csv")

    # 筛选指定中心度的数据
    df_cent = df[df['Centrality'] == centrality]

    # 计算18q和18r的平均值
    df_avg = df_cent.groupby('Particle').agg({
        'Primary/Total': 'mean',
        'Primary/Total Error': lambda x: np.sqrt(np.sum(x**2)) / len(x)  # 误差传播
    }).reset_index()

    # 获取质子和反质子的数据
    p_data = df_avg[df_avg['Particle'] == 'proton'].iloc[0]
    ap_data = df_avg[df_avg['Particle'] == 'antiproton'].iloc[0]

    # 计算最终的平均值
    fp = (p_data['Primary/Total'] + ap_data['Primary/Total']) / 2
    fp_err = np.sqrt(p_data['Primary/Total Error']**2 + ap_data['Primary/Total Error']**2) / 2

    return fp, fp_err

def calculate_ratio_lpsec_ll():
    """
    Return the ratio between O_lpsec and O_ll, and its error.
    Default: ratio = 0.9033 ± 0.0436
    """
    ratio = 1.
    ratio_err = 0.
    return ratio, ratio_err

def feeddown_correction(O_data, sigma_D, O_ll, sigma_ll, fp, sigma_f):
    """
    Apply feeddown correction using:
    O_lpsec = ratio * O_ll
    O_lp = (O_data - (1 - fp) * O_lpsec) / fp

    Propagate error including ratio, O_ll, fp, and O_data.
    """
    if pd.isna(O_data) or pd.isna(O_ll) or fp == 0:
        return np.nan, np.nan

    # Get ratio and its error
    ratio, ratio_err = calculate_ratio_lpsec_ll()

    # Compute O_lpsec
    O_lpsec = ratio * O_ll

    # Main value of O_lp
    O_lp = (O_data - (1 - fp) * O_lpsec) / fp

    # Error propagation
    term1 = sigma_D**2
    term2 = ((1 - fp) * ratio * sigma_ll)**2
    term3 = ((1 - fp) * O_ll * ratio_err)**2
    term4 = ((O_lpsec - O_lp)**2) * sigma_f**2

    sigma_Olp = (1 / fp) * np.sqrt(term1 + term2 + term3 + term4)

    return O_lp, sigma_Olp

def load_lambda_data(task_name, lambda_dir="../dataset_merger/Merged/Lambda"):
    """
    Load Lambda data for the given task. If not found, use default.
    """
    lambda_file = f"{lambda_dir}/finalise_{task_name}.csv"

    if os.path.exists(lambda_file):
        print(f"    Loading Lambda data: {lambda_file}")
        df = pd.read_csv(lambda_file)
    else:
        default_file = f"{lambda_dir}/finalise_default.csv"
        print(f"    Warning: {lambda_file} not found, using {default_file}")
        if os.path.exists(default_file):
            df = pd.read_csv(default_file)
        else:
            print(f"    Error: Default Lambda file not found: {default_file}")
            return None

    # Standardize column names if needed
    if "gamma_err" in df.columns:
        df = df.rename(columns={"gamma_err": "gamma_err"})

    # Remove resolution column if present
    if "resolution" in df.columns:
        df = df.drop(columns=["resolution"])

    return df

def get_lambda_values(lambda_df, centrality, pair_type):
    """
    Get Lambda values for given centrality and pair_type.
    """
    if lambda_df is None:
        return np.nan, np.nan, np.nan, np.nan

    mask = (lambda_df["centrality"] == centrality) & (lambda_df["pair_type"] == pair_type)
    matching_rows = lambda_df[mask]

    if len(matching_rows) == 0:
        return np.nan, np.nan, np.nan, np.nan

    row = matching_rows.iloc[0]
    return row["delta"], row["delta_err"], row["gamma"], row["gamma_err"]

def process_proton_data(proton_df, lambda_df, fp, fp_err):
    """
    Process Proton data with feeddown correction.
    """
    if proton_df is None or lambda_df is None:
        return None

    # Create a copy for processing
    result_df = proton_df.copy()

    # Process each row (except Del rows)
    for idx, row in result_df.iterrows():
        if row["pair_type"] == "Del":
            continue  # Skip Del rows, will be calculated later

        centrality = row["centrality"]
        pair_type = row["pair_type"]

        # Get Lambda secondary values
        lambda_delta, lambda_delta_err, lambda_gamma, lambda_gamma_err = get_lambda_values(
            lambda_df, centrality, pair_type
        )

        if pd.isna(lambda_delta) or pd.isna(lambda_gamma):
            print(f"      Warning: No Lambda data for centrality={centrality}, pair_type={pair_type}")
            continue

        # Apply feeddown correction for delta
        corrected_delta, corrected_delta_err = feeddown_correction(
            row["delta"], row["delta_err"],
            lambda_delta, lambda_delta_err,
            fp, fp_err
        )

        # Apply feeddown correction for gamma
        corrected_gamma, corrected_gamma_err = feeddown_correction(
            row["gamma"], row["gamma_err"],
            lambda_gamma, lambda_gamma_err,
            fp, fp_err
        )

        # Update the row
        result_df.at[idx, "delta"] = corrected_delta
        result_df.at[idx, "delta_err"] = corrected_delta_err
        result_df.at[idx, "gamma"] = corrected_gamma
        result_df.at[idx, "gamma_err"] = corrected_gamma_err

    # Recalculate Del rows (OS - SS)
    recalculate_del_rows(result_df)

    return result_df

def recalculate_del_rows(df):
    """
    Recalculate Del rows as OS - SS.
    """
    # Group by centrality, diff_type, diff_bin
    group_cols = ["centrality", "diff_type", "diff_bin"]

    for group_vals, group_df in df.groupby(group_cols):
        os_rows = group_df[group_df["pair_type"] == "OS"]
        ss_rows = group_df[group_df["pair_type"] == "SS"]
        del_rows = group_df[group_df["pair_type"] == "Del"]

        if len(os_rows) == 1 and len(ss_rows) == 1 and len(del_rows) == 1:
            os_row = os_rows.iloc[0]
            ss_row = ss_rows.iloc[0]
            del_idx = del_rows.index[0]

            # Calculate Del = OS - SS
            del_delta = os_row["delta"] - ss_row["delta"]
            del_delta_err = np.sqrt(os_row["delta_err"]**2 + ss_row["delta_err"]**2)

            del_gamma = os_row["gamma"] - ss_row["gamma"]
            del_gamma_err = np.sqrt(os_row["gamma_err"]**2 + ss_row["gamma_err"]**2)

            # Update Del row
            df.at[del_idx, "delta"] = del_delta
            df.at[del_idx, "delta_err"] = del_delta_err
            df.at[del_idx, "gamma"] = del_gamma
            df.at[del_idx, "gamma_err"] = del_gamma_err

def get_task_names(proton_dir="../dataset_merger/Merged/Proton"):
    """
    Get all task names from Proton finalise files.
    """
    pattern = f"{proton_dir}/finalise_*.csv"
    files = glob.glob(pattern)
    task_names = []

    for file_path in files:
        filename = os.path.basename(file_path)
        task_name = filename.replace('finalise_', '').replace('.csv', '')
        task_names.append(task_name)

    return task_names

def process_all_tasks(proton_dir="../dataset_merger/Merged/Proton",
                     lambda_dir="../dataset_merger/Merged/Lambda",
                     output_dir="./Proton",
                     specific_task=None):
    """
    Process all tasks or a specific task.
    """
    # 获取任务名称
    if specific_task:
        task_names = [specific_task]
        print(f"Processing specific task: {specific_task}")
    else:
        task_names = get_task_names(proton_dir)
        print(f"Found {len(task_names)} tasks: {task_names}")

    if not task_names:
        print("No tasks found!")
        return

    # 创建输出目录
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    successful_processes = 0

    for task in task_names:
        print(f"\n  Processing task: {task}")

        # 从任务名称中提取中心度
        try:
            centrality = int(task.split('_')[0])  # 假设任务名称格式为 "centX_..."
        except (ValueError, IndexError):
            print(f"    Warning: Could not extract centrality from task name {task}, using default")
            centrality = 25  # 默认使用25%中心度

        # 计算该中心度的fp
        fp, fp_err = calculate_fp(centrality)
        print(f"    Using fp = {fp:.6f} ± {fp_err:.6f} for centrality {centrality}%")

        # 加载Proton数据
        proton_file = f"{proton_dir}/finalise_{task}.csv"
        if not os.path.exists(proton_file):
            print(f"    Error: Proton file not found: {proton_file}")
            continue

        proton_df = pd.read_csv(proton_file)
        print(f"    Loaded Proton data: {len(proton_df)} rows")

        # 加载Lambda数据
        lambda_df = load_lambda_data(task, lambda_dir)

        if lambda_df is None:
            print(f"    Skipping {task} due to missing Lambda data")
            continue

        # 处理数据
        result_df = process_proton_data(proton_df, lambda_df, fp, fp_err)

        if result_df is not None:
            # 保存结果
            output_file = f"{output_dir}/finalise_feeddown_dispose_{task}.csv"
            try:
                result_df.to_csv(output_file, index=False)
                print(f"    ✓ Saved: {output_file}")
                successful_processes += 1
            except Exception as e:
                print(f"    Error saving {output_file}: {e}")
        else:
            print(f"    Failed to process {task}")

    print(f"\nCompleted: {successful_processes}/{len(task_names)} tasks processed successfully")

def main():
    """
    Main function with command line argument support.
    """
    parser = argparse.ArgumentParser(
        description="Apply feeddown correction to Proton data using Lambda secondary data",
        epilog="Example: python feeddown_dispose.py -t default"
    )
    parser.add_argument(
        "--task", "-t",
        type=str,
        help="Specify task name to process. If not provided, all tasks will be processed."
    )
    parser.add_argument(
        "--proton-dir",
        type=str,
        default="../dataset_merger/Merged/Proton",
        help="Directory containing Proton data files"
    )
    parser.add_argument(
        "--lambda-dir",
        type=str,
        default="../dataset_merger/Merged/Lambda",
        help="Directory containing Lambda data files"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./Proton",
        help="Output directory for processed files"
    )
    parser.add_argument(
        "--ratio-lpsec-ll",
        type=float,
        default=0.9033,
        help="Ratio between O_lpsec and O_ll (default: 0.9033)"
    )
    parser.add_argument(
        "--ratio-lpsec-ll-err",
        type=float,
        default=0.0436,
        help="Error of ratio between O_lpsec and O_ll (default: 0.0436)"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="Feeddown Dispose 1.0.0"
    )

    args = parser.parse_args()

    print("Feeddown Dispose Script")
    print("=" * 50)
    print(f"Using O_lpsec/O_ll ratio: {args.ratio_lpsec_ll} ± {args.ratio_lpsec_ll_err}")

    process_all_tasks(
        proton_dir=args.proton_dir,
        lambda_dir=args.lambda_dir,
        output_dir=args.output_dir,
        specific_task=args.task
    )

    print("\n" + "=" * 50)
    print("Feeddown correction completed!")

if __name__ == "__main__":
    main()
