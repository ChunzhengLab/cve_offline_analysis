#!/usr/bin/env python3
"""
Plotting script for LHC particle correlation data.
This script generates various plots for Lambda-Lambda, Lambda-Proton,
Lambda-Hadron, and Lambda-Pion correlations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from pathlib import Path

# Set up modern style for plots
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif'],
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.figsize': (10, 8),
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 2.0,
    'errorbar.capsize': 3,
})

# Define color palettes (modern and colorblind-friendly)
COLORS = {
    'merged': '#1f77b4',  # blue
    '18q': '#ff7f0e',     # orange
    '18r': '#2ca02c',     # green
    'feeddown': '#d62728', # red

    'Lambda-Lambda': '#1f77b4',  # blue
    'Lambda-Proton': '#ff7f0e',  # orange
    'Lambda-Hadron': '#2ca02c',  # green
    'Lambda-Pion': '#9467bd',    # purple
}

# Data path definitions
DATA_PATHS = {
    'Lambda-Lambda': {
        'merged': '../dataset_merger/Merged/Lambda/finalise_default.csv',
        '18q': '../obv_fit_lambda2/LHC18q/Lambda/finalise_default.csv',
        '18r': '../obv_fit_lambda2/LHC18r/Lambda/finalise_default.csv'
    },
    'Lambda-Proton': {
        'merged': '../dataset_merger/Merged/Proton/finalise_default.csv',
        '18q': '../obv_fit/LHC18q/Proton/finalise_default.csv',
        '18r': '../obv_fit/LHC18r/Proton/finalise_default.csv',
        'feeddown': '../feeddown_dispose/Proton/finalise_feeddown_dispose_default.csv'
    },
    'Lambda-Hadron': {
        'merged': '../dataset_merger/Merged/Hadron/finalise_default.csv',
        '18q': '../obv_fit/LHC18q/Hadron/finalise_default.csv',
        '18r': '../obv_fit/LHC18r/Hadron/finalise_default.csv'
    },
    'Lambda-Pion': {
        'merged': '../dataset_merger/Merged/Pion/finalise_default.csv',
        '18q': '../obv_fit/LHC18q/Pion/finalise_default.csv',
        '18r': '../obv_fit/LHC18r/Pion/finalise_default.csv'
    }
}

def load_data(file_path, is_lambda_lambda=False):
    """
    Load data from CSV file with appropriate handling of Lambda-Lambda format.

    Args:
        file_path: Path to the CSV file
        is_lambda_lambda: Boolean indicating if data is in Lambda-Lambda format

    Returns:
        DataFrame or None if file not found
    """
    try:
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}")
            return None

        df = pd.read_csv(file_path)

        # Remove resolution column if present
        if 'resolution' in df.columns:
            df = df.drop(columns=['resolution'])

        # Filter centrality range (0-60)
        df = df[df['centrality'] <= 60]

        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def prepare_output_directory(directory):
    """Create output directory if it doesn't exist"""
    Path(directory).mkdir(parents=True, exist_ok=True)
    print(f"Output directory prepared: {directory}")

def plot_particle_comparison(particle_type, value_type, diff_type=None, diff_bin=None):
    """
    Plot comparison of 18q, 18r, and merged data for a given particle type.

    Args:
        particle_type: Type of particle combination (e.g. 'Lambda-Lambda')
        value_type: 'delta' or 'gamma'
        diff_type: Differential type (e.g. 'DEta', 'SPt', 'Intg')
        diff_bin: Differential bin value
    """
    is_lambda_lambda = (particle_type == 'Lambda-Lambda')
    paths = DATA_PATHS[particle_type]

    # Load data
    data = {}
    for source, path in paths.items():
        if source == 'feeddown':  # Skip feeddown for this plot
            continue
        data[source] = load_data(path, is_lambda_lambda)

    # Check if we have any data to plot
    if not any(df is not None and not df.empty for df in data.values()):
        print(f"Skipping {particle_type} plot due to missing data")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"{particle_type} Correlation: {value_type.capitalize()} Comparison", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Plot data
    for source, df in data.items():
        if df is None:
            continue

        if is_lambda_lambda:
            # For Lambda-Lambda, filter by pair_type directly
            ss_data = df[df['pair_type'] == 'SS']
            os_data = df[df['pair_type'] == 'OS']
            del_data = df[df['pair_type'] == 'Del']
        else:
            # For other particle types, also filter by diff_type and diff_bin
            filter_cond = True
            if diff_type:
                filter_cond = (df['diff_type'] == diff_type)
            if diff_bin is not None:
                filter_cond = filter_cond & (df['diff_bin'] == diff_bin)

            ss_data = df[filter_cond & (df['pair_type'] == 'SS')]
            os_data = df[filter_cond & (df['pair_type'] == 'OS')]
            del_data = df[filter_cond & (df['pair_type'] == 'Del')]

        # Plot SS and OS in first subplot
        ax1.errorbar(
            ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[source], alpha=0.7,
            label=f"{source} SS", markerfacecolor='none'
        )

        ax1.errorbar(
            os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
            marker='s', linestyle='--', color=COLORS[source],
            label=f"{source} OS"
        )

        # Plot Del in second subplot
        ax2.errorbar(
            del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[source],
            label=f"{source} Del"
        )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Add grid
    ax1.grid(True, alpha=0.3)
    ax2.grid(True, alpha=0.3)

    # Prepare plot title and filename details
    plot_title = particle_type
    if diff_type:
        plot_title += f" {diff_type}"
    if diff_bin is not None:
        plot_title += f" {diff_bin}"

    # Create output directory
    output_dir = f"./plots/1_period_comparison/{particle_type}"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"{particle_type}_{value_type}"
    if diff_type:
        filename += f"_{diff_type}"
    if diff_bin is not None:
        filename += f"_{diff_bin}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def plot_proton_feeddown_comparison(value_type, diff_type=None, diff_bin=None):
    """
    Plot comparison of merged and feeddown-corrected Lambda-Proton data.

    Args:
        value_type: 'delta' or 'gamma'
        diff_type: Differential type (e.g. 'DEta', 'SPt', 'Intg')
        diff_bin: Differential bin value
    """
    particle_type = 'Lambda-Proton'
    paths = DATA_PATHS[particle_type]

    # Load data
    merged_data = load_data(paths['merged'])
    feeddown_data = load_data(paths['feeddown'])

    if merged_data is None or feeddown_data is None:
        print(f"Skipping feeddown comparison plot due to missing data")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"{particle_type} Correlation: {value_type.capitalize()} with Feeddown Correction", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Apply filters
    filter_cond = True
    if diff_type:
        filter_cond = (merged_data['diff_type'] == diff_type) & (feeddown_data['diff_type'] == diff_type)
    if diff_bin is not None:
        filter_cond = filter_cond & (merged_data['diff_bin'] == diff_bin) & (feeddown_data['diff_bin'] == diff_bin)

    # Filter data
    merged_ss = merged_data[filter_cond & (merged_data['pair_type'] == 'SS')]
    merged_os = merged_data[filter_cond & (merged_data['pair_type'] == 'OS')]
    merged_del = merged_data[filter_cond & (merged_data['pair_type'] == 'Del')]

    feeddown_ss = feeddown_data[filter_cond & (feeddown_data['pair_type'] == 'SS')]
    feeddown_os = feeddown_data[filter_cond & (feeddown_data['pair_type'] == 'OS')]
    feeddown_del = feeddown_data[filter_cond & (feeddown_data['pair_type'] == 'Del')]

    # Plot SS and OS in first subplot
    ax1.errorbar(
        merged_ss['centrality'], merged_ss[value_type], merged_ss[f"{value_type}_err"],
        marker='o', linestyle='-', color=COLORS['merged'], alpha=0.7,
        label="Merged SS", markerfacecolor='none'
    )

    ax1.errorbar(
        merged_os['centrality'], merged_os[value_type], merged_os[f"{value_type}_err"],
        marker='s', linestyle='--', color=COLORS['merged'],
        label="Merged OS"
    )

    ax1.errorbar(
        feeddown_ss['centrality'], feeddown_ss[value_type], feeddown_ss[f"{value_type}_err"],
        marker='o', linestyle='-', color=COLORS['feeddown'], alpha=0.7,
        label="Feeddown-corrected SS", markerfacecolor='none'
    )

    ax1.errorbar(
        feeddown_os['centrality'], feeddown_os[value_type], feeddown_os[f"{value_type}_err"],
        marker='s', linestyle='--', color=COLORS['feeddown'],
        label="Feeddown-corrected OS"
    )

    # Plot Del in second subplot
    ax2.errorbar(
        merged_del['centrality'], merged_del[value_type], merged_del[f"{value_type}_err"],
        marker='o', linestyle='-', color=COLORS['merged'],
        label="Merged Del"
    )

    ax2.errorbar(
        feeddown_del['centrality'], feeddown_del[value_type], feeddown_del[f"{value_type}_err"],
        marker='o', linestyle='-', color=COLORS['feeddown'],
        label="Feeddown-corrected Del"
    )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Prepare plot title and filename details
    plot_title = particle_type
    if diff_type:
        plot_title += f" {diff_type}"
    if diff_bin is not None:
        plot_title += f" {diff_bin}"

    # Create output directory
    output_dir = f"./plots/2_feeddown_comparison"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"{particle_type}_{value_type}_feeddown"
    if diff_type:
        filename += f"_{diff_type}"
    if diff_bin is not None:
        filename += f"_{diff_bin}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def plot_all_particles_comparison(value_type, diff_type='Intg', diff_bin=0.5):
    """
    Plot comparison of all particle types (merged data).

    Args:
        value_type: 'delta' or 'gamma'
        diff_type: Differential type (default: 'Intg')
        diff_bin: Differential bin value (default: 0.5)
    """
    # Load data for each particle type
    data = {}
    for particle_type in DATA_PATHS.keys():
        is_lambda_lambda = (particle_type == 'Lambda-Lambda')
        data_path = DATA_PATHS[particle_type]['merged']
        df = load_data(data_path, is_lambda_lambda)

        if df is not None:
            data[particle_type] = df

    if not data:
        print("No data available for all-particle comparison plot")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"All Particles Comparison: {value_type.capitalize()}", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Plot data for each particle type
    for particle_type, df in data.items():
        is_lambda_lambda = (particle_type == 'Lambda-Lambda')

        if is_lambda_lambda:
            # For Lambda-Lambda, filter by pair_type directly
            ss_data = df[df['pair_type'] == 'SS']
            os_data = df[df['pair_type'] == 'OS']
            del_data = df[df['pair_type'] == 'Del']
        else:
            # For other particle types, filter by diff_type and diff_bin
            filter_cond = (df['diff_type'] == diff_type) & (df['diff_bin'] == diff_bin)
            ss_data = df[filter_cond & (df['pair_type'] == 'SS')]
            os_data = df[filter_cond & (df['pair_type'] == 'OS')]
            del_data = df[filter_cond & (df['pair_type'] == 'Del')]

        # Plot SS and OS in first subplot
        ax1.errorbar(
            ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[particle_type], alpha=0.7,
            label=f"{particle_type} SS", markerfacecolor='none'
        )

        ax1.errorbar(
            os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
            marker='s', linestyle='--', color=COLORS[particle_type],
            label=f"{particle_type} OS"
        )

        # Plot Del in second subplot
        ax2.errorbar(
            del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[particle_type],
            label=f"{particle_type} Del"
        )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Create output directory
    output_dir = f"./plots/3_all_particles"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"all_particles_{value_type}"
    if diff_type and diff_type != 'Intg':
        filename += f"_{diff_type}"
    if diff_bin is not None and diff_bin != 0.5:
        filename += f"_{diff_bin}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def plot_proton_different_spt_comparison(value_type):
    """
    Plot comparison of Lambda-Proton results for different SPt values,
    including Intg for comparison.

    Args:
        value_type: 'delta' or 'gamma'
    """
    particle_type = 'Lambda-Proton'
    data_path = DATA_PATHS[particle_type]['merged']

    # Load data
    df = load_data(data_path)

    if df is None:
        print(f"Skipping SPt comparison plot due to missing data")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"{particle_type}: {value_type.capitalize()} for Different SPt Values", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Define SPt values to plot
    spt_bins = [0.5, 1.5, 2.5]

    # Define custom colors for different SPt values
    spt_colors = {
        'Intg_0.5': '#1f77b4',  # blue
        'SPt_0.5': '#ff7f0e',   # orange
        'SPt_1.5': '#2ca02c',   # green
        'SPt_2.5': '#d62728',   # red
    }

    # Plot data for Intg first
    intg_filter = (df['diff_type'] == 'Intg') & (df['diff_bin'] == 0.5)

    # SS data for Intg
    ss_data = df[intg_filter & (df['pair_type'] == 'SS')]
    ax1.errorbar(
        ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
        marker='o', linestyle='-', color=spt_colors['Intg_0.5'], alpha=0.7,
        label=f"Intg SS", markerfacecolor='none'
    )

    # OS data for Intg
    os_data = df[intg_filter & (df['pair_type'] == 'OS')]
    ax1.errorbar(
        os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
        marker='s', linestyle='--', color=spt_colors['Intg_0.5'],
        label=f"Intg OS"
    )

    # Del data for Intg
    del_data = df[intg_filter & (df['pair_type'] == 'Del')]
    ax2.errorbar(
        del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
        marker='o', linestyle='-', color=spt_colors['Intg_0.5'],
        label=f"Intg Del"
    )

    # Plot data for each SPt value
    for spt_bin in spt_bins:
        spt_filter = (df['diff_type'] == 'SPt') & (df['diff_bin'] == spt_bin)

        # SS data
        ss_data = df[spt_filter & (df['pair_type'] == 'SS')]
        ax1.errorbar(
            ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=spt_colors[f'SPt_{spt_bin}'], alpha=0.7,
            label=f"SPt {spt_bin} SS", markerfacecolor='none'
        )

        # OS data
        os_data = df[spt_filter & (df['pair_type'] == 'OS')]
        ax1.errorbar(
            os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
            marker='s', linestyle='--', color=spt_colors[f'SPt_{spt_bin}'],
            label=f"SPt {spt_bin} OS"
        )

        # Del data
        del_data = df[spt_filter & (df['pair_type'] == 'Del')]
        ax2.errorbar(
            del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=spt_colors[f'SPt_{spt_bin}'],
            label=f"SPt {spt_bin} Del"
        )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Create output directory
    output_dir = f"./plots/5_proton_spt_comparison"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"Lambda_Proton_SPt_comparison_{value_type}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def plot_proton_different_deta_comparison(value_type):
    """
    Plot comparison of Lambda-Proton results for different DEta values,
    including Intg for comparison.

    Args:
        value_type: 'delta' or 'gamma'
    """
    particle_type = 'Lambda-Proton'
    data_path = DATA_PATHS[particle_type]['merged']

    # Load data
    df = load_data(data_path)

    if df is None:
        print(f"Skipping DEta comparison plot due to missing data")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"{particle_type}: {value_type.capitalize()} for Different DEta Values", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Define DEta values to plot
    deta_bins = [0.5, 1.5]

    # Define custom colors for different DEta values
    deta_colors = {
        'Intg_0.5': '#1f77b4',  # blue
        'DEta_0.5': '#ff7f0e',   # orange
        'DEta_1.5': '#2ca02c',   # green
    }

    # Plot data for Intg first
    intg_filter = (df['diff_type'] == 'Intg') & (df['diff_bin'] == 0.5)

    # SS data for Intg
    ss_data = df[intg_filter & (df['pair_type'] == 'SS')]
    ax1.errorbar(
        ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
        marker='o', linestyle='-', color=deta_colors['Intg_0.5'], alpha=0.7,
        label=f"Intg SS", markerfacecolor='none'
    )

    # OS data for Intg
    os_data = df[intg_filter & (df['pair_type'] == 'OS')]
    ax1.errorbar(
        os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
        marker='s', linestyle='--', color=deta_colors['Intg_0.5'],
        label=f"Intg OS"
    )

    # Del data for Intg
    del_data = df[intg_filter & (df['pair_type'] == 'Del')]
    ax2.errorbar(
        del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
        marker='o', linestyle='-', color=deta_colors['Intg_0.5'],
        label=f"Intg Del"
    )

    # Plot data for each DEta value
    for deta_bin in deta_bins:
        deta_filter = (df['diff_type'] == 'DEta') & (df['diff_bin'] == deta_bin)

        # SS data
        ss_data = df[deta_filter & (df['pair_type'] == 'SS')]
        ax1.errorbar(
            ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=deta_colors[f'DEta_{deta_bin}'], alpha=0.7,
            label=f"DEta {deta_bin} SS", markerfacecolor='none'
        )

        # OS data
        os_data = df[deta_filter & (df['pair_type'] == 'OS')]
        ax1.errorbar(
            os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
            marker='s', linestyle='--', color=deta_colors[f'DEta_{deta_bin}'],
            label=f"DEta {deta_bin} OS"
        )

        # Del data
        del_data = df[deta_filter & (df['pair_type'] == 'Del')]
        ax2.errorbar(
            del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=deta_colors[f'DEta_{deta_bin}'],
            label=f"DEta {deta_bin} Del"
        )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Create output directory
    output_dir = f"./plots/6_proton_deta_comparison"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"Lambda_Proton_DEta_comparison_{value_type}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def plot_all_particles_feeddown_proton_comparison(value_type, diff_type='Intg', diff_bin=0.5):
    """
    Plot comparison of all particle types with feeddown-corrected proton data.

    Args:
        value_type: 'delta' or 'gamma'
        diff_type: Differential type (default: 'Intg')
        diff_bin: Differential bin value (default: 0.5)
    """
    # Load data for each particle type
    data = {}
    for particle_type in DATA_PATHS.keys():
        is_lambda_lambda = (particle_type == 'Lambda-Lambda')

        # Use feeddown-corrected data for Lambda-Proton
        if particle_type == 'Lambda-Proton':
            data_path = DATA_PATHS[particle_type]['feeddown']
        else:
            data_path = DATA_PATHS[particle_type]['merged']

        df = load_data(data_path, is_lambda_lambda)

        if df is not None:
            data[particle_type] = df

    if not data:
        print("No data available for feeddown-proton comparison plot")
        return

    # Prepare figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"All Particles with Feeddown-corrected Proton: {value_type.capitalize()}", fontsize=18)

    # Set up axes limits and labels
    for ax in [ax1, ax2]:
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel(f"{value_type.capitalize()} Value", fontsize=14)
        ax.set_xlim(0, 60)

    ax1.set_title("SS and OS", fontsize=16)
    ax2.set_title("Del (OS-SS)", fontsize=16)

    # Plot data for each particle type
    for particle_type, df in data.items():
        is_lambda_lambda = (particle_type == 'Lambda-Lambda')

        # Adjust label for Proton to indicate feeddown correction
        particle_label = f"{particle_type} (Feeddown-corrected)" if particle_type == 'Lambda-Proton' else particle_type

        if is_lambda_lambda:
            # For Lambda-Lambda, filter by pair_type directly
            ss_data = df[df['pair_type'] == 'SS']
            os_data = df[df['pair_type'] == 'OS']
            del_data = df[df['pair_type'] == 'Del']
        else:
            # For other particle types, filter by diff_type and diff_bin
            filter_cond = (df['diff_type'] == diff_type) & (df['diff_bin'] == diff_bin)
            ss_data = df[filter_cond & (df['pair_type'] == 'SS')]
            os_data = df[filter_cond & (df['pair_type'] == 'OS')]
            del_data = df[filter_cond & (df['pair_type'] == 'Del')]

        # Plot SS and OS in first subplot
        ax1.errorbar(
            ss_data['centrality'], ss_data[value_type], ss_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[particle_type], alpha=0.7,
            label=f"{particle_label} SS", markerfacecolor='none'
        )

        ax1.errorbar(
            os_data['centrality'], os_data[value_type], os_data[f"{value_type}_err"],
            marker='s', linestyle='--', color=COLORS[particle_type],
            label=f"{particle_label} OS"
        )

        # Plot Del in second subplot
        ax2.errorbar(
            del_data['centrality'], del_data[value_type], del_data[f"{value_type}_err"],
            marker='o', linestyle='-', color=COLORS[particle_type],
            label=f"{particle_label} Del"
        )

    # Add legends
    ax1.legend(loc='best', frameon=True, framealpha=0.9)
    ax2.legend(loc='best', frameon=True, framealpha=0.9)

    # Create output directory
    output_dir = f"./plots/4_all_particles_feeddown_proton"
    prepare_output_directory(output_dir)

    # Save plot
    filename = f"all_particles_feeddown_proton_{value_type}"
    if diff_type and diff_type != 'Intg':
        filename += f"_{diff_type}"
    if diff_bin is not None and diff_bin != 0.5:
        filename += f"_{diff_bin}"

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_dir}/{filename}.pdf")

def main():
    """
    Main function to generate all plots.
    """
    print("Generating plots for particle correlation data...")

    # 1. Period comparison plots (18q, 18r, merged)
    print("\n1. Generating period comparison plots...")

    # Lambda-Lambda
    for value_type in ['delta', 'gamma']:
        plot_particle_comparison('Lambda-Lambda', value_type)

    # Lambda-Proton (only Intg)
    for value_type in ['delta', 'gamma']:
        # Intg plots
        plot_particle_comparison('Lambda-Proton', value_type, 'Intg', 0.5)

    # Lambda-Hadron
    for value_type in ['delta', 'gamma']:
        plot_particle_comparison('Lambda-Hadron', value_type, 'Intg', 0.5)

    # Lambda-Pion (only if data exists)
    for value_type in ['delta', 'gamma']:
        plot_particle_comparison('Lambda-Pion', value_type, 'Intg', 0.5)

    # 2. Feeddown comparison plots (merged vs. feeddown-corrected)
    print("\n2. Generating feeddown comparison plots...")

    for value_type in ['delta', 'gamma']:
        # Only Intg plots for Lambda-Proton
        plot_proton_feeddown_comparison(value_type, 'Intg', 0.5)

    # 3. All particles comparison (merged data)
    print("\n3. Generating all particles comparison plots...")

    for value_type in ['delta', 'gamma']:
        plot_all_particles_comparison(value_type)

    # 4. All particles with feeddown-corrected proton
    print("\n4. Generating all particles with feeddown-corrected proton plots...")

    for value_type in ['delta', 'gamma']:
        plot_all_particles_feeddown_proton_comparison(value_type)

    # 5. Lambda-Proton different SPt values comparison
    print("\n5. Generating Lambda-Proton different SPt values comparison plots...")

    for value_type in ['delta', 'gamma']:
        plot_proton_different_spt_comparison(value_type)

    # 6. Lambda-Proton different DEta values comparison
    print("\n6. Generating Lambda-Proton different DEta values comparison plots...")

    for value_type in ['delta', 'gamma']:
        plot_proton_different_deta_comparison(value_type)

    print("\nAll plots generated successfully!")

if __name__ == "__main__":
    main()
