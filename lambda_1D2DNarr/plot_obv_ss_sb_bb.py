import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import seaborn as sns
import argparse

# 设置样式
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 集中管理所有方法/文件/区间的映射
METHOD_FILES = [
    {'filename': './LHC18q/Lambda/fit_obvs_2Dfit_default.csv', 'method': '2D Fit', 'range': 'Default'},
    {'filename': './LHC18q/Lambda/fit_obvs_1Dfit_default.csv', 'method': '1D Fit', 'range': 'Default'},
    {'filename': './LHC18q/Lambda/fit_obvs_narrMass_default.csv', 'method': 'Narrow Mass', 'range': 'Default'},
    {'filename': './outputs_scan/fit_obvs_default_narrMass_range0p020001.csv', 'method': 'Narrow Mass Range', 'range': '0.020002'},
    # 可继续添加其它range
]

# 图表中使用的颜色方案
COLORS = {
    '2D Fit': '#e74c3c',
    '1D Fit': '#3498db',
    'Narrow Mass': '#2ecc71',
    'Narrow Mass Range': '#9b59b6'
}

# 定义不同中心度下的 fs/fb 比例
FS_FB_RATIOS = {
    "00-10": 1.6,
    "10-20": 2.25,
    "20-30": 3.4,
    "30-40": 4.4,
    "40-50": 5.5,
    "50-60": 6.0
}

def load_all_data():
    """统一加载所有方法/区间的数据，生成标准化DataFrame"""
    all_data = []
    for entry in METHOD_FILES:
        try:
            df = pd.read_csv(entry['filename'])
            df['method'] = entry['method']
            df['range'] = entry['range']
            # 筛选中心度 < 60的数据
            df = df[df['centrality'] < 60].copy()
            # 对于1D fit，将sb重新映射为bb
            if entry['method'] == '1D Fit':
                df['delta_bb'] = df['delta_sb'].copy()
                df['rawgamma_bb'] = df['rawgamma_sb'].copy()
                df['delta_bb_err'] = df['delta_sb_err'].copy()
                df['rawgamma_bb_err'] = df['rawgamma_sb_err'].copy()
            # 如果是Narrow Mass Range，确保所有必要的列都存在
            if entry['method'] == 'Narrow Mass Range':
                if 'delta_bb' not in df.columns:
                    df['delta_bb'] = df['delta_sb'].copy()
                    df['delta_bb_err'] = df['delta_sb_err'].copy()
                if 'rawgamma_bb' not in df.columns:
                    df['rawgamma_bb'] = df['rawgamma_sb'].copy()
                    df['rawgamma_bb_err'] = df['rawgamma_sb_err'].copy()
            all_data.append(df)
            print(f"Successfully loaded {entry['filename']}: {len(df)} rows")
        except FileNotFoundError:
            print(f"Warning: File {entry['filename']} not found")
            continue
    if not all_data:
        raise FileNotFoundError("No data files found! Please ensure files are in current directory.")
    # 合并所有数据
    combined_data = pd.concat(all_data, ignore_index=True)
    return combined_data

def create_delta_plots(data, ss_only=False):
    """Create Delta parameters comparison plots"""
    # Get unique pair_types
    pair_types = data['pair_type'].unique()
    methods = data['method'].unique()

    # 使用全局定义的颜色方案
    colors = COLORS

    # Create 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Delta Parameters vs Centrality: Method Comparison',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    for i, pair_type in enumerate(pair_types):
        if i >= 4:  # Only process first 4 pair_types
            break

        ax = axes[i]

        # Filter data for current pair_type
        pair_data = data[data['pair_type'] == pair_type]

        # Plot data for each method
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue

            centrality = method_data['centrality']
            color = colors.get(method, '#95a5a6')

            # 为不同方法选择不同的线型
            linestyle = '-'
            if method == 'Narrow Mass Range':
                linestyle = '--'
            elif method == 'Narrow Mass':
                linestyle = '-.'

            # Plot delta parameters with error bars
            ax.errorbar(centrality, method_data['delta_ss'], yerr=method_data['delta_ss_err'],
                   fmt='o'+linestyle, color=color, alpha=0.8, linewidth=2.5, markersize=7,
                   label=f'{method} - δ_ss', capsize=4)

            # Only plot sb and bb if ss_only is False
            if not ss_only:
                # For 1D fit, plot bb instead of sb
                if method == '1D Fit':
                    ax.errorbar(centrality, method_data['delta_bb'], yerr=method_data['delta_bb_err'],
                           fmt='s--', color=color, alpha=0.8, linewidth=2.5, markersize=7,
                           label=f'{method} - δ_bb', capsize=4)
                else:
                    ax.errorbar(centrality, method_data['delta_sb'], yerr=method_data['delta_sb_err'],
                           fmt='s--', color=color, alpha=0.8, linewidth=2.5, markersize=7,
                           label=f'{method} - δ_sb', capsize=4)

                # Plot delta_bb if not all zeros (for 2D and Narrow Mass)
                if method != '1D Fit' and not (method_data['delta_bb'] == 0).all():
                    ax.errorbar(centrality, method_data['delta_bb'], yerr=method_data['delta_bb_err'],
                           fmt='^:', color=color, alpha=0.8, linewidth=2, markersize=6,
                           label=f'{method} - δ_bb', capsize=4)

        # Set subplot properties
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Delta Parameter Value', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=0.8)

        # Set legend
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    output_filename = 'delta_parameters_comparison_ss_only.png' if ss_only else 'delta_parameters_comparison.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

def create_gamma_plots(data, ss_only=False):
    """Create Gamma parameters comparison plots"""
    # Get unique pair_types
    pair_types = data['pair_type'].unique()
    methods = data['method'].unique()

    # Use global color scheme
    colors = COLORS

    # Create 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Gamma Parameters vs Centrality: Method Comparison',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    for i, pair_type in enumerate(pair_types):
        if i >= 4:  # Only process first 4 pair_types
            break

        ax = axes[i]

        # Filter data for current pair_type
        pair_data = data[data['pair_type'] == pair_type]

        # Plot data for each method
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue

            centrality = method_data['centrality']
            color = colors.get(method, '#95a5a6')

            # 为不同方法选择不同的线型
            linestyle = '-'
            if method == 'Narrow Mass Range':
                linestyle = '--'
            elif method == 'Narrow Mass':
                linestyle = '-.'

            # Plot gamma parameters
            ax.plot(centrality, method_data['rawgamma_ss'],
                   'o'+linestyle, color=color, alpha=0.8, linewidth=2.5, markersize=7,
                   label=f'{method} - γ_ss')

            # Only plot sb and bb if ss_only is False
            if not ss_only:
                # For 1D fit, plot bb instead of sb
                if method == '1D Fit':
                    ax.plot(centrality, method_data['rawgamma_bb'],
                           's--', color=color, alpha=0.8, linewidth=2.5, markersize=7,
                           label=f'{method} - γ_bb')
                else:
                    ax.plot(centrality, method_data['rawgamma_sb'],
                           's--', color=color, alpha=0.8, linewidth=2.5, markersize=7,
                           label=f'{method} - γ_sb')

                # Plot rawgamma_bb if not all zeros (for 2D and Narrow Mass)
                if method != '1D Fit' and not (method_data['rawgamma_bb'] == 0).all():
                    ax.plot(centrality, method_data['rawgamma_bb'],
                           '^:', color=color, alpha=0.8, linewidth=2, markersize=6,
                           label=f'{method} - γ_bb')

        # Set subplot properties
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Gamma Parameter Value', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=0.8)

        # Set legend
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    output_filename = 'gamma_parameters_comparison_ss_only.png' if ss_only else 'gamma_parameters_comparison.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.show()

def create_delta_ratio_plots(data):
    """Create Delta SS ratio plots: (2D or 1D) / Narrow Mass and Narrow Mass Range"""
    pair_types = data['pair_type'].unique()

    # Create 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Delta SS Parameter Ratios vs Centrality: Methods / Reference',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    for i, pair_type in enumerate(pair_types):
        if i >= 4:
            break

        ax = axes[i]
        pair_data = data[data['pair_type'] == pair_type]

        # Get data for each method
        narrow_data = pair_data[pair_data['method'] == 'Narrow Mass']
        narrow_range_data = pair_data[pair_data['method'] == 'Narrow Mass Range']
        twod_data = pair_data[pair_data['method'] == '2D Fit']
        oned_data = pair_data[pair_data['method'] == '1D Fit']

        # Calculate ratios for 2D/Narrow Mass
        if len(twod_data) > 0 and len(narrow_data) > 0:
            # Merge on centrality to ensure matching points
            merged_2d = pd.merge(twod_data, narrow_data, on='centrality', suffixes=('_2d', '_narrow'))
            if len(merged_2d) > 0:
                ratio_2d = merged_2d['delta_ss_2d'] / merged_2d['delta_ss_narrow']
                # Handle division by zero or near zero
                ratio_2d = ratio_2d.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_2d['centrality'], ratio_2d,
                       'o-', color='#e74c3c', alpha=0.8, linewidth=2.5, markersize=7,
                       label='2D Fit / Narrow Mass')

        # Calculate ratios for 2D/Narrow Mass Range
        if len(twod_data) > 0 and len(narrow_range_data) > 0:
            # Merge on centrality to ensure matching points
            merged_2d_range = pd.merge(twod_data, narrow_range_data, on='centrality', suffixes=('_2d', '_narrow_range'))
            if len(merged_2d_range) > 0:
                ratio_2d_range = merged_2d_range['delta_ss_2d'] / merged_2d_range['delta_ss_narrow_range']
                # Handle division by zero or near zero
                ratio_2d_range = ratio_2d_range.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_2d_range['centrality'], ratio_2d_range,
                       'o--', color='#9b59b6', alpha=0.8, linewidth=2.5, markersize=7,
                       label='2D Fit / Narrow Mass Range')

        # Calculate ratios for 1D/Narrow Mass
        if len(oned_data) > 0 and len(narrow_data) > 0:
            merged_1d = pd.merge(oned_data, narrow_data, on='centrality', suffixes=('_1d', '_narrow'))
            if len(merged_1d) > 0:
                ratio_1d = merged_1d['delta_ss_1d'] / merged_1d['delta_ss_narrow']
                ratio_1d = ratio_1d.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_1d['centrality'], ratio_1d,
                       's--', color='#3498db', alpha=0.8, linewidth=2.5, markersize=7,
                       label='1D Fit / Narrow Mass')

        # Calculate ratios for 1D/Narrow Mass Range
        if len(oned_data) > 0 and len(narrow_range_data) > 0:
            merged_1d_range = pd.merge(oned_data, narrow_range_data, on='centrality', suffixes=('_1d', '_narrow_range'))
            if len(merged_1d_range) > 0:
                ratio_1d_range = merged_1d_range['delta_ss_1d'] / merged_1d_range['delta_ss_narrow_range']
                ratio_1d_range = ratio_1d_range.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_1d_range['centrality'], ratio_1d_range,
                       's-.', color='#f39c12', alpha=0.8, linewidth=2.5, markersize=7,
                       label='1D Fit / Narrow Mass Range')

        # Set subplot properties
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Delta SS Ratio', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=1, color='black', linestyle='-', alpha=0.5, linewidth=1)

        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig('delta_ratio_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_gamma_ratio_plots(data):
    """Create Gamma SS ratio plots: (2D or 1D) / Narrow Mass and Narrow Mass Range"""
    pair_types = data['pair_type'].unique()

    # Create 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Gamma SS Parameter Ratios vs Centrality: Methods / Reference',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    for i, pair_type in enumerate(pair_types):
        if i >= 4:
            break

        ax = axes[i]
        pair_data = data[data['pair_type'] == pair_type]

        # Get data for each method
        narrow_data = pair_data[pair_data['method'] == 'Narrow Mass']
        narrow_range_data = pair_data[pair_data['method'] == 'Narrow Mass Range']
        twod_data = pair_data[pair_data['method'] == '2D Fit']
        oned_data = pair_data[pair_data['method'] == '1D Fit']

        # Calculate ratios for 2D/Narrow Mass
        if len(twod_data) > 0 and len(narrow_data) > 0:
            # Merge on centrality to ensure matching points
            merged_2d = pd.merge(twod_data, narrow_data, on='centrality', suffixes=('_2d', '_narrow'))
            if len(merged_2d) > 0:
                ratio_2d = merged_2d['rawgamma_ss_2d'] / merged_2d['rawgamma_ss_narrow']
                # Handle division by zero or near zero
                ratio_2d = ratio_2d.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_2d['centrality'], ratio_2d,
                       'o-', color='#e74c3c', alpha=0.8, linewidth=2.5, markersize=7,
                       label='2D Fit / Narrow Mass')

        # Calculate ratios for 2D/Narrow Mass Range
        if len(twod_data) > 0 and len(narrow_range_data) > 0:
            # Merge on centrality to ensure matching points
            merged_2d_range = pd.merge(twod_data, narrow_range_data, on='centrality', suffixes=('_2d', '_narrow_range'))
            if len(merged_2d_range) > 0:
                ratio_2d_range = merged_2d_range['rawgamma_ss_2d'] / merged_2d_range['rawgamma_ss_narrow_range']
                # Handle division by zero or near zero
                ratio_2d_range = ratio_2d_range.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_2d_range['centrality'], ratio_2d_range,
                       'o--', color='#9b59b6', alpha=0.8, linewidth=2.5, markersize=7,
                       label='2D Fit / Narrow Mass Range')

        # Calculate ratios for 1D/Narrow Mass
        if len(oned_data) > 0 and len(narrow_data) > 0:
            merged_1d = pd.merge(oned_data, narrow_data, on='centrality', suffixes=('_1d', '_narrow'))
            if len(merged_1d) > 0:
                ratio_1d = merged_1d['rawgamma_ss_1d'] / merged_1d['rawgamma_ss_narrow']
                ratio_1d = ratio_1d.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_1d['centrality'], ratio_1d,
                       's--', color='#3498db', alpha=0.8, linewidth=2.5, markersize=7,
                       label='1D Fit / Narrow Mass')

        # Calculate ratios for 1D/Narrow Mass Range
        if len(oned_data) > 0 and len(narrow_range_data) > 0:
            merged_1d_range = pd.merge(oned_data, narrow_range_data, on='centrality', suffixes=('_1d', '_narrow_range'))
            if len(merged_1d_range) > 0:
                ratio_1d_range = merged_1d_range['rawgamma_ss_1d'] / merged_1d_range['rawgamma_ss_narrow_range']
                ratio_1d_range = ratio_1d_range.replace([np.inf, -np.inf], np.nan)
                ax.plot(merged_1d_range['centrality'], ratio_1d_range,
                       's-.', color='#f39c12', alpha=0.8, linewidth=2.5, markersize=7,
                       label='1D Fit / Narrow Mass Range')

        # Set subplot properties
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Gamma SS Ratio', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=1, color='black', linestyle='-', alpha=0.5, linewidth=1)

        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig('gamma_ratio_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_method_comparison_summary(data):
    """Create method comparison summary plots"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Method Comparison Summary: Parameter Averages', fontsize=18, fontweight='bold')

    methods = data['method'].unique()
    colors = ['#e74c3c', '#3498db', '#2ecc71']

    parameters = [
        ('delta_ss', 'Delta SS'),
        ('delta_sb', 'Delta SB'),
        ('delta_bb', 'Delta BB'),
        ('rawgamma_ss', 'Raw Gamma SS'),
        ('rawgamma_sb', 'Raw Gamma SB'),
        ('rawgamma_bb', 'Raw Gamma BB')
    ]

    for idx, (param, title) in enumerate(parameters):
        ax = axes[idx // 3, idx % 3]

        for i, method in enumerate(methods):
            method_data = data[data['method'] == method]
            if len(method_data) == 0:
                continue

            # Group by pair_type and calculate average
            avg_by_pair = method_data.groupby('pair_type')[param].mean()

            x_pos = np.arange(len(avg_by_pair)) + i * 0.25
            bars = ax.bar(x_pos, avg_by_pair.values, width=0.25,
                         color=colors[i], alpha=0.8, label=method)

            # Add value labels
            for bar in bars:
                height = bar.get_height()
                if abs(height) > 1e-10:  # Only show non-zero values
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{height:.4f}', ha='center', va='bottom',
                           fontsize=8, rotation=90)

        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('Pair Type', fontsize=10)
        ax.set_ylabel('Average Value', fontsize=10)
        ax.set_xticks(np.arange(len(avg_by_pair)) + 0.25)
        ax.set_xticklabels(avg_by_pair.index, rotation=45, ha='right')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5)

    plt.tight_layout()
    plt.savefig('method_comparison_summary.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_centrality_trend_analysis(data):
    """Create centrality trend analysis plots"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle('Centrality Trend Analysis: Parameter Changes Across Methods', fontsize=16, fontweight='bold')

    pair_types = data['pair_type'].unique()[:4]
    methods = data['method'].unique()

    colors = {'2D Fit': '#e74c3c', '1D Fit': '#3498db', 'Narrow Mass': '#2ecc71'}

    for i, pair_type in enumerate(pair_types):
        ax = axes[i//2, i%2]

        pair_data = data[data['pair_type'] == pair_type]

        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue

            centrality = method_data['centrality']

            # Calculate delta_ss and delta_sb difference as trend indicator
            if method == '1D Fit':
                trend_indicator = method_data['delta_ss'] - method_data['delta_bb']
            else:
                trend_indicator = method_data['delta_ss'] - method_data['delta_sb']

            ax.plot(centrality, trend_indicator, 'o-',
                   color=colors.get(method, '#95a5a6'),
                   linewidth=2.5, markersize=7, alpha=0.8,
                   label=f'{method}')

        ax.set_title(f'{pair_type}\n(δ_ss - δ_sb/bb Trend)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Centrality (%)', fontsize=11)
        ax.set_ylabel('Trend Indicator', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('centrality_trend_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def print_data_summary(data):
    """Print data summary"""
    print("\n" + "="*50)
    print("Data Summary")
    print("="*50)

    print(f"Total data rows: {len(data)}")
    print(f"Method types: {', '.join(data['method'].unique())}")
    print(f"Pair types: {', '.join(data['pair_type'].unique())}")
    print(f"Centrality range: {data['centrality'].min():.1f} - {data['centrality'].max():.1f}")

    print("\nData distribution by method:")
    method_counts = data['method'].value_counts()
    for method, count in method_counts.items():
        print(f"  {method}: {count} rows")

    print("\nData distribution by pair type:")
    pair_counts = data['pair_type'].value_counts()
    for pair_type, count in pair_counts.items():
        print(f"  {pair_type}: {count} rows")

def create_weighted_delta_plots(data):
    """创建加权的delta组件图 (delta_ss * fsfs + delta_bs * 2*fsfb + delta_bb * fbfb)"""
    # 分析2D Fit方法和Narrow Mass Range方法
    method_data = data[(data['method'] == '2D Fit') | (data['method'] == 'Narrow Mass Range')].copy()

    if len(method_data) == 0:
        print("警告: 没有找到2D Fit或Narrow Mass Range方法的数据")
        return

    # 获取唯一的粒子对类型
    pair_types = method_data['pair_type'].unique()

    # 创建2x2子图
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Weighted Delta Components by Centrality',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    # 定义组件颜色
    component_colors = {
        'fs²·δ_ss': '#e74c3c',
        '2fsfb·δ_sb': '#3498db',
        'fb²·δ_bb': '#2ecc71',
        'Weighted Sum 2D': '#8e44ad',
        'Weighted Sum NMR': '#ff7f0e'  # 为Narrow Mass Range的加权和使用不同颜色
    }

    # 定义方法线型
    method_linestyles = {
        '2D Fit': '-',
        'Narrow Mass Range': '--'
    }

    # 将中心度范围映射到实际的fs/fb比例
    def map_centrality_to_ratio(cent):
        if cent < 10:
            return "00-10"
        elif cent < 20:
            return "10-20"
        elif cent < 30:
            return "20-30"
        elif cent < 40:
            return "30-40"
        elif cent < 50:
            return "40-50"
        else:
            return "50-60"

    # 计算fs和fb的值
    def calculate_fs_fb(cent_range):
        ratio = FS_FB_RATIOS[cent_range]
        fs = ratio / (1 + ratio)
        fb = 1 / (1 + ratio)
        return fs, fb

    # 存储所有处理过的数据，用于创建堆叠图
    processed_data = {}

    for i, pair_type in enumerate(pair_types):
        if i >= 4:  # 只处理前4种粒子对类型
            break

        ax = axes[i]
        pair_processed_data = {}

        # 为每个方法单独处理
        for method in ['2D Fit', 'Narrow Mass Range']:
            # 筛选当前粒子对类型和方法的数据
            method_pair_data = method_data[(method_data['pair_type'] == pair_type) &
                                          (method_data['method'] == method)].copy()

            if len(method_pair_data) == 0:
                print(f"警告: 未找到 {method} 方法的 {pair_type} 数据，跳过此组合")
                continue

            # 按中心度排序
            method_pair_data = method_pair_data.sort_values('centrality')

            # 计算每个中心度下的fs和fb
            fs_values = []
            fb_values = []
            for cent in method_pair_data['centrality']:
                cent_range = map_centrality_to_ratio(cent)
                fs, fb = calculate_fs_fb(cent_range)
                fs_values.append(fs)
                fb_values.append(fb)

            method_pair_data['fs'] = fs_values
            method_pair_data['fb'] = fb_values

            # 计算加权组件
            method_pair_data['weighted_ss'] = method_pair_data['delta_ss'] * (method_pair_data['fs'] ** 2)
            method_pair_data['weighted_sb'] = method_pair_data['delta_sb'] * (2 * method_pair_data['fs'] * method_pair_data['fb'])
            method_pair_data['weighted_bb'] = method_pair_data['delta_bb'] * (method_pair_data['fb'] ** 2)
            method_pair_data['weighted_sum'] = method_pair_data['weighted_ss'] + method_pair_data['weighted_sb'] + method_pair_data['weighted_bb']

            # 保存处理后的数据
            pair_processed_data[method] = method_pair_data

            # 选择线型
            linestyle = method_linestyles[method]

            # 如果是2D Fit，显示所有组件
            if method == '2D Fit':
                # 绘制各个加权组件
                ax.plot(method_pair_data['centrality'], method_pair_data['weighted_ss'],
                       'o'+linestyle, color=component_colors['fs²·δ_ss'], alpha=0.8, linewidth=2.5, markersize=7,
                       label=f'{method} fs²·δ_ss')

                ax.plot(method_pair_data['centrality'], method_pair_data['weighted_sb'],
                       's'+linestyle, color=component_colors['2fsfb·δ_sb'], alpha=0.8, linewidth=2.5, markersize=7,
                       label=f'{method} 2fsfb·δ_sb')

                ax.plot(method_pair_data['centrality'], method_pair_data['weighted_bb'],
                       '^'+linestyle, color=component_colors['fb²·δ_bb'], alpha=0.8, linewidth=2.5, markersize=7,
                       label=f'{method} fb²·δ_bb')

                # 绘制加权和
                ax.plot(method_pair_data['centrality'], method_pair_data['weighted_sum'],
                       'D'+linestyle, color=component_colors['Weighted Sum 2D'], alpha=1.0, linewidth=3, markersize=9,
                       label=f'{method} Weighted Sum')
            # 如果是Narrow Mass Range，只显示加权和
            else:
                # 只绘制加权和
                ax.plot(method_pair_data['centrality'], method_pair_data['weighted_sum'],
                       'D'+linestyle, color=component_colors['Weighted Sum NMR'], alpha=1.0, linewidth=3, markersize=9,
                       label=f'{method} Weighted Sum')

        # 保存当前粒子对类型的处理数据
        processed_data[pair_type] = pair_processed_data

        # 设置子图属性
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Weighted Delta Value', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=0.8)

        # 设置图例
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    # 调整布局
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig('weighted_delta_components.png', dpi=300, bbox_inches='tight')
    plt.show()

    # 创建堆叠面积图 - 只对2D Fit方法
    if processed_data:
        create_stacked_area_plots(processed_data)
    else:
        print("警告: 未能创建堆叠图，没有有效的处理数据")

    # 返回处理后的数据
    return processed_data

def create_stacked_area_plots(processed_data):
    """创建堆叠面积图，显示不同组分的贡献"""
    if not processed_data or all(not pair_data for pair_data in processed_data.values()):
        print("警告: 无法创建堆叠图，缺少处理后的数据")
        return

    # 创建2x2子图
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Stacked Delta Components by Centrality (2D Fit)',
                 fontsize=20, fontweight='bold', y=0.98)

    axes = axes.flatten()

    # 定义组件颜色 - 保持与线图相同的颜色方案
    component_colors = {
        'fs²·δ_ss': '#e74c3c',
        '2fsfb·δ_sb': '#3498db',
        'fb²·δ_bb': '#2ecc71'
    }

    for i, pair_type in enumerate(processed_data.keys()):
        if i >= 4:  # 只处理前4种粒子对类型
            break

        ax = axes[i]

        # 获取2D Fit方法的数据
        if '2D Fit' not in processed_data[pair_type]:
            # 如果没有2D Fit数据，跳过此粒子对类型
            ax.text(0.5, 0.5, f"No 2D Fit data for {pair_type}",
                   horizontalalignment='center', verticalalignment='center',
                   transform=ax.transAxes, fontsize=12)
            continue

        method_data = processed_data[pair_type]['2D Fit']

        # 准备堆叠图数据
        x = method_data['centrality'].values
        y1 = method_data['weighted_ss'].values  # fs²·δ_ss
        y2 = method_data['weighted_sb'].values  # 2fsfb·δ_sb
        y3 = method_data['weighted_bb'].values  # fb²·δ_bb

        # 创建堆叠面积图
        ax.fill_between(x, 0, y1, alpha=0.7, color=component_colors['fs²·δ_ss'], label='fs²·δ_ss')
        ax.fill_between(x, y1, y1+y2, alpha=0.7, color=component_colors['2fsfb·δ_sb'], label='2fsfb·δ_sb')
        ax.fill_between(x, y1+y2, y1+y2+y3, alpha=0.7, color=component_colors['fb²·δ_bb'], label='fb²·δ_bb')

        # 绘制总和线
        ax.plot(x, y1+y2+y3, 'k-', linewidth=2.5, label='Total')

        # 如果有Narrow Mass Range数据，也绘制其总和线
        if 'Narrow Mass Range' in processed_data[pair_type]:
            nmr_data = processed_data[pair_type]['Narrow Mass Range']
            if len(nmr_data) > 0:
                ax.plot(nmr_data['centrality'], nmr_data['weighted_sum'],
                      'k--', linewidth=2.5, label='Narrow Mass Range Total')

        # 设置子图属性
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Centrality (%)', fontsize=12)
        ax.set_ylabel('Weighted Delta Components', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=0.8)

        # 设置图例
        if i == 0:
            ax.legend(loc='best', fontsize=9)

    # 调整布局
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig('stacked_weighted_delta_components.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """主函数"""
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='绘制参数对比图表')
    parser.add_argument('--ss-only', action='store_true',
                        help='只绘制SS图形，跳过SB和BB图形')
    args = parser.parse_args()

    try:
        print("正在加载数据...")
        data = load_all_data()
        # 检查是否成功加载了Narrow Mass Range数据
        has_nmr = 'Narrow Mass Range' in data['method'].unique()
        if not has_nmr:
            print("警告: 未找到Narrow Mass Range方法的数据，部分分析可能不完整")
            print("提示: 请确保'fit_obvs_default_narrMass_range0p020002.csv'文件存在")
        print_data_summary(data)
        print("\n正在生成Delta参数对比图表...")
        create_delta_plots(data, ss_only=args.ss_only)
        print("正在生成Gamma参数对比图表...")
        create_gamma_plots(data, ss_only=args.ss_only)
        print("正在生成Delta SS比值分析图表...")
        create_delta_ratio_plots(data)
        print("正在生成Gamma SS比值分析图表...")
        create_gamma_ratio_plots(data)
        try:
            print("正在生成方法对比总结...")
            create_method_comparison_summary(data)
            print("✓ 已生成方法对比总结")
        except Exception as e:
            print(f"警告: 无法生成方法对比总结: {e}")
        try:
            print("正在生成中心度趋势分析...")
            create_centrality_trend_analysis(data)
            print("✓ 已生成中心度趋势分析")
        except Exception as e:
            print(f"警告: 无法生成中心度趋势分析: {e}")
        try:
            print("正在生成加权Delta组件图表...")
            weighted_data = create_weighted_delta_plots(data)
            if weighted_data:
                print("✓ 加权分析完成，成功生成加权组件图")
            else:
                print("警告: 加权分析未生成有效数据")
        except Exception as e:
            print(f"警告: 无法生成加权Delta组件图表: {e}")
            print("提示: 确保至少有一个2D Fit或Narrow Mass Range方法的数据文件")
        print("\n分析完成！生成的图片文件:")
        if args.ss_only:
            print("- delta_parameters_comparison_ss_only.png")
            print("- gamma_parameters_comparison_ss_only.png")
        else:
            print("- delta_parameters_comparison.png")
            print("- gamma_parameters_comparison.png")
        print("- delta_ratio_analysis.png")
        print("- gamma_ratio_analysis.png")
        # 检查文件是否存在
        import os
        if os.path.exists("method_comparison_summary.png"):
            print("- method_comparison_summary.png")
        if os.path.exists("centrality_trend_analysis.png"):
            print("- centrality_trend_analysis.png")
        if os.path.exists("weighted_delta_components.png"):
            print("- weighted_delta_components.png")
        if os.path.exists("stacked_weighted_delta_components.png"):
            print("- stacked_weighted_delta_components.png")
        if has_nmr:
            print("(注: 包含了Narrow Mass Range方法的分析结果)")
        else:
            print("(注: 未检测到Narrow Mass Range方法的数据，部分高级分析可能不完整)")
    except Exception as e:
        print(f"错误: {e}")
        print("\n请确保以下文件在当前目录下:")
        print("必需文件:")
        print("- fit_obvs_2Dfit_default.csv")
        print("- fit_obvs_1Dfit_default.csv")
        print("- fit_obvs_narrMass_default.csv")
        print("\n可选文件:")
        print("- fit_obvs_default_narrMass_range0p020002.csv (用于加权组件分析和Narrow Mass Range比较)")
        import traceback
        print("\n详细错误信息:")
        traceback.print_exc()

if __name__ == "__main__":
    main()
