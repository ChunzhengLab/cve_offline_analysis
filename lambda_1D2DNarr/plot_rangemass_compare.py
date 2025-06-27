import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import re

# 查找所有扫描结果文件
result_files = glob.glob('outputs_scan_all/fit_obvs_default_*.csv')
print(f"找到的结果文件数量: {len(result_files)}")
print("文件列表:")
for f in result_files:
    print(f"  {os.path.basename(f)}")

# 解析文件名，提取方法和range信息
def parse_filename(filename):
    """解析文件名，提取方法和range_mass"""
    basename = os.path.basename(filename)
    # 匹配格式: fit_obvs_default_<method>_range<value>.csv
    match = re.search(r'fit_obvs_default_(\w+)_range([\d\.p]+)\.csv', basename)
    if match:
        method = match.group(1)
        range_str = match.group(2)
        # 将 0p020000 转换回 0.020000
        range_value = float(range_str.replace('p', '.'))
        return method, range_value
    return None, None

# 读取所有数据并组织
data_dict = {}  # {method: {range_value: dataframe}}

# 读取单独的2D Fit参考文件（如果存在）
reference_2Dfit_file = 'fit_obvs_2Dfit_default.csv'
reference_2Dfit_data = None
if os.path.exists(reference_2Dfit_file):
    reference_2Dfit_data = pd.read_csv(reference_2Dfit_file)
    print(f"读取2D Fit参考数据: {reference_2Dfit_file}")
else:
    print(f"警告: 未找到2D Fit参考文件 {reference_2Dfit_file}")

for file in result_files:
    method, range_value = parse_filename(file)
    if method and range_value is not None:
        if method not in data_dict:
            data_dict[method] = {}
        data_dict[method][range_value] = pd.read_csv(file)
        print(f"读取: {method}, range={range_value:.6f}")

if not data_dict:
    raise FileNotFoundError("未找到任何有效的扫描结果文件")

print(f"\n可用的拟合方法: {list(data_dict.keys())}")

# 检查是否有narrMass方法的数据
if 'narrMass' not in data_dict:
    print("警告: 没有找到Narrow Mass方法的数据，无法计算比值")
    has_narrMass = False
else:
    has_narrMass = True
    print(f"找到Narrow Mass方法的数据，可用的range值: {sorted(data_dict['narrMass'].keys())}")

# 定义粒子对类型、颜色和标记
pair_types = ['LambdaLambda', 'LambdaLambdaBar', 'LambdaBarLambda', 'LambdaBarLambdaBar']
method_colors = {'narrMass': 'red', '1Dfit': 'blue', '2Dfit': 'green'}
method_markers = {'narrMass': 'o', '1Dfit': 's', '2Dfit': '^'}
method_labels = {'narrMass': 'Narrow Mass Window', '1Dfit': '1D Fit (fixed mass2)', '2Dfit': '2D Template Fit'}

# Centrality过滤函数
def filter_centrality(df, max_cent=60):
    """过滤掉centrality > max_cent的数据"""
    return df[df['centrality'] <= max_cent]

# Ratio计算函数
def calculate_ratio(scan_value, ref_value, scan_err, ref_err):
    """计算ratio和其误差"""
    if ref_value == 0:
        return np.nan, np.nan
    ratio = scan_value / ref_value
    # 误差传播: σ(A/B) = |A/B| * sqrt((σA/A)² + (σB/B)²)
    ratio_err = abs(ratio) * np.sqrt((scan_err/abs(scan_value))**2 + (ref_err/abs(ref_value))**2) if scan_value != 0 else np.nan
    return ratio, ratio_err

# 获取特定range下的Narrow Mass数据作为参考
def get_narrmass_reference(range_val, pair_type, centrality):
    # 获取特定range下的Narrow Mass参考数据
    if not has_narrMass or range_val not in data_dict['narrMass']:
        return None

    df = filter_centrality(data_dict['narrMass'][range_val])
    subset = df[(df['pair_type'] == pair_type) & (df['centrality'] == centrality)]

    if len(subset) > 0:
        return subset
    return None

# 如果有2D Fit参考数据，添加2D Fit与不同range下的Narrow Mass比较图
if reference_2Dfit_data is not None and 'narrMass' in data_dict and data_dict['narrMass']:
    print("\n创建2D Fit与不同range的Narrow Mass比较图...")

    # 创建Delta SS ratio图
    fig_2d_delta, axes_2d_delta = plt.subplots(2, 2, figsize=(16, 12))
    fig_2d_delta.suptitle('2D Fit Delta SS Ratio (vs Narrow Mass) with Range on X-axis', fontsize=16)

    # 创建Rawgamma SS ratio图
    fig_2d_gamma, axes_2d_gamma = plt.subplots(2, 2, figsize=(16, 12))
    fig_2d_gamma.suptitle('2D Fit Rawgamma SS Ratio (vs Narrow Mass) with Range on X-axis', fontsize=16)

    # 获取所有narrMass的range值
    narr_range_values = sorted(data_dict['narrMass'].keys())

    # 筛选可用的中心度值
    available_centralities = set()
    for range_val in narr_range_values:
        df = filter_centrality(data_dict['narrMass'][range_val])
        for pair in pair_types:
            subset = df[df['pair_type'] == pair]
            available_centralities.update(subset['centrality'].unique())

    # 选择一些代表性的中心度值（最多6个，避免图表过于拥挤）
    available_centralities = sorted(list(available_centralities))
    if len(available_centralities) > 6:
        selected_centralities = np.linspace(min(available_centralities),
                                           max(available_centralities),
                                           6, endpoint=True).astype(int)
    else:
        selected_centralities = available_centralities

    # 为不同中心度定义不同的颜色和标记
    centrality_colors = plt.cm.tab10(np.linspace(0, 1, len(selected_centralities)))
    centrality_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    for i, pair_type in enumerate(pair_types):
        ax_delta = axes_2d_delta[i//2, i%2]
        ax_gamma = axes_2d_gamma[i//2, i%2]

        # 对于每个选定的中心度
        for idx, centrality in enumerate(selected_centralities):
            # 准备数据
            x_vals = []  # range值
            y_delta_vals = []  # Delta SS比值
            y_delta_errs = []  # Delta SS比值误差
            y_gamma_vals = []  # Rawgamma SS比值
            y_gamma_errs = []  # Rawgamma SS比值误差

            # 过滤2D Fit数据
            ref_df = filter_centrality(reference_2Dfit_data)
            ref_subset = ref_df[(ref_df['pair_type'] == pair_type) &
                              (ref_df['centrality'] == centrality)]

            if len(ref_subset) == 0:
                continue

            ref_val_delta = ref_subset['delta_ss'].iloc[0]
            ref_err_delta = ref_subset['delta_ss_err'].iloc[0]
            ref_val_gamma = ref_subset['rawgamma_ss'].iloc[0]
            ref_err_gamma = ref_subset['rawgamma_ss_err'].iloc[0]

            # 遍历所有range值
            for range_val in narr_range_values:
                narr_df = filter_centrality(data_dict['narrMass'][range_val])
                narr_subset = narr_df[(narr_df['pair_type'] == pair_type) &
                                    (narr_df['centrality'] == centrality)]

                if len(narr_subset) == 0:
                    continue

                narr_val_delta = narr_subset['delta_ss'].iloc[0]
                narr_err_delta = narr_subset['delta_ss_err'].iloc[0]
                narr_val_gamma = narr_subset['rawgamma_ss'].iloc[0]
                narr_err_gamma = narr_subset['rawgamma_ss_err'].iloc[0]

                # 计算比值
                ratio_delta, ratio_err_delta = calculate_ratio(ref_val_delta, narr_val_delta,
                                                            ref_err_delta, narr_err_delta)
                ratio_gamma, ratio_err_gamma = calculate_ratio(ref_val_gamma, narr_val_gamma,
                                                            ref_err_gamma, narr_err_gamma)

                if not np.isnan(ratio_delta):
                    x_vals.append(range_val)
                    y_delta_vals.append(ratio_delta)
                    y_delta_errs.append(ratio_err_delta)

                if not np.isnan(ratio_gamma):
                    y_gamma_vals.append(ratio_gamma)
                    y_gamma_errs.append(ratio_err_gamma)

            # 如果有数据，绘制曲线
            if x_vals:
                color = centrality_colors[idx % len(centrality_colors)]
                marker = centrality_markers[idx % len(centrality_markers)]

                # 确保数据按range值排序
                sorted_indices = np.argsort(x_vals)
                x_sorted = [x_vals[i] for i in sorted_indices]

                # Delta SS比值图
                y_delta_sorted = [y_delta_vals[i] for i in sorted_indices]
                y_delta_err_sorted = [y_delta_errs[i] for i in sorted_indices]
                ax_delta.errorbar(x_sorted, y_delta_sorted, yerr=y_delta_err_sorted,
                              label=f'Centrality {centrality}%',
                              color=color, marker=marker, markersize=5,
                              linewidth=1.5, capsize=3)

                # Rawgamma SS比值图
                if len(y_gamma_vals) == len(x_vals):  # 确保数据长度一致
                    y_gamma_sorted = [y_gamma_vals[i] for i in sorted_indices]
                    y_gamma_err_sorted = [y_gamma_errs[i] for i in sorted_indices]
                    ax_gamma.errorbar(x_sorted, y_gamma_sorted, yerr=y_gamma_err_sorted,
                                  label=f'Centrality {centrality}%',
                                  color=color, marker=marker, markersize=5,
                                  linewidth=1.5, capsize=3)

        # 添加比值=1的参考线
        ax_delta.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Ratio = 1')
        ax_gamma.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Ratio = 1')

        # 设置标题和标签
        ax_delta.set_title(f'{pair_type}', fontsize=12, fontweight='bold')
        ax_delta.set_xlabel('Range Mass Value', fontsize=10)
        ax_delta.set_ylabel('2D Fit / Narrow Mass Ratio', fontsize=10)
        ax_delta.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax_delta.grid(True, alpha=0.3)

        ax_gamma.set_title(f'{pair_type}', fontsize=12, fontweight='bold')
        ax_gamma.set_xlabel('Range Mass Value', fontsize=10)
        ax_gamma.set_ylabel('2D Fit / Narrow Mass Ratio', fontsize=10)
        ax_gamma.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax_gamma.grid(True, alpha=0.3)

    # 保存图表
    plt.figure(fig_2d_delta.number)
    plt.tight_layout()
    plt.savefig('2Dfit_vs_narrMass_delta_ss_ratio_by_range.pdf', bbox_inches='tight')
    plt.close(fig_2d_delta)

    plt.figure(fig_2d_gamma.number)
    plt.tight_layout()
    plt.savefig('2Dfit_vs_narrMass_rawgamma_ss_ratio_by_range.pdf', bbox_inches='tight')
    plt.close(fig_2d_gamma)

# 为每个方法创建range_mass扫描图
for method in data_dict.keys():
    if not data_dict[method]:
        continue

    # 获取所有range值并排序
    range_values = sorted(data_dict[method].keys())
    print(f"\n{method} 方法的range值: {range_values}")

    # 创建Delta SS ratio vs centrality图
    fig1, axes1 = plt.subplots(2, 2, figsize=(16, 12))
    fig1.suptitle(f'Delta SS Ratio (vs Narrow Mass) - {method_labels.get(method, method)}', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes1[i//2, i%2]

        # 为每个range值创建一条线
        for range_val in range_values:
            df = filter_centrality(data_dict[method][range_val])
            subset = df[df['pair_type'] == pair_type].sort_values('centrality')

            if len(subset) > 0:
                x_vals, y_vals, y_errs = [], [], []

                for _, row in subset.iterrows():
                    centrality = row['centrality']
                    scan_val = row['delta_ss']
                    scan_err = row['delta_ss_err']

                    # 获取相同range下的Narrow Mass参考值
                    ref_subset = get_narrmass_reference(range_val, pair_type, centrality)

                    if method == 'narrMass':
                        # 如果当前方法就是Narrow Mass，ratio为1
                        ratio, ratio_err = 1.0, 0.0
                        x_vals.append(centrality)
                        y_vals.append(ratio)
                        y_errs.append(ratio_err)
                    elif ref_subset is not None:
                        # 其他方法使用相同range下的Narrow Mass作为参考
                        ref_val = ref_subset['delta_ss'].iloc[0]
                        ref_err = ref_subset['delta_ss_err'].iloc[0]
                        ratio, ratio_err = calculate_ratio(scan_val, ref_val, scan_err, ref_err)

                        if not np.isnan(ratio):
                            x_vals.append(centrality)
                            y_vals.append(ratio)
                            y_errs.append(ratio_err)

                if x_vals:
                    ax.errorbar(x_vals, y_vals, yerr=y_errs,
                               label=f'Range {range_val:.6f}',
                               marker='o', markersize=4, linewidth=1.5, capsize=3)

        # 添加ratio=1的参考线
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label="Narrow Mass reference")

        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Delta SS Ratio')
        ax.set_title(f'{pair_type}')
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig(f'delta_ss_ratio_cent_scan_{method}_vs_narrMass.pdf', bbox_inches='tight')
    plt.close()

    # 创建Rawgamma SS ratio vs centrality图
    fig2, axes2 = plt.subplots(2, 2, figsize=(16, 12))
    fig2.suptitle(f'Rawgamma SS Ratio (vs Narrow Mass) - {method_labels.get(method, method)}', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes2[i//2, i%2]

        # 为每个range值创建一条线
        for range_val in range_values:
            df = filter_centrality(data_dict[method][range_val])
            subset = df[df['pair_type'] == pair_type].sort_values('centrality')

            if len(subset) > 0:
                x_vals, y_vals, y_errs = [], [], []

                for _, row in subset.iterrows():
                    centrality = row['centrality']
                    scan_val = row['rawgamma_ss']
                    scan_err = row['rawgamma_ss_err']

                    # 获取相同range下的Narrow Mass参考值
                    ref_subset = get_narrmass_reference(range_val, pair_type, centrality)

                    if method == 'narrMass':
                        # 如果当前方法就是Narrow Mass，ratio为1
                        ratio, ratio_err = 1.0, 0.0
                        x_vals.append(centrality)
                        y_vals.append(ratio)
                        y_errs.append(ratio_err)
                    elif ref_subset is not None:
                        # 其他方法使用相同range下的Narrow Mass作为参考
                        ref_val = ref_subset['rawgamma_ss'].iloc[0]
                        ref_err = ref_subset['rawgamma_ss_err'].iloc[0]
                        ratio, ratio_err = calculate_ratio(scan_val, ref_val, scan_err, ref_err)

                        if not np.isnan(ratio):
                            x_vals.append(centrality)
                            y_vals.append(ratio)
                            y_errs.append(ratio_err)

                if x_vals:
                    ax.errorbar(x_vals, y_vals, yerr=y_errs,
                               label=f'Range {range_val:.6f}',
                               marker='o', markersize=4, linewidth=1.5, capsize=3)

        # 添加ratio=1的参考线
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label="Narrow Mass reference")

        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Rawgamma SS Ratio')
        ax.set_title(f'{pair_type}')
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig(f'rawgamma_ss_ratio_cent_scan_{method}_vs_narrMass.pdf', bbox_inches='tight')
    plt.close()

# 如果有多个方法，创建总的对比图
if len(data_dict) > 1:
    print(f"\n创建总的方法对比图...")

    # 创建Delta SS ratio总对比图
    fig_total1, axes_total1 = plt.subplots(2, 2, figsize=(20, 12))
    fig_total1.suptitle('Delta SS Ratio Comparison: Methods vs Reference', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes_total1[i//2, i%2]

        for method in data_dict.keys():
            # 获取该方法所有range值
            method_ranges = sorted(data_dict[method].keys())

            for range_val in method_ranges:
                df = filter_centrality(data_dict[method][range_val])
                subset = df[df['pair_type'] == pair_type].sort_values('centrality')

                if len(subset) > 0:
                    x_vals, y_vals, y_errs = [], [], []

                    for _, row in subset.iterrows():
                        centrality = row['centrality']
                        scan_val = row['delta_ss']
                        scan_err = row['delta_ss_err']

                        # 获取相同range下的Narrow Mass参考值
                        ref_subset = get_narrmass_reference(range_val, pair_type, centrality)

                        if method == 'narrMass':
                            # 如果当前方法就是Narrow Mass，ratio为1
                            ratio, ratio_err = 1.0, 0.0
                            x_vals.append(centrality)
                            y_vals.append(ratio)
                            y_errs.append(ratio_err)
                        elif ref_subset is not None:
                            # 其他方法使用相同range下的Narrow Mass作为参考
                            ref_val = ref_subset['delta_ss'].iloc[0]
                            ref_err = ref_subset['delta_ss_err'].iloc[0]
                            ratio, ratio_err = calculate_ratio(scan_val, ref_val, scan_err, ref_err)

                            if not np.isnan(ratio):
                                x_vals.append(centrality)
                                y_vals.append(ratio)
                                y_errs.append(ratio_err)

                    if x_vals:
                        # 使用不同的颜色和线型区分方法和range
                        base_color = method_colors.get(method, 'black')
                        marker = method_markers.get(method, 'o')

                        # 为不同range使用不同的透明度
                        alpha = 0.3 + 0.7 * (method_ranges.index(range_val) / len(method_ranges))

                        ax.errorbar(x_vals, y_vals, yerr=y_errs,
                                   label=f'{method} Range {range_val:.6f}',
                                   color=base_color, marker=marker, alpha=alpha,
                                   markersize=3, linewidth=1.2, capsize=2)

        # 添加ratio=1的参考线
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.8, linewidth=2, label='Narrow Mass reference')

        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Delta SS Ratio')
        ax.set_title(f'{pair_type}')
        ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig('delta_ss_ratio_total_comparison_vs_narrMass.pdf', bbox_inches='tight')
    plt.close()

    # 创建Rawgamma SS ratio总对比图
    fig_total2, axes_total2 = plt.subplots(2, 2, figsize=(20, 12))
    fig_total2.suptitle('Rawgamma SS Ratio Comparison: Methods vs Reference', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes_total2[i//2, i%2]

        for method in data_dict.keys():
            # 获取该方法所有range值
            method_ranges = sorted(data_dict[method].keys())

            for range_val in method_ranges:
                df = filter_centrality(data_dict[method][range_val])
                subset = df[df['pair_type'] == pair_type].sort_values('centrality')

                if len(subset) > 0:
                    x_vals, y_vals, y_errs = [], [], []

                    for _, row in subset.iterrows():
                        centrality = row['centrality']
                        scan_val = row['rawgamma_ss']
                        scan_err = row['rawgamma_ss_err']

                        # 获取相同range下的Narrow Mass参考值
                        ref_subset = get_narrmass_reference(range_val, pair_type, centrality)

                        if method == 'narrMass':
                            # 如果当前方法就是Narrow Mass，ratio为1
                            ratio, ratio_err = 1.0, 0.0
                            x_vals.append(centrality)
                            y_vals.append(ratio)
                            y_errs.append(ratio_err)
                        elif ref_subset is not None:
                            # 其他方法使用相同range下的Narrow Mass作为参考
                            ref_val = ref_subset['rawgamma_ss'].iloc[0]
                            ref_err = ref_subset['rawgamma_ss_err'].iloc[0]
                            ratio, ratio_err = calculate_ratio(scan_val, ref_val, scan_err, ref_err)

                            if not np.isnan(ratio):
                                x_vals.append(centrality)
                                y_vals.append(ratio)
                                y_errs.append(ratio_err)

                    if x_vals:
                        # 使用不同的颜色和线型区分方法和range
                        base_color = method_colors.get(method, 'black')
                        marker = method_markers.get(method, 'o')

                        # 为不同range使用不同的透明度
                        alpha = 0.3 + 0.7 * (method_ranges.index(range_val) / len(method_ranges))

                        ax.errorbar(x_vals, y_vals, yerr=y_errs,
                                   label=f'{method} Range {range_val:.6f}',
                                   color=base_color, marker=marker, alpha=alpha,
                                   markersize=3, linewidth=1.2, capsize=2)

        # 添加ratio=1的参考线
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.8, linewidth=2, label='Reference')

        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Rawgamma SS Ratio')
        ax.set_title(f'{pair_type}')
        ax.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 60)

    plt.tight_layout()
    plt.savefig('rawgamma_ss_ratio_total_comparison_vs_narrMass.pdf', bbox_inches='tight')
    plt.close()

# 如果有多个方法，创建方法对比图（在固定range下）
if len(data_dict) > 1:
    # 找到所有方法都有的range值
    common_ranges = set(data_dict[list(data_dict.keys())[0]].keys())
    for method in list(data_dict.keys())[1:]:
        common_ranges &= set(data_dict[method].keys())

    if common_ranges:
        # 选择一个中等的range值进行对比
        comparison_range = sorted(list(common_ranges))[len(common_ranges)//2]
        print(f"\n使用 range={comparison_range:.6f} 进行方法对比")

        # 创建方法对比图
        fig3, axes3 = plt.subplots(2, 2, figsize=(16, 12))
        fig3.suptitle(f'Method Comparison at Range Mass = {comparison_range:.6f}', fontsize=16)

        for i, pair_type in enumerate(pair_types):
            ax = axes3[i//2, i%2]

            for method in data_dict.keys():
                if comparison_range in data_dict[method]:
                    df = data_dict[method][comparison_range]
                    data = df[df['pair_type'] == pair_type].sort_values('centrality')
                    if len(data) > 0:
                        ax.errorbar(data['centrality'], data['delta_ss'], yerr=data['delta_ss_err'],
                                   label=method_labels.get(method, method),
                                   color=method_colors.get(method, 'black'),
                                   marker=method_markers.get(method, 'o'),
                                   markersize=6, linewidth=2, capsize=4)

            ax.set_xlabel('Centrality (%)')
            ax.set_ylabel('Delta SS')
            ax.set_title(f'{pair_type}')
            ax.legend()
            ax.set_xlim(0, 70)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'method_comparison_vs_narrMass_range{comparison_range:.6f}.pdf')
        plt.close()

# 打印统计信息
print("\n=== 扫描结果统计 ===")
for method in data_dict.keys():
    print(f"\n{method} 方法:")
    print(f"  扫描的range数量: {len(data_dict[method])}")
    print(f"  Range范围: {min(data_dict[method].keys()):.6f} - {max(data_dict[method].keys()):.6f}")

    # 统计每个range的数据点数
    for range_val in sorted(data_dict[method].keys()):
        df = data_dict[method][range_val]
        print(f"    Range {range_val:.6f}: {len(df)} 数据点")

print(f"\n生成的图片文件:")
for method in data_dict.keys():
    print(f"- delta_ss_ratio_cent_scan_{method}_vs_narrMass.pdf")
    print(f"- rawgamma_ss_ratio_cent_scan_{method}_vs_narrMass.pdf")

if len(data_dict) > 1:
    print("- delta_ss_ratio_total_comparison_vs_narrMass.pdf")
    print("- rawgamma_ss_ratio_total_comparison_vs_narrMass.pdf")

if len(data_dict) > 1 and 'common_ranges' in locals() and common_ranges:
    print(f"- method_comparison_range{comparison_range:.6f}.pdf")

if reference_2Dfit_data is not None and 'narrMass' in data_dict and data_dict['narrMass']:
    print("- 2Dfit_vs_narrMass_delta_ss_ratio_by_range.pdf")
    print("- 2Dfit_vs_narrMass_rawgamma_ss_ratio_by_range.pdf")

if not has_narrMass:
    print("\n注意: 未找到Narrow Mass方法的数据，无法计算ratio")
