import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import re

# 查找所有扫描结果文件
result_files = glob.glob('outputs_scan/fit_obvs_default_*.csv')
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

# 读取2D Fit参考数据
reference_2Dfit_file = 'fit_obvs_2Dfit_default.csv'
reference_2Dfit_data = None

if os.path.exists(reference_2Dfit_file):
    reference_2Dfit_data = pd.read_csv(reference_2Dfit_file)
    print(f"读取2D Fit参考数据: {reference_2Dfit_file}")
else:
    print(f"警告: 未找到2D Fit参考文件 {reference_2Dfit_file}")

# 读取1D Fit参考数据
reference_1Dfit_file = 'fit_obvs_1Dfit_default.csv'
reference_1Dfit_data = None

if os.path.exists(reference_1Dfit_file):
    reference_1Dfit_data = pd.read_csv(reference_1Dfit_file)
    print(f"读取1D Fit参考数据: {reference_1Dfit_file}")
else:
    print(f"警告: 未找到1D Fit参考文件 {reference_1Dfit_file}")

# 检查是否有narrMass方法的数据
if 'narrMass' not in data_dict:
    print("警告: 没有找到Narrow Mass方法的数据，无法计算比值")
    has_narrMass = False
else:
    has_narrMass = True
    print(f"找到Narrow Mass方法的数据，可用的range值: {sorted(data_dict['narrMass'].keys())}")

# 定义粒子对类型、颜色和标记
pair_types = ['LambdaLambda', 'LambdaLambdaBar', 'LambdaBarLambda', 'LambdaBarLambdaBar']
method_colors = {'2Dfit': 'green', '1Dfit': 'blue', 'narrMass': 'red'}
method_markers = {'2Dfit': '^', '1Dfit': 's', 'narrMass': 'o'}
method_labels = {'2Dfit': '2D Template Fit', '1Dfit': '1D Fit (fixed mass2)', 'narrMass': 'Narrow Mass Window'}

# Centrality过滤函数
def filter_centrality(df, max_cent=60):
    """过滤掉centrality > max_cent的数据"""
    return df[df['centrality'] <= max_cent]

# Ratio计算函数
def calculate_ratio(num_value, denom_value, num_err, denom_err):
    """计算ratio和其误差"""
    if denom_value == 0:
        return np.nan, np.nan
    ratio = num_value / denom_value
    # 误差传播: σ(A/B) = |A/B| * sqrt((σA/A)² + (σB/B)²)
    ratio_err = abs(ratio) * np.sqrt((num_err/abs(num_value))**2 + (denom_err/abs(denom_value))**2) if num_value != 0 else np.nan
    return ratio, ratio_err

# 获取特定range下的Narrow Mass数据作为参考
def get_narrmass_reference(range_val, pair_type, centrality):
    """获取特定range下的Narrow Mass参考数据"""
    if not has_narrMass or range_val not in data_dict['narrMass']:
        return None

    df = filter_centrality(data_dict['narrMass'][range_val])
    subset = df[(df['pair_type'] == pair_type) & (df['centrality'] == centrality)]

    if len(subset) > 0:
        return subset
    return None

# 主要分析函数：创建以Range为X轴的比值图
def create_range_ratio_plots(method_name, param_name, param_label):
    """创建以Range为X轴的比值图，分子分母都是同一range下的2Dfit/1Dfit和narrMass"""
    if not has_narrMass or method_name not in data_dict:
        print(f"警告: 无法创建{method_name}与Narrow Mass的比值图，缺少必要数据")
        return []

    # 获取所有range值（两者都要有）
    narr_range_values = sorted(set(data_dict['narrMass'].keys()) & set(data_dict[method_name].keys()))

    # 筛选可用的中心度值
    available_centralities = set()
    for range_val in narr_range_values:
        df = filter_centrality(data_dict['narrMass'][range_val])
        for pair in pair_types:
            subset = df[df['pair_type'] == pair]
            available_centralities.update(subset['centrality'].unique())
    available_centralities = sorted(list(available_centralities))
    if len(available_centralities) > 6:
        indices = np.linspace(0, len(available_centralities)-1, 6, dtype=int)
        selected_centralities = [available_centralities[i] for i in indices]
    else:
        selected_centralities = available_centralities

    centrality_colors = plt.cm.tab10(np.linspace(0, 1, len(selected_centralities)))
    centrality_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'{method_name} {param_label} Ratio (vs Narrow Mass) with Range on X-axis', fontsize=16)
    all_data = []

    for i, pair_type in enumerate(pair_types):
        ax = axes[i//2, i%2]
        for idx, centrality in enumerate(selected_centralities):
            x_vals, y_vals, y_errs = [], [], []
            for range_val in narr_range_values:
                narr_df = filter_centrality(data_dict['narrMass'][range_val])
                method_df = filter_centrality(data_dict[method_name][range_val])
                narr_subset = narr_df[(narr_df['pair_type'] == pair_type) & (narr_df['centrality'] == centrality)]
                method_subset = method_df[(method_df['pair_type'] == pair_type) & (method_df['centrality'] == centrality)]
                if len(narr_subset) == 0 or len(method_subset) == 0:
                    continue
                narr_val = narr_subset[param_name].iloc[0]
                narr_err = narr_subset[f'{param_name}_err'].iloc[0]
                method_val = method_subset[param_name].iloc[0]
                method_err = method_subset[f'{param_name}_err'].iloc[0]
                ratio, ratio_err = calculate_ratio(method_val, narr_val, method_err, narr_err)
                if not np.isnan(ratio):
                    x_vals.append(range_val)
                    y_vals.append(ratio)
                    y_errs.append(ratio_err)
                    all_data.append({
                        'method': method_name,
                        'pair_type': pair_type,
                        'centrality': centrality,
                        'range': range_val,
                        'ratio': ratio,
                        'ratio_err': ratio_err
                    })
            if x_vals:
                color = centrality_colors[idx % len(centrality_colors)]
                marker = centrality_markers[idx % len(centrality_markers)]
                sorted_indices = np.argsort(x_vals)
                x_sorted = [x_vals[i] for i in sorted_indices]
                y_sorted = [y_vals[i] for i in sorted_indices]
                y_err_sorted = [y_errs[i] for i in sorted_indices]
                ax.errorbar(x_sorted, y_sorted, yerr=y_err_sorted,
                          label=f'Centrality {centrality}%',
                          color=color, marker=marker, markersize=5,
                          linewidth=1.5, capsize=3)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Ratio = 1')
        ax.set_title(f'{pair_type}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Range Mass Value', fontsize=10)
        ax.set_ylabel(f'{method_name} / Narrow Mass Ratio', fontsize=10)
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{method_name.lower()}_vs_narrMass_{param_name}_ratio_by_range.pdf', bbox_inches='tight')
    plt.close()
    return all_data

# 创建所有方法和所有粒子对类型的汇总图
def create_all_methods_all_pairs_plot(all_data_2D, all_data_1D, param_name, param_label):
    """创建2D Fit和1D Fit对Narrow Mass比值的比较图，自动用all_data里实际有的pair_type和中心度"""
    if not all_data_2D and not all_data_1D:
        print(f"警告: 无法创建汇总图，缺少必要数据")
        return
    import pandas as pd
    all_data = pd.DataFrame(all_data_2D + all_data_1D)
    # 自动用实际有的pair_type
    pair_types_in_data = sorted(all_data['pair_type'].unique())
    # 自动选择实际有的中心度
    centralities = sorted(all_data['centrality'].unique())
    if len(centralities) >= 3:
        actual_centralities = [centralities[0], centralities[len(centralities)//2], centralities[-1]]
    else:
        actual_centralities = centralities
    print('实际用于画图的pair_type:', pair_types_in_data)
    print('实际用于画图的中心度:', actual_centralities)
    fig, axes = plt.subplots(2, 2, figsize=(20, 14))
    fig.suptitle(f'All Methods & Pair Types: {param_label} Ratio vs Narrow Mass', fontsize=18, fontweight='bold')
    method_colors = {
        '2Dfit': '#1f77b4',
        '1Dfit': '#d62728'
    }
    cent_linestyles = {
        actual_centralities[0]: '-',
        actual_centralities[1]: '--',
        actual_centralities[2]: '-.'
    } if len(actual_centralities) == 3 else {}
    cent_markers = {
        actual_centralities[0]: 'o',
        actual_centralities[1]: 's',
        actual_centralities[2]: '^'
    } if len(actual_centralities) == 3 else {}
    for i, pair_type in enumerate(pair_types_in_data):
        ax = axes[i//2, i%2]
        for method_idx, method in enumerate(['2Dfit', '1Dfit']):
            for cent_idx, centrality in enumerate(actual_centralities):
                method_data = all_data[(all_data['method'] == method) &
                                      (all_data['pair_type'] == pair_type) &
                                      (np.isclose(all_data['centrality'], centrality, atol=1e-3))]
                print(f"方法: {method}, pair_type: {pair_type}, centrality: {centrality}, 数据点数: {len(method_data)}")
                if len(method_data) > 0:
                    x_vals = method_data['range'].values
                    y_vals = method_data['ratio'].values
                    y_errs = method_data['ratio_err'].values
                    sorted_indices = np.argsort(x_vals)
                    x_sorted = x_vals[sorted_indices]
                    y_sorted = y_vals[sorted_indices]
                    y_err_sorted = y_errs[sorted_indices]
                    current_color = method_colors[method]
                    current_linestyle = cent_linestyles.get(centrality, '-')
                    current_marker = cent_markers.get(centrality, 'o')
                    # 图例美化
                    method_label = method.replace('fit', ' Fit')
                    ax.errorbar(x_sorted, y_sorted, yerr=y_err_sorted,
                                label=f"{method_label} - Cent {centrality}%",
                                color=current_color, marker=current_marker, markersize=8,
                                linestyle=current_linestyle,
                                linewidth=2.5, capsize=4, elinewidth=1.5)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Range Mass Value', fontsize=12)
        ax.set_ylabel('Method / Narrow Mass Ratio', fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.grid(True, alpha=0.3)
        legend = ax.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left',
                           frameon=True, facecolor='white', edgecolor='gray')
        ax.minorticks_on()
        ax.tick_params(which='minor', length=4, width=1)
        ax.tick_params(which='major', length=7, width=1.5)
    plt.tight_layout()
    plt.savefig(f'all_methods_all_pairs_{param_name}_ratio.pdf',
                bbox_inches='tight', dpi=300)
    plt.savefig(f'all_methods_all_pairs_{param_name}_ratio.png',
                bbox_inches='tight', dpi=300)
    plt.close()

def main():
    """主函数"""
    if not has_narrMass:
        print("错误: 未找到Narrow Mass方法的数据，无法进行分析")
        return
    print("\n创建Delta SS比值图...")
    all_data_2D_delta = create_range_ratio_plots('2Dfit', 'delta_ss', 'Delta SS')
    all_data_1D_delta = create_range_ratio_plots('1Dfit', 'delta_ss', 'Delta SS')
    print("\n创建Rawgamma SS比值图...")
    all_data_2D_gamma = create_range_ratio_plots('2Dfit', 'rawgamma_ss', 'Rawgamma SS')
    all_data_1D_gamma = create_range_ratio_plots('1Dfit', 'rawgamma_ss', 'Rawgamma SS')
    print("\n创建所有方法、所有粒子对类型和所有指定中心度的汇总图...")
    create_all_methods_all_pairs_plot(all_data_2D_delta, all_data_1D_delta, 'delta_ss', 'Delta SS')
    create_all_methods_all_pairs_plot(all_data_2D_gamma, all_data_1D_gamma, 'rawgamma_ss', 'Rawgamma SS')

if __name__ == "__main__":
    main()
