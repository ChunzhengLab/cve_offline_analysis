import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# 查找所有拟合结果文件
result_files = glob.glob('LHC18q/Lambda/fit_obvs_*.csv')
print(f"找到的结果文件: {result_files}")

# 读取所有方法的数据
data_dict = {}
for file in result_files:
    base = os.path.basename(file)  # e.g. "fit_obvs_1Dfit_default.csv"
    method = base.replace("fit_obvs_", "").replace(".csv", "")  # "1Dfit_default"
    method = method.split("_")[0]  # 只留 "1Dfit"
    print(f"读取方法: {method}")
    data_dict[method] = pd.read_csv(file)

# 如果没有找到分方法的文件，尝试读取原始文件
if not data_dict and os.path.exists('fit_obvs.csv'):
    print("读取原始fit_obvs.csv文件")
    df_orig = pd.read_csv('fit_obvs.csv')
    if 'fit_mode' in df_orig.columns:
        for mode in df_orig['fit_mode'].unique():
            data_dict[mode] = df_orig[df_orig['fit_mode'] == mode]
    else:
        data_dict['2Dfit'] = df_orig

if not data_dict:
    raise FileNotFoundError("未找到任何拟合结果文件")

print(f"可用的拟合方法: {list(data_dict.keys())}")

# 定义粒子对类型、颜色和标记
pair_types = ['LambdaLambda', 'LambdaLambdaBar', 'LambdaBarLambda', 'LambdaBarLambdaBar']
method_colors = {'2Dfit': 'red', '1Dfit': 'blue', 'narrMass': 'green'}
method_markers = {'2Dfit': 'o', '1Dfit': 's', 'narrMass': '^'}
method_labels = {'2Dfit': '2D Template Fit', '1Dfit': '1D Fit (fixed mass2)', 'narrMass': 'Narrow Mass Window'}

# 创建Delta SS对比图
fig1, axes1 = plt.subplots(2, 2, figsize=(16, 12))
fig1.suptitle('Delta SS Parameter Comparison - Different Methods', fontsize=16)

for i, pair_type in enumerate(pair_types):
    ax = axes1[i//2, i%2]
    
    for method, df in data_dict.items():
        if method not in method_colors:
            continue
        data = df[df['pair_type'] == pair_type].sort_values('centrality')
        data = data[data['centrality'] <= 60]  # 只保留centrality<=60
        if len(data) > 0:
            ax.errorbar(data['centrality'], data['delta_ss'], yerr=data['delta_ss_err'], 
                       label=method_labels.get(method, method), 
                       color=method_colors[method], marker=method_markers[method], 
                       markersize=6, linewidth=2, capsize=4)
    
    ax.set_xlabel('Centrality (%)')
    ax.set_ylabel('Delta SS')
    ax.set_title(f'{pair_type}')
    ax.legend()
    ax.set_xlim(0, 60) 
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delta_ss_method_comparison.pdf')

# 创建Rawgamma SS对比图
fig2, axes2 = plt.subplots(2, 2, figsize=(16, 12))
fig2.suptitle('Rawgamma SS Parameter Comparison - Different Methods', fontsize=16)

for i, pair_type in enumerate(pair_types):
    ax = axes2[i//2, i%2]
    
    for method, df in data_dict.items():
        if method not in method_colors:
            continue
        data = df[df['pair_type'] == pair_type].sort_values('centrality')
        data = data[data['centrality'] <= 60]  # 只保留centrality<=60
        if len(data) > 0:
            ax.errorbar(data['centrality'], data['rawgamma_ss'], yerr=data['rawgamma_ss_err'], 
                       label=method_labels.get(method, method), 
                       color=method_colors[method], marker=method_markers[method], 
                       markersize=6, linewidth=2, capsize=4)
    
    ax.set_xlabel('Centrality (%)')
    ax.set_ylabel('Rawgamma SS')
    ax.set_title(f'{pair_type}')
    ax.legend()
    ax.set_xlim(0, 60)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('rawgamma_ss_method_comparison.pdf')

# 打印方法对比统计
print("\n=== 方法对比统计 ===")
for method, df in data_dict.items():
    print(f"\n{method} 方法:")
    print(f"  数据点数: {len(df)}")
    if 'delta_ss' in df.columns:
        print(f"  Delta SS 平均值: {df['delta_ss'].mean():.6f} ± {df['delta_ss'].std():.6f}")
    if 'rawgamma_ss' in df.columns:
        print(f"  Rawgamma SS 平均值: {df['rawgamma_ss'].mean():.6f} ± {df['rawgamma_ss'].std():.6f}")

print(f"\n生成的图片文件:")
print("- delta_ss_method_comparison.png")
print("- rawgamma_ss_method_comparison.png")

def plot_method_ratio(data_dict, pair_types):
    if not all(k in data_dict for k in ['2Dfit', '1Dfit', 'narrMass']):
        print("缺少2Dfit、1Dfit或narrMass数据，无法画比值图")
        return

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Method / Narrow Mass Ratio vs Centrality', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes[i//2, i%2]
        # merge三种方法，保证中心度一一对应
        df_2d = data_dict['2Dfit'][data_dict['2Dfit']['pair_type'] == pair_type]
        df_1d = data_dict['1Dfit'][data_dict['1Dfit']['pair_type'] == pair_type]
        df_nm = data_dict['narrMass'][data_dict['narrMass']['pair_type'] == pair_type]
        df_2d = df_2d[df_2d['centrality'] <= 60]
        df_1d = df_1d[df_1d['centrality'] <= 60]
        df_nm = df_nm[df_nm['centrality'] <= 60]

        # 只保留三者都存在的中心度
        merged_2d = pd.merge(df_2d, df_nm, on='centrality', suffixes=('_2d', '_nm'))
        merged_1d = pd.merge(df_1d, df_nm, on='centrality', suffixes=('_1d', '_nm'))

        # 2Dfit/narrMass
        if len(merged_2d) > 0:
            ratio_2d = merged_2d['delta_ss_2d'] / merged_2d['delta_ss_nm']
            # 误差传播
            err_2d = np.abs(ratio_2d) * np.sqrt(
                (merged_2d['delta_ss_err_2d'] / merged_2d['delta_ss_2d'])**2 +
                (merged_2d['delta_ss_err_nm'] / merged_2d['delta_ss_nm'])**2
            )
            ax.errorbar(merged_2d['centrality'], ratio_2d, yerr=err_2d,
                        label='2D Fit / narrMass', color='red', marker='o', linestyle='-', linewidth=2, capsize=4)

        # 1Dfit/narrMass
        if len(merged_1d) > 0:
            ratio_1d = merged_1d['delta_ss_1d'] / merged_1d['delta_ss_nm']
            err_1d = np.abs(ratio_1d) * np.sqrt(
                (merged_1d['delta_ss_err_1d'] / merged_1d['delta_ss_1d'])**2 +
                (merged_1d['delta_ss_err_nm'] / merged_1d['delta_ss_nm'])**2
            )
            ax.errorbar(merged_1d['centrality'], ratio_1d, yerr=err_1d,
                        label='1D Fit / narrMass', color='blue', marker='s', linestyle='--', linewidth=2, capsize=4)

        ax.axhline(1, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Delta SS Ratio')
        ax.set_title(pair_type)
        ax.legend()
        ax.set_xlim(0, 60)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('delta_ss_method_ratio_vs_narrMass.pdf')
    print('已生成比值图: delta_ss_method_ratio_vs_narrMass.pdf')

plot_method_ratio(data_dict, pair_types)

def plot_rawgamma_method_ratio(data_dict, pair_types):
    if not all(k in data_dict for k in ['2Dfit', '1Dfit', 'narrMass']):
        print("缺少2Dfit、1Dfit或narrMass数据，无法画rawgamma比值图")
        return

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Rawgamma SS: Method / Narrow Mass Ratio vs Centrality', fontsize=16)

    for i, pair_type in enumerate(pair_types):
        ax = axes[i//2, i%2]
        df_2d = data_dict['2Dfit'][data_dict['2Dfit']['pair_type'] == pair_type]
        df_1d = data_dict['1Dfit'][data_dict['1Dfit']['pair_type'] == pair_type]
        df_nm = data_dict['narrMass'][data_dict['narrMass']['pair_type'] == pair_type]
        df_2d = df_2d[df_2d['centrality'] <= 60]
        df_1d = df_1d[df_1d['centrality'] <= 60]
        df_nm = df_nm[df_nm['centrality'] <= 60]

        merged_2d = pd.merge(df_2d, df_nm, on='centrality', suffixes=('_2d', '_nm'))
        merged_1d = pd.merge(df_1d, df_nm, on='centrality', suffixes=('_1d', '_nm'))

        # 2Dfit/narrMass
        if len(merged_2d) > 0:
            ratio_2d = merged_2d['rawgamma_ss_2d'] / merged_2d['rawgamma_ss_nm']
            err_2d = np.abs(ratio_2d) * np.sqrt(
                (merged_2d['rawgamma_ss_err_2d'] / merged_2d['rawgamma_ss_2d'])**2 +
                (merged_2d['rawgamma_ss_err_nm'] / merged_2d['rawgamma_ss_nm'])**2
            )
            ax.errorbar(merged_2d['centrality'], ratio_2d, yerr=err_2d,
                        label='2D Fit / narrMass', color='red', marker='o', linestyle='-', linewidth=2, capsize=4)

        # 1Dfit/narrMass
        if len(merged_1d) > 0:
            ratio_1d = merged_1d['rawgamma_ss_1d'] / merged_1d['rawgamma_ss_nm']
            err_1d = np.abs(ratio_1d) * np.sqrt(
                (merged_1d['rawgamma_ss_err_1d'] / merged_1d['rawgamma_ss_1d'])**2 +
                (merged_1d['rawgamma_ss_err_nm'] / merged_1d['rawgamma_ss_nm'])**2
            )
            ax.errorbar(merged_1d['centrality'], ratio_1d, yerr=err_1d,
                        label='1D Fit / narrMass', color='blue', marker='s', linestyle='--', linewidth=2, capsize=4)

        ax.axhline(1, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Centrality (%)')
        ax.set_ylabel('Rawgamma SS Ratio')
        ax.set_title(pair_type)
        ax.legend()
        ax.set_xlim(0, 60)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('rawgamma_ss_method_ratio_vs_narrMass.pdf')
    print('已生成比值图: rawgamma_ss_method_ratio_vs_narrMass.pdf')

plot_rawgamma_method_ratio(data_dict, pair_types)