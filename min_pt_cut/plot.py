import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# 读取所有 finalise_*.csv 文件
csv_files = glob.glob('./LHC18q/Proton/finalise_*.csv')

# Debug: 打印找到的文件
if not csv_files:
    raise FileNotFoundError("No CSV files matching './LHC18q/Proton/finalise_*.csv' found.")

# 读取并合并所有CSV文件
dfs = [pd.read_csv(f) for f in csv_files]

# 添加文件名列以便区分不同来源，只保留简短名称
simplified_names = {
    'finalise_default.csv': 'default',
    'finalise_LambdaProtonPtCut.csv': 'LambdaProtonPtCut',
    'finalise_LambdaPtCut.csv': 'LambdaPtCut',
    'finalise_ProtonPtCut.csv': 'ProtonPtCut'
}

dataframes = []
for f, df in zip(csv_files, dfs):
    basename = os.path.basename(f)
    df['source_file'] = simplified_names.get(basename, basename)
    dataframes.append(df)

data = pd.concat(dataframes, ignore_index=True)

# 过滤数据
data_filtered = data[
    (data['diff_type'] == 'Intg') &
    (data['diff_bin'] == 0.5) &
    (data['pair_type'].isin(['SS', 'OS'])) &
    (data['centrality'] <= 60)
]

# 设置输出目录
output_dir = './plots'
os.makedirs(output_dir, exist_ok=True)

# 获取不同来源文件名
different_sources = data_filtered['source_file'].unique()
markers = ['o', 's', 'D', '^', 'v', '<', '>']
colors = plt.cm.tab10.colors  # 使用tab10配色，最多10种颜色

# 1. 画总的图
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for i, source in enumerate(different_sources):
    source_data = data_filtered[data_filtered['source_file'] == source]
    ss_data = source_data[source_data['pair_type'] == 'SS'].sort_values('centrality')
    os_data = source_data[source_data['pair_type'] == 'OS'].sort_values('centrality')

    color = colors[i % len(colors)]

    axes[0].errorbar(ss_data['centrality'], ss_data['delta'], yerr=ss_data['delta_err'], fmt=markers[i % len(markers)]+'-', color=color, linewidth=2)
    axes[0].errorbar(os_data['centrality'], os_data['delta'], yerr=os_data['delta_err'], fmt=markers[i % len(markers)]+'--', color=color, linewidth=2, label=source)

    axes[1].errorbar(ss_data['centrality'], ss_data['gamma'], yerr=ss_data['gamma_err'], fmt=markers[i % len(markers)]+'-', color=color, linewidth=2)
    axes[1].errorbar(os_data['centrality'], os_data['gamma'], yerr=os_data['gamma_err'], fmt=markers[i % len(markers)]+'--', color=color, linewidth=2, label=source)

axes[0].set_xlabel('Centrality')
axes[0].set_ylabel('Delta')
axes[0].set_title('Delta vs Centrality')
axes[0].legend()
axes[0].grid(True)

axes[1].set_xlabel('Centrality')
axes[1].set_ylabel('Gamma')
axes[1].set_title('Gamma vs Centrality')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig(f'{output_dir}/all_sources.pdf')

# 保存当前轴范围，后面比照用
xlim0, ylim0 = axes[0].get_xlim(), axes[0].get_ylim()
xlim1, ylim1 = axes[1].get_xlim(), axes[1].get_ylim()
plt.close()

# 2. 分别和default对比
base_name = 'default'
comparison_files = [src for src in different_sources if src != base_name]

for comp_file in comparison_files:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, source in enumerate([base_name, comp_file]):
        source_data = data_filtered[data_filtered['source_file'] == source]
        ss_data = source_data[source_data['pair_type'] == 'SS'].sort_values('centrality')
        os_data = source_data[source_data['pair_type'] == 'OS'].sort_values('centrality')

        color = colors[idx % len(colors)]
        style = '-' if source == base_name else '--'

        # 只在第一条曲线加label
        axes[0].errorbar(ss_data['centrality'], ss_data['delta'], yerr=ss_data['delta_err'], fmt=style, color=color, linewidth=2, label=source)
        axes[0].errorbar(os_data['centrality'], os_data['delta'], yerr=os_data['delta_err'], fmt=style, color=color, linewidth=2)

        axes[1].errorbar(ss_data['centrality'], ss_data['gamma'], yerr=ss_data['gamma_err'], fmt=style, color=color, linewidth=2, label=source)
        axes[1].errorbar(os_data['centrality'], os_data['gamma'], yerr=os_data['gamma_err'], fmt=style, color=color, linewidth=2)

    # 应用之前保存的坐标范围
    axes[0].set_xlim(xlim0)
    axes[0].set_ylim(ylim0)
    axes[1].set_xlim(xlim1)
    axes[1].set_ylim(ylim1)

    axes[0].set_xlabel('Centrality')
    axes[0].set_ylabel('Delta')
    axes[0].set_title('Delta vs Centrality')
    axes[0].legend()
    axes[0].grid(True)

    axes[1].set_xlabel('Centrality')
    axes[1].set_ylabel('Gamma')
    axes[1].set_title('Gamma vs Centrality')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    output_path = f'{output_dir}/compare_{comp_file}_vs_default.pdf'
    plt.savefig(output_path)
    plt.close()