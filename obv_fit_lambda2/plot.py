import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
df = pd.read_csv('fit_obv.csv')

# 定义粒子对类型和颜色
pair_types = ['LambdaLambda', 'LambdaLambdaBar', 'LambdaBarLambda', 'LambdaBarLambdaBar']
colors = ['red', 'blue', 'green', 'orange']
markers = ['o', 's', '^', 'D']

# 创建delta参数图
fig1, axes1 = plt.subplots(2, 2, figsize=(15, 12))
fig1.suptitle('Delta Parameters vs Centrality', fontsize=16)

# Delta SS
ax = axes1[0, 0]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['delta_ss'], yerr=data['delta_ss_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Delta SS')
ax.set_title('Delta Signal-Signal')
ax.legend()
ax.grid(True, alpha=0.3)

# Delta SB
ax = axes1[0, 1]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['delta_sb'], yerr=data['delta_sb_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Delta SB')
ax.set_title('Delta Signal-Background')
ax.legend()
ax.grid(True, alpha=0.3)

# Delta BB
ax = axes1[1, 0]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['delta_bb'], yerr=data['delta_bb_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Delta BB')
ax.set_title('Delta Background-Background')
ax.legend()
ax.grid(True, alpha=0.3)

# Delta Chi2/NDF
ax = axes1[1, 1]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    chi2_ndf = data['delta_chi2'] / data['delta_ndf']
    ax.plot(data['centrality'], chi2_ndf, 
            label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Delta χ²/NDF')
ax.set_title('Delta Fit Quality')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delta_parameters_vs_centrality.png', dpi=300, bbox_inches='tight')
plt.show()

# 创建rawgamma参数图
fig2, axes2 = plt.subplots(2, 2, figsize=(15, 12))
fig2.suptitle('Rawgamma Parameters vs Centrality', fontsize=16)

# Rawgamma SS
ax = axes2[0, 0]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['rawgamma_ss'], yerr=data['rawgamma_ss_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Rawgamma SS')
ax.set_title('Rawgamma Signal-Signal')
ax.legend()
ax.grid(True, alpha=0.3)

# Rawgamma SB
ax = axes2[0, 1]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['rawgamma_sb'], yerr=data['rawgamma_sb_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Rawgamma SB')
ax.set_title('Rawgamma Signal-Background')
ax.legend()
ax.grid(True, alpha=0.3)

# Rawgamma BB
ax = axes2[1, 0]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    ax.errorbar(data['centrality'], data['rawgamma_bb'], yerr=data['rawgamma_bb_err'], 
                label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Rawgamma BB')
ax.set_title('Rawgamma Background-Background')
ax.legend()
ax.grid(True, alpha=0.3)

# Rawgamma Chi2/NDF
ax = axes2[1, 1]
for i, pair_type in enumerate(pair_types):
    data = df[df['pair_type'] == pair_type].sort_values('centrality')
    chi2_ndf = data['gamma_chi2'] / data['gamma_ndf']
    ax.plot(data['centrality'], chi2_ndf, 
            label=pair_type, color=colors[i], marker=markers[i], markersize=6, linewidth=2)
ax.set_xlabel('Centrality')
ax.set_ylabel('Rawgamma χ²/NDF')
ax.set_title('Rawgamma Fit Quality')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('rawgamma_parameters_vs_centrality.png', dpi=300, bbox_inches='tight')
plt.show()

# 打印一些统计信息
print("数据统计:")
print(f"中心度范围: {df['centrality'].min()} - {df['centrality'].max()}")
print(f"粒子对类型: {df['pair_type'].unique()}")
print(f"总数据点数: {len(df)}")

# 打印每个粒子对的delta_ss平均值
print("\nDelta SS 平均值:")
for pair_type in pair_types:
    avg_delta_ss = df[df['pair_type'] == pair_type]['delta_ss'].mean()
    print(f"{pair_type}: {avg_delta_ss:.6f}")

# 打印每个粒子对的rawgamma_ss平均值
print("\nRawgamma SS 平均值:")
for pair_type in pair_types:
    avg_gamma_ss = df[df['pair_type'] == pair_type]['rawgamma_ss'].mean()
    print(f"{pair_type}: {avg_gamma_ss:.6f}")