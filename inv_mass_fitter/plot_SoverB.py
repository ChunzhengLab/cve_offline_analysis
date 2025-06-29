import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 获取当前目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 读取18q和18r的数据
q_csv_centrality = os.path.join(current_dir, 'ratio_sigbkg_centrality_18q.csv')
r_csv_centrality = os.path.join(current_dir, 'ratio_sigbkg_centrality_18r.csv')
q_csv_pt = os.path.join(current_dir, 'ratio_sigbkg_centrality_pT_18q.csv')
r_csv_pt = os.path.join(current_dir, 'ratio_sigbkg_centrality_pT_18r.csv')

# 检查文件是否存在
files_to_check = [
    (q_csv_centrality, "18q中心度数据"),
    (r_csv_centrality, "18r中心度数据"),
    (q_csv_pt, "18q微分数据"),
    (r_csv_pt, "18r微分数据")
]

missing_files = []
for file_path, desc in files_to_check:
    if not os.path.exists(file_path):
        missing_files.append(f"{desc}({file_path})")

if missing_files:
    print(f"警告：以下文件不存在，某些图表可能无法生成：")
    for missing in missing_files:
        print(f"  - {missing}")

# 1. 绘制中心度图表
if os.path.exists(q_csv_centrality) or os.path.exists(r_csv_centrality):
    print("正在绘制信噪比随中心度变化的图表...")
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig1.suptitle('Signal/Background Ratio vs Centrality', fontsize=16)

    # 左边子图：18q数据
    if os.path.exists(q_csv_centrality):
        df_18q = pd.read_csv(q_csv_centrality)
        for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
            sub = df_18q[df_18q['particle'] == particle]
            # 去除NaN值
            sub = sub.dropna(subset=['signal_bkg_ratio'])
            # 确保按centrality排序
            sub = sub.sort_values('centrality')
            ax1.plot(sub['centrality'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)

        ax1.set_xlabel('Centrality (%)')
        ax1.set_ylabel('Signal/Background Ratio')
        ax1.set_title('18q Data')
        ax1.legend()
        ax1.grid(alpha=0.3)
    else:
        ax1.text(0.5, 0.5, "18q数据文件不存在", ha='center', va='center')
        ax1.set_title("缺少数据")

    # 右边子图：18r数据
    if os.path.exists(r_csv_centrality):
        df_18r = pd.read_csv(r_csv_centrality)
        for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
            sub = df_18r[df_18r['particle'] == particle]
            # 去除NaN值
            sub = sub.dropna(subset=['signal_bkg_ratio'])
            # 确保按centrality排序
            sub = sub.sort_values('centrality')
            ax2.plot(sub['centrality'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)

        ax2.set_xlabel('Centrality (%)')
        ax2.set_title('18r Data')
        ax2.legend()
        ax2.grid(alpha=0.3)
    else:
        ax2.text(0.5, 0.5, "18r数据文件不存在", ha='center', va='center')
        ax2.set_title("缺少数据")

    # 调整布局并保存
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(current_dir, 'SoverB_vs_centrality_18q_18r.pdf'))
    plt.savefig(os.path.join(current_dir, 'SoverB_vs_centrality_18q_18r.png'), dpi=300)
    plt.close()
    print("中心度图表已保存")

# 2. 绘制微分图表（信噪比随pT变化）
if os.path.exists(q_csv_pt) or os.path.exists(r_csv_pt):
    print("正在绘制信噪比随pT变化的图表...")
    
    # 读取微分数据
    df_18q_pt = pd.read_csv(q_csv_pt) if os.path.exists(q_csv_pt) else pd.DataFrame()
    df_18r_pt = pd.read_csv(r_csv_pt) if os.path.exists(r_csv_pt) else pd.DataFrame()
    
    # 获取所有可用的中心度值
    centrality_values = set()
    if not df_18q_pt.empty:
        centrality_values.update(df_18q_pt['centrality'].unique())
    if not df_18r_pt.empty:
        centrality_values.update(df_18r_pt['centrality'].unique())
    centrality_values = sorted(centrality_values)
    
    if not centrality_values:
        print("警告：找不到中心度值，无法绘制微分图表")
    else:
        # 选择一个中心度值进行展示，优先选择35.0或接近的值
        target_centrality = 35.0
        if target_centrality not in centrality_values:
            centrality_to_plot = min(centrality_values, key=lambda x: abs(x - target_centrality))
            print(f"注意：中心度{target_centrality}不存在，使用最接近的值{centrality_to_plot}")
        else:
            centrality_to_plot = target_centrality
        
        fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
        fig2.suptitle(f'Signal/Background Ratio vs pT (Centrality={centrality_to_plot}%)', fontsize=16)
        
        # 左边子图：18q数据
        if not df_18q_pt.empty:
            for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
                sub = df_18q_pt[(df_18q_pt['particle'] == particle) & (df_18q_pt['centrality'] == centrality_to_plot)]
                if not sub.empty:
                    sub = sub.sort_values('pT_mean')
                    ax1.plot(sub['pT_mean'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)
            
            ax1.set_xlabel('pT (GeV/c)')
            ax1.set_ylabel('Signal/Background Ratio')
            ax1.set_title('18q Data')
            ax1.legend()
            ax1.grid(alpha=0.3)
        else:
            ax1.text(0.5, 0.5, "18q微分数据文件不存在", ha='center', va='center')
            ax1.set_title("缺少数据")
        
        # 右边子图：18r数据
        if not df_18r_pt.empty:
            for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
                sub = df_18r_pt[(df_18r_pt['particle'] == particle) & (df_18r_pt['centrality'] == centrality_to_plot)]
                if not sub.empty:
                    sub = sub.sort_values('pT_mean')
                    ax2.plot(sub['pT_mean'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)
            
            ax2.set_xlabel('pT (GeV/c)')
            ax2.set_title('18r Data')
            ax2.legend()
            ax2.grid(alpha=0.3)
        else:
            ax2.text(0.5, 0.5, "18r微分数据文件不存在", ha='center', va='center')
            ax2.set_title("缺少数据")
        
        # 调整布局并保存
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(os.path.join(current_dir, f'SoverB_vs_pT_cent{int(centrality_to_plot)}_18q_18r.pdf'))
        plt.savefig(os.path.join(current_dir, f'SoverB_vs_pT_cent{int(centrality_to_plot)}_18q_18r.png'), dpi=300)
        plt.close()
        print(f"pT微分图表（中心度={centrality_to_plot}%）已保存")

# 3. 打印最大信噪比值
if os.path.exists(q_csv_centrality):
    df_18q = pd.read_csv(q_csv_centrality)
    print("\n18q 数据最大信噪比:")
    for particle in ['Lambda', 'LambdaBar']:
        particle_data = df_18q[df_18q['particle'] == particle]
        if not particle_data.empty and not particle_data['signal_bkg_ratio'].isna().all():
            max_value = particle_data['signal_bkg_ratio'].max()
            max_cent = particle_data[particle_data['signal_bkg_ratio'] == max_value]['centrality'].values[0]
            print(f"{particle}: {max_value:.2f} (中心度 {max_cent}%)")
        else:
            print(f"{particle}: 数据不可用")

if os.path.exists(r_csv_centrality):
    df_18r = pd.read_csv(r_csv_centrality)
    print("\n18r 数据最大信噪比:")
    for particle in ['Lambda', 'LambdaBar']:
        particle_data = df_18r[df_18r['particle'] == particle]
        if not particle_data.empty and not particle_data['signal_bkg_ratio'].isna().all():
            max_value = particle_data['signal_bkg_ratio'].max()
            max_cent = particle_data[particle_data['signal_bkg_ratio'] == max_value]['centrality'].values[0]
            print(f"{particle}: {max_value:.2f} (中心度 {max_cent}%)")
        else:
            print(f"{particle}: 数据不可用")

print("\n绘图完成！")