import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse
import glob

# 解析命令行参数
parser = argparse.ArgumentParser(description='拟合不变质量分布并提取信号分数（新格式）')
parser.add_argument('-i', '--input-csv', type=str, required=True, help='输入CSV文件，由flatten_data.py生成')
parser.add_argument('-o', '--output-csv', type=str, required=True, help='输出CSV文件，含信号分数')
parser.add_argument('--signal-min', type=float, default=1.103, help='信号区域的最小不变质量值')
parser.add_argument('--signal-max', type=float, default=1.128, help='信号区域的最大不变质量值')
parser.add_argument('--poly-order', type=int, default=2, help='背景拟合使用的多项式阶数')
parser.add_argument('--fs-force-non-neg', action='store_true', default=False, help='是否强制信号分数为非负值,默认不强制')
parser.add_argument('--maxfev', type=int, default=1000000, help='curve_fit的最大函数评估次数')
parser.add_argument('--plot-signal-sum', action='store_true', default=False,
                   help='画出信号总量和信号/背景比值随中心度变化的图')
args = parser.parse_args()

# 判断数据集后缀
input_filename = os.path.basename(args.input_csv)
if '18q' in input_filename:
    dataset_suffix = '_18q'
elif '18r' in input_filename:
    dataset_suffix = '_18r'
else:
    dataset_suffix = '_unknowndata'

# 定义拟合函数
def exp_poly_func(x, *params):
    A, k = params[0], params[1]
    poly_params = params[2:]
    poly = sum(c * (x ** i) for i, c in enumerate(poly_params))
    return A * np.exp(k * x) + poly

SIGNAL_MIN = args.signal_min
SIGNAL_MAX = args.signal_max

# 读取数据
print(f"读取输入文件: {args.input_csv}")
df = pd.read_csv(args.input_csv)
df['__rowid__'] = np.arange(len(df))

# 分组
# 每个分组是 (centrality, particle, pT_mean)
group_cols = ['centrality', 'particle', 'pT_mean']
all_groups = df.groupby(group_cols)
print(f"找到 {len(all_groups)} 个分组")

fs_mass = np.zeros(len(df))
fit_count = 0
fail_count = 0
fit_param_dict = {}  # 新增：保存每组的拟合参数

for group_keys, group_df in all_groups:
    side_mask = (group_df['inv_mass_mean'] < SIGNAL_MIN) | (group_df['inv_mass_mean'] > SIGNAL_MAX)
    x_fit = group_df.loc[side_mask, 'inv_mass_mean'].values
    y_fit = group_df.loc[side_mask, 'inv_mass_count'].values

    if len(x_fit) < args.poly_order + 3 or np.all(y_fit == 0):
        params = [0.0] * (args.poly_order + 3)
        fail_count += 1
    else:
        try:
            p0 = [np.max(y_fit), -5] + [0.0] * (args.poly_order + 1)
            popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)
            params = popt
            fit_count += 1
        except Exception as e:
            print(f"拟合失败 {group_keys}: {str(e)}")
            params = [0.0] * (args.poly_order + 3)
            fail_count += 1

    fit_param_dict[group_keys] = params  # 新增：保存参数

    group_indices = group_df['__rowid__'].values
    x_all = group_df['inv_mass_mean'].values
    y_all = group_df['inv_mass_count'].values
    bkg_all = exp_poly_func(x_all, *params)
    signal = y_all - bkg_all
    fs_group = np.zeros_like(signal)
    valid_mask = y_all > 0
    if args.fs_force_non_neg:
        signal_non_neg = np.maximum(signal, 0)
        fs_group[valid_mask] = signal_non_neg[valid_mask] / y_all[valid_mask]
    else:
        fs_group[valid_mask] = signal[valid_mask] / y_all[valid_mask]
    fs_mass[group_indices] = fs_group

print(f"拟合成功: {fit_count}, 拟合失败: {fail_count}")

df['fs'] = fs_mass
df.drop(columns=['__rowid__'], inplace=True)

# 只保留需要的列
out_cols = ['centrality', 'particle', 'pT_mean', 'inv_mass_mean', 'fs']
# 获取粒子名称，用于文件命名
particle_list = sorted(df['particle'].unique())
particle_str = '_'.join(particle_list)
# 使用输出目录
output_dir = os.path.dirname(args.output_csv)
# 保存主输出文件 - 只生成一个合并的文件
output_csv = f"inv_mass_frac_sig{dataset_suffix}.csv"
output_path = os.path.join(output_dir, output_csv)
df[out_cols].to_csv(output_path, index=False)
print(f"已生成 {output_path}（含信号分数）")

# 保存用于后续处理的完整路径
output_base_path = os.path.join(output_dir, f"inv_mass_frac_sig{dataset_suffix}")

# 统计信号总量和信号/背景比值
signal_sum_records = []
for group_keys, group_df in all_groups:
    centrality, particle, pT_mean = group_keys
    # 将centrality元组的第一个元素作为浮点数
    cent_value = centrality[0] if isinstance(centrality, tuple) else centrality
    params = fit_param_dict.get(group_keys, [0.0] * (args.poly_order + 3))
    in_signal_mask = (group_df['inv_mass_mean'] >= SIGNAL_MIN) & (group_df['inv_mass_mean'] <= SIGNAL_MAX)
    y_obs = group_df.loc[in_signal_mask, 'inv_mass_count'].values
    x = group_df.loc[in_signal_mask, 'inv_mass_mean'].values
    y_bkg = exp_poly_func(x, *params)
    signal_sum = np.sum(y_obs - y_bkg)
    bkg_sum = np.sum(y_bkg)
    signal_bkg_ratio = signal_sum / bkg_sum if bkg_sum != 0 else np.nan
    signal_sum_records.append({
        'centrality': cent_value,
        'particle': particle,
        'pT_mean': pT_mean,
        'signal_sum': signal_sum,
        'bkg_sum': bkg_sum,
        'signal_bkg_ratio': signal_bkg_ratio
    })
signal_sum_df = pd.DataFrame(signal_sum_records)

# 保存微分数据文件 (包含pT信息)
diff_signal_sum_csv = os.path.join(output_dir, f'ratio_sigbkg_centrality_pT{dataset_suffix}.csv')
signal_sum_df.to_csv(diff_signal_sum_csv, index=False)
print(f'已生成微分信号总量/比值文件 (包含pT): {diff_signal_sum_csv}')

# 按中心度和粒子类型进行汇总，去除pT维度
agg_columns = ['centrality', 'particle']
agg_signal_sum_df = signal_sum_df.groupby(agg_columns, as_index=False).agg({
    'signal_sum': 'sum',
    'bkg_sum': 'sum'
})
# 重新计算汇总后的信噪比
agg_signal_sum_df['signal_bkg_ratio'] = agg_signal_sum_df.apply(
    lambda row: row['signal_sum'] / row['bkg_sum'] if row['bkg_sum'] != 0 else np.nan, 
    axis=1
)

# 保存汇总数据文件 (不包含pT信息)
agg_signal_sum_csv = os.path.join(output_dir, f'ratio_sigbkg_centrality{dataset_suffix}.csv')
agg_signal_sum_df.to_csv(agg_signal_sum_csv, index=False)
print(f'已生成汇总信号总量/比值文件: {agg_signal_sum_csv}')

# 可选画图
if args.plot_signal_sum:
    for ycol, ylab in zip(['signal_sum', 'signal_bkg_ratio'], ['Signal Sum', 'Signal/Background Ratio']):
        plt.figure(figsize=(8,6))
        for (centrality, particle, pT_mean), sub in signal_sum_df.groupby(['centrality', 'particle', 'pT_mean']):
            plt.plot(sub['pT_mean'], sub[ycol], 'o-', label=f'{centrality}-{particle}-{pT_mean}')
        plt.xlabel('pT (GeV/c)')
        plt.ylabel(ylab)
        plt.title(f'{ylab} vs pT (Centrality%)')
        plt.legend(fontsize=8)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        outfig = os.path.join(output_dir, f'{ycol}_vs_pT_mean_ALL{dataset_suffix}.pdf')
        plt.savefig(outfig)
        print(f'已生成信号总量/比值图: {outfig}')
        plt.close()

# 分组输出PDF，每个(centrality, particle)一个PDF，3x6子图
print(f"开始生成拟合图像...")
plot_dir = os.path.join(output_dir, f'plot{dataset_suffix}')
os.makedirs(plot_dir, exist_ok=True)
pdf_count = 0

# 获取所有unique (centrality, particle)
cent_part_groups = df.groupby(['centrality', 'particle'])

for (centrality, particle), subdf in cent_part_groups:
    # 按pT排序
    pt_bins = sorted(subdf['pT_mean'].unique())
    n_pt = len(pt_bins)
    ncols = 6
    nrows = 3
    pdfname = os.path.join(plot_dir, f"fitplots_cent{centrality}_{particle}.pdf")
    with PdfPages(pdfname) as pdf:
        # ----------- pT分布图（无公式） -----------
        for i in range(0, n_pt, nrows * ncols):
            pt_this_page = pt_bins[i:i + nrows * ncols]
            fig, axes = plt.subplots(nrows, ncols, figsize=(18, 9), sharex=False, sharey=False)
            axes = axes.flatten()
            for j, pT in enumerate(pt_this_page):
                ax = axes[j]
                g = subdf[subdf['pT_mean'] == pT]
                x = g['inv_mass_mean'].values
                y = g['inv_mass_count'].values
                params = fit_param_dict.get((centrality, particle, pT), [0] * (args.poly_order + 3))
                ax.plot(x, y, 'o-', color='red', label='Data')
                x_plot = np.linspace(x.min(), x.max(), 200)
                y_fit = exp_poly_func(x_plot, *params)
                ax.plot(x_plot, y_fit, '-', color='blue', label='Background')
                ax.set_title(f"pT={pT:.2f}")
                # 不再显示公式
                ax.grid(True, alpha=0.3)
                if j == 0:
                    ax.legend()
            # 隐藏多余子图
            for j in range(len(pt_this_page), nrows * ncols):
                axes[j].set_visible(False)
            fig.suptitle(f"centrality={centrality}, particle={particle}", fontsize=16)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)
        # ----------- pT积分后的拟合图（2x3，带公式，单独PDF，英文） -----------
        # 先收集所有particle
        all_particles = sorted(df['particle'].unique())
        ptint_groups = df.groupby(['centrality', 'particle', 'inv_mass_mean'], as_index=False)['inv_mass_count'].sum()
        cent_list = sorted(df['centrality'].unique())[:6]
        for particle in all_particles:
            fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=False, sharey=False)
            axes = axes.flatten()
            for idx, cent in enumerate(cent_list):
                ax = axes[idx]
                group_df = ptint_groups[(ptint_groups['centrality'] == cent) & (ptint_groups['particle'] == particle)]
                if group_df.empty:
                    ax.set_visible(False)
                    continue
                x = group_df['inv_mass_mean'].values
                y = group_df['inv_mass_count'].values
                # 侧带用于拟合
                side_mask = (x < SIGNAL_MIN) | (x > SIGNAL_MAX)
                x_fit = x[side_mask]
                y_fit = y[side_mask]
                if len(x_fit) < args.poly_order + 3 or np.all(y_fit == 0):
                    params = [0.0] * (args.poly_order + 3)
                else:
                    try:
                        p0 = [np.max(y_fit), -5] + [0.0] * (args.poly_order + 1)
                        popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)
                        params = popt
                    except Exception as e:
                        print(f"pT-integrated fit failed {(cent, particle)}: {str(e)}")
                        params = [0.0] * (args.poly_order + 3)
                ax.plot(x, y, 'o-', color='red', label='Data')
                x_plot = np.linspace(x.min(), x.max(), 200)
                y_fit = exp_poly_func(x_plot, *params)
                ax.plot(x_plot, y_fit, '-', color='blue', label='Background')
                ax.set_title(f"Centrality={cent}")
                # 显示公式
                param_text = f"({params[0]:.2g})e^{{{params[1]:.2g}x}}"
                for k, p in enumerate(params[2:]):
                    if abs(p) > 1e-10:
                        if k == 0:
                            param_text += f" + ({p:.2g})"
                        else:
                            param_text += f" + ({p:.2g})x^{k}"
                param_text = f"${param_text}$"
                ax.text(0.97, 0.95, param_text, transform=ax.transAxes, ha='right', va='top', fontsize=8,
                        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
                ax.grid(True, alpha=0.3)
                if idx == 0:
                    ax.legend()
            for idx in range(len(cent_list), 6):
                axes[idx].set_visible(False)
            fig.suptitle(f"pT-integrated fit: particle={particle}", fontsize=16)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            ptint_pdf = os.path.join(plot_dir, f"fitplots_ptint_{particle}.pdf")
            fig.savefig(ptint_pdf)
            plt.close(fig)
            print(f"Output pT-integrated PDF: {ptint_pdf}")

        # ----------- 生成pT积分后的数据文件 -----------
        # 收集所有粒子的数据到一个列表中
        all_ptint_signal_sum_records = []
        
        for particle_type in all_particles:
            # 筛选当前粒子类型的记录
            particle_groups = ptint_groups[ptint_groups['particle'] == particle_type].groupby(['centrality'])
            
            for centrality, group_df in particle_groups:
                x = group_df['inv_mass_mean'].values
                y = group_df['inv_mass_count'].values
                side_mask = (x < SIGNAL_MIN) | (x > SIGNAL_MAX)
                x_fit = x[side_mask]
                y_fit = y[side_mask]
                if len(x_fit) < args.poly_order + 3 or np.all(y_fit == 0):
                    params = [0.0] * (args.poly_order + 3)
                else:
                    try:
                        p0 = [np.max(y_fit), -5] + [0.0] * (args.poly_order + 1)
                        popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)
                        params = popt
                    except Exception as e:
                        print(f"pT-integrated fit failed {(centrality, particle_type)}: {str(e)}")
                        params = [0.0] * (args.poly_order + 3)
                in_signal_mask = (x >= SIGNAL_MIN) & (x <= SIGNAL_MAX)
                y_obs = y[in_signal_mask]
                x_sig = x[in_signal_mask]
                y_bkg = exp_poly_func(x_sig, *params)
                signal_sum = np.sum(y_obs - y_bkg)
                bkg_sum = np.sum(y_bkg)
                signal_bkg_ratio = signal_sum / bkg_sum if bkg_sum != 0 else np.nan
                # 将centrality元组的第一个元素作为浮点数
                cent_value = centrality[0] if isinstance(centrality, tuple) else centrality
                all_ptint_signal_sum_records.append({
                    'centrality': cent_value,
                    'particle': particle_type,
                    'signal_sum': signal_sum,
                    'bkg_sum': bkg_sum,
                    'signal_bkg_ratio': signal_bkg_ratio
                })
        
        # 仅生成一个包含所有粒子的CSV文件
        if all_ptint_signal_sum_records:  # 确保有数据再生成文件
            all_ptint_signal_sum_df = pd.DataFrame(all_ptint_signal_sum_records)
            # 这部分数据与前面的汇总文件重复，所以不再输出到文件
            # 只用于生成PDF图表
            pass
    pdf_count += 1
    print(f"已输出 {pdfname}")

print(f"所有处理完成！生成了 {pdf_count} 个PDF文件")
print(f"处理模式: {'强制fs非负' if args.fs_force_non_neg else '保留原始fs值'}")

# 新功能：对pT整体积分后再fit
print("开始对pT整体积分后的分布进行拟合和信号分数计算...")
ptint_groups = df.groupby(['centrality', 'particle', 'inv_mass_mean'], as_index=False)['inv_mass_count'].sum()

# 重新分组：每个(centrality, particle)一组
ptint_all_groups = ptint_groups.groupby(['centrality', 'particle'])

ptint_records = []
for (centrality, particle), group_df in ptint_all_groups:
    # 将centrality元组的第一个元素作为浮点数
    cent_value = centrality[0] if isinstance(centrality, tuple) else centrality
    x = group_df['inv_mass_mean'].values
    y = group_df['inv_mass_count'].values
    # 侧带用于拟合
    side_mask = (x < SIGNAL_MIN) | (x > SIGNAL_MAX)
    x_fit = x[side_mask]
    y_fit = y[side_mask]
    if len(x_fit) < args.poly_order + 3 or np.all(y_fit == 0):
        params = [0.0] * (args.poly_order + 3)
    else:
        try:
            p0 = [np.max(y_fit), -5] + [0.0] * (args.poly_order + 1)
            popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)
            params = popt
        except Exception as e:
            print(f"pT整体积分拟合失败 {(centrality, particle)}: {str(e)}")
            params = [0.0] * (args.poly_order + 3)
    bkg_all = exp_poly_func(x, *params)
    signal = y - bkg_all
    fs = np.zeros_like(signal)
    valid_mask = y > 0
    if args.fs_force_non_neg:
        signal_non_neg = np.maximum(signal, 0)
        fs[valid_mask] = signal_non_neg[valid_mask] / y[valid_mask]
    else:
        fs[valid_mask] = signal[valid_mask] / y[valid_mask]
    for inv_mass_mean, fs_val in zip(x, fs):
        ptint_records.append({
            'centrality': cent_value,
            'particle': particle,
            'inv_mass_mean': inv_mass_mean,
            'fs': fs_val
        })

ptint_df = pd.DataFrame(ptint_records)

# 只保存一个包含所有粒子的pT整合文件
ptint_csv = os.path.join(output_dir, f'inv_mass_frac_sig{dataset_suffix}_ptint.csv')
ptint_df.to_csv(ptint_csv, index=False)
print(f"已生成pT整体积分后的信号分数表: {ptint_csv}")
