import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse
import glob

# 解析命令行参数
parser = argparse.ArgumentParser(description='拟合不变质量分布并提取信号分数')
parser.add_argument('-i', '--input-dir', type=str, default="./",
                   help='输入目录，包含按数据集和粒子类型组织的子目录')
parser.add_argument('-d', '--dataset', type=str, default="LHC18q",
                   help='数据集名称（如LHC18q，LHC18r）')
parser.add_argument('-p', '--particle', type=str, default="Proton",
                   help='粒子类型（Proton，Pion，Hadron）')
parser.add_argument('-t', '--task', type=str, default="default",
                   help='任务名称，用于定位输入文件和命名输出文件')
parser.add_argument('-o', '--output-dir', type=str, default=None,
                   help='输出目录，默认与输入目录相同')
parser.add_argument('--signal-min', type=float, default=1.103,
                   help='信号区域的最小不变质量值')
parser.add_argument('--signal-max', type=float, default=1.128,
                   help='信号区域的最大不变质量值')
parser.add_argument('--poly-order', type=int, default=2,
                   help='背景拟合使用的多项式阶数')
parser.add_argument('--fs-force-non-neg', action='store_true', default=False,
                   help='是否强制信号分数为非负值,默认不强制')
parser.add_argument('--maxfev', type=int, default=1000000,
                   help='curve_fit的最大函数评估次数')
parser.add_argument('--plot-signal-sum', action='store_true', default=False,
                   help='画出信号总量和信号/背景比值随中心度变化的图')
args = parser.parse_args()

# 设置输出目录
output_dir = args.output_dir if args.output_dir is not None else args.input_dir

# 构建输入输出路径
input_path = os.path.join(args.input_dir, args.dataset, args.particle)
input_csv = os.path.join(input_path, f"flatten_data_{args.task}.csv")

# 检查输入文件是否存在
if not os.path.exists(input_csv):
    print(f"错误：输入文件 {input_csv} 不存在")
    exit(1)

# 创建输出目录
output_path = os.path.join(output_dir, args.dataset, args.particle)
os.makedirs(output_path, exist_ok=True)

# 输出文件路径
output_csv = os.path.join(output_path, f"fit_inv_mass_{args.task}.csv")
params_csv = os.path.join(output_path, f"fit_inv_mass_pars_{args.task}.csv")
img_dir = os.path.join(output_path, f"fitplots_mass_{args.task}")
os.makedirs(img_dir, exist_ok=True)

print(f"处理数据集: {args.dataset}, 粒子类型: {args.particle}, 任务: {args.task}")
print(f"输入文件: {input_csv}")
print(f"输出目录: {output_path}")

# 定义拟合函数
def exp_poly_func(x, *params):
    """指数 + 多项式组合函数
    形式: (A)e^(kx) + (c₀) + (c₁)x + (c₂)x² + ... + (cₙ)xⁿ
    """
    A, k = params[0], params[1]  # A作为指数项的常数系数
    poly_params = params[2:]
    poly = sum(c * (x ** i) for i, c in enumerate(poly_params))
    return A * np.exp(k * x) + poly

# 信号区域
SIGNAL_MIN = args.signal_min
SIGNAL_MAX = args.signal_max

# 读取数据
print(f"读取输入文件: {input_csv}")
df = pd.read_csv(input_csv)
df['__rowid__'] = np.arange(len(df))

fit_param_records = []
fs_mass = np.zeros(len(df))

# 分组
group_cols = ['centrality', 'diff_type', 'diff_bin', 'pair_type']
full_group_cols = ['diff_type', 'diff_bin', 'pair_type']

all_groups = df.groupby(group_cols)
pdf_groups = df.groupby(full_group_cols)

print(f"找到 {len(all_groups)} 个分组")

# 拟合和fs计算
fit_param_dict = {}  # 保存参数方便后面画图调用
fit_count = 0
fail_count = 0

for group_keys, group_df in all_groups:
    side_mask = (group_df['inv_mass_center'] < SIGNAL_MIN) | (group_df['inv_mass_center'] > SIGNAL_MAX)
    x_fit = group_df.loc[side_mask, 'inv_mass_center'].values
    y_fit = group_df.loc[side_mask, 'inv_mass_counts'].values

    if len(x_fit) < args.poly_order + 3 or np.all(y_fit == 0):  # +3 是因为需要额外的两个指数参数
        params = [0.0] * (args.poly_order + 3)  # 指数参数 + 多项式参数
        fail_count += 1
    else:
        try:
            # 直接进行指数+多项式拟合
            p0 = [np.max(y_fit), -5] + [0.0] * (args.poly_order + 1)  # 初始值：A=最大值，k=-5，多项式系数全为0
            popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)
            params = popt
            fit_count += 1
        except Exception as e:
            print(f"拟合失败 {group_keys}: {str(e)}")
            params = [0.0] * (args.poly_order + 3)
            fail_count += 1

    # 保存拟合参数
    param_dict = {
        "exp_A": params[0],
        "exp_k": params[1],
        **{f"poly{i}": params[i+2] if i+2 < len(params) else 0.0 for i in range(args.poly_order + 1)}
    }
    record = dict(
        centrality=group_keys[0],
        diff_type=group_keys[1],
        diff_bin=group_keys[2],
        pair_type=group_keys[3],
        **param_dict
    )
    fit_param_records.append(record)

    # 保存拟合参数供绘图使用
    fit_param_dict[group_keys] = params

    # 计算信号分数
    group_indices = group_df['__rowid__'].values
    x_all = group_df['inv_mass_center'].values
    y_all = group_df['inv_mass_counts'].values
    bkg_all = exp_poly_func(x_all, *params)
    signal = y_all - bkg_all

    # 计算原始信号分数
    fs_group = np.zeros_like(signal)
    valid_mask = y_all > 0

    # 根据参数决定是否强制信号为非负
    if args.fs_force_non_neg:
        signal_non_neg = np.maximum(signal, 0)  # 信号不能为负
        fs_group[valid_mask] = signal_non_neg[valid_mask] / y_all[valid_mask]
    else:
        fs_group[valid_mask] = signal[valid_mask] / y_all[valid_mask]

    fs_mass[group_indices] = fs_group

print(f"拟合成功: {fit_count}, 拟合失败: {fail_count}")

# 保存结果
df['fs_mass'] = fs_mass
df.drop(columns=['__rowid__'], inplace=True)

# 记录fs处理方式
print(f"信号分数处理: {'强制非负' if args.fs_force_non_neg else '保留原始值（可能包含负值）'}")
df.to_csv(output_csv, index=False)
fit_param_df = pd.DataFrame(fit_param_records)
fit_param_df.to_csv(params_csv, index=False)
print(f"已生成 {output_csv}（含信号分数） 和 {params_csv}（拟合参数）")

# 统计信号总量和信号/背景比值
signal_sum_records = []
for group_keys, group_df in all_groups:
    centrality, diff_type, diff_bin, pair_type = group_keys
    params = fit_param_dict[group_keys]
    in_signal_mask = (group_df['inv_mass_center'] >= SIGNAL_MIN) & (group_df['inv_mass_center'] <= SIGNAL_MAX)
    y_obs = group_df.loc[in_signal_mask, 'inv_mass_counts'].values
    x = group_df.loc[in_signal_mask, 'inv_mass_center'].values
    y_bkg = exp_poly_func(x, *params)
    signal_sum = np.sum(y_obs - y_bkg)
    bkg_sum = np.sum(y_bkg)
    signal_bkg_ratio = signal_sum / bkg_sum if bkg_sum != 0 else np.nan
    signal_sum_records.append({
        'centrality': centrality,
        'diff_type': diff_type,
        'diff_bin': diff_bin,
        'pair_type': pair_type,
        'signal_sum': signal_sum,
        'bkg_sum': bkg_sum,
        'signal_bkg_ratio': signal_bkg_ratio
    })
signal_sum_df = pd.DataFrame(signal_sum_records)
signal_sum_csv = os.path.join(output_path, f'signal_sum_vs_centrality_{args.task}.csv')
signal_sum_df.to_csv(signal_sum_csv, index=False)
print(f'已生成信号总量/比值文件: {signal_sum_csv}')

# 可选画图
if args.plot_signal_sum:
    for ycol, ylab in zip(['signal_sum', 'signal_bkg_ratio'], ['Signal Sum', 'Signal/Bkg Ratio']):
        plt.figure(figsize=(8,6))
        for (diff_type, diff_bin, pair_type), sub in signal_sum_df.groupby(['diff_type', 'diff_bin', 'pair_type']):
            plt.plot(sub['centrality'], sub[ycol], 'o-', label=f'{diff_type}-{diff_bin}-{pair_type}')
        plt.xlabel('Centrality')
        plt.ylabel(ylab)
        plt.title(f'{ylab} vs Centrality')
        plt.legend(fontsize=8)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        outfig = os.path.join(output_path, f'{ycol}_vs_centrality_{args.task}.pdf')
        plt.savefig(outfig)
        print(f'已生成信号总量/比值图: {outfig}')
        plt.close()

# 绘图函数
def plot_grid_centralities(centrality_list, cendf_list, fit_param_dict, diff_type, diff_bin, pair_type):
    n = len(centrality_list)
    ncols = 3
    nrows = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 12), sharex=False, sharey=False)
    axes = axes.flatten()
    for i, (centrality, cendf) in enumerate(zip(centrality_list, cendf_list)):
        if i >= len(axes):
            break

        x = cendf['inv_mass_center'].values
        y = cendf['inv_mass_counts'].values

        # 获取拟合参数
        params = fit_param_dict.get((centrality, diff_type, diff_bin, pair_type), [0] * (args.poly_order + 3))

        ax = axes[i]
        ax.plot(x, y, 'o-', color='red', label='Data')

        # 绘制拟合曲线
        x_plot = np.linspace(x.min(), x.max(), 200)
        y_fit = exp_poly_func(x_plot, *params)
        ax.plot(x_plot, y_fit, '-', color='blue', label='Background')

        ax.set_title(f"centrality={centrality}")

        # 显示拟合参数
        param_text = f"({params[0]:.2g})e^{{{params[1]:.2g}x}}"
        for j, p in enumerate(params[2:]):
            if abs(p) > 1e-10:
                if j == 0:
                    param_text += f" + ({p:.2g})"  # 常数项
                else:
                    param_text += f" + ({p:.2g})x^{j}"
        param_text = f"${param_text}$"  # 整个公式用LaTeX格式
        ax.text(0.97, 0.95, param_text, transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend()

    # 隐藏未用子图
    for j in range(i+1, nrows*ncols):
        if j < len(axes):
            axes[j].set_visible(False)

    fig.suptitle(f"{args.dataset}/{args.particle}: diff_type={diff_type}, diff_bin={diff_bin}, pair_type={pair_type}", fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig

# 分组输出PDF，每9个centrality一页
print(f"开始生成拟合图像...")
pdf_count = 0

for (diff_type, diff_bin, pair_type), subdf in pdf_groups:
    pdfname = os.path.join(img_dir, f"fitplots_{diff_type}_{diff_bin}_{pair_type}.pdf")
    centralities = sorted(subdf['centrality'].unique())

    with PdfPages(pdfname) as pdf:
        # 分页，每9个一组
        for i in range(0, len(centralities), 9):
            centers = centralities[i:i+9]
            cendf_list = [subdf[subdf['centrality'] == c] for c in centers]
            fig = plot_grid_centralities(centers, cendf_list, fit_param_dict, diff_type, diff_bin, pair_type)
            pdf.savefig(fig)
            plt.close(fig)

    pdf_count += 1
    print(f"已输出 {pdfname}")

print(f"所有处理完成！生成了 {pdf_count} 个PDF文件")
print(f"处理模式: {'强制fs非负' if args.fs_force_non_neg else '保留原始fs值'}")
