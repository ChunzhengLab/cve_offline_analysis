import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='拟合不变质量谱并计算信号分数')
parser.add_argument('-i', '--input-dir', type=str, default='./',
                   help='输入目录')
parser.add_argument('-d', '--dataset', type=str, default=None,
                   help='数据集名称（可选，自动推断）')
parser.add_argument('-t', '--task', type=str, default='default',
                   help='任务名称')
parser.add_argument('-o', '--output-dir', type=str, default='./',
                   help='输出目录')
parser.add_argument('--fs-force-non-neg', action='store_true',
                   help='强制信号分数大于0（默认不强制）')
parser.add_argument('--signal-min', type=float, default=1.103,
                   help='信号区域的最小不变质量值')
parser.add_argument('--signal-max', type=float, default=1.128,
                   help='信号区域的最大不变质量值')
parser.add_argument('--poly-order', type=int, default=2,
                   help='背景拟合使用的多项式阶数')
parser.add_argument('--maxfev', type=int, default=1000000,
                   help='curve_fit的最大函数评估次数')
parser.add_argument('--plot-signal-sum', action='store_true',
                   help='画出信号总量随中心度变化的图')
parser.add_argument('--read-frac', action='store_true',
                   help='直接从mass_to_fs.csv读取fs_mass_1和fs_mass_2，不做拟合')
args = parser.parse_args()

# 固定粒子类型为Lambda
particle_type = 'Lambda'

# 如果没有指定数据集，尝试自动推断
if args.dataset is None:
    # 查找可能的数据集目录
    for potential_dataset in ['LHC18q', 'LHC18r']:
        if os.path.exists(os.path.join(args.input_dir, potential_dataset, particle_type)):
            args.dataset = potential_dataset
            break
    if args.dataset is None:
        args.dataset = 'unknown'

# 构建输入输出路径
input_path = os.path.join(args.input_dir, args.dataset, particle_type)
output_path = os.path.join(args.output_dir, args.dataset, particle_type)
os.makedirs(output_path, exist_ok=True)

# 设置文件路径
input_csv = os.path.join(input_path, f"flatten_data_{args.task}.csv")
output_csv = os.path.join(output_path, f"fit_inv_mass_{args.task}.csv")
params_csv = os.path.join(output_path, f"fit_inv_mass_pars_{args.task}.csv")
img_dir = os.path.join(output_path, f"fitplots_mass_{args.task}")
os.makedirs(img_dir, exist_ok=True)

print(f"输入目录: {input_path}")
print(f"输出目录: {output_path}")
print(f"数据集: {args.dataset}")
print(f"任务: {args.task}")
print(f"参数设置：强制信号分数大于0 = {args.fs_force_non_neg}")

# 定义拟合函数
def exp_poly_func(x, *params):
    """指数 + 多项式组合函数
    形式: (A)e^(kx) + (c₀) + (c₁)x + (c₂)x² + ... + (cₙ)xⁿ
    """
    A, k = params[0], params[1]
    poly_params = params[2:]
    poly = sum(c * (x ** i) for i, c in enumerate(poly_params))
    return A * np.exp(k * x) + poly

SIGNAL_MIN = args.signal_min
SIGNAL_MAX = args.signal_max

df = pd.read_csv(input_csv)
df['__rowid__'] = np.arange(len(df))

# 分组：只按centrality和pair_type分组
group_cols = ['centrality', 'pair_type']
full_group_cols = ['pair_type']

all_groups = df.groupby(group_cols)
pdf_groups = df.groupby(full_group_cols)

# 存储积分结果和拟合参数
integrated_data = {}  # 存储积分后的数据
fit_param_dict = {}  # 存储拟合参数
fit_param_records = []

# 1. 首先计算积分结果
for group_keys, group_df in all_groups:
    centrality, pair_type = group_keys
    # 确保centrality是float类型，pair_type是字符串
    centrality = float(centrality)
    pair_type = str(pair_type)

    # 对inv_mass_2进行积分，得到inv_mass_1的分布
    mass1_df = group_df.groupby('inv_mass_1_center').agg({
        'inv_mass_counts': 'sum'
    }).reset_index()
    mass1_df.columns = ['inv_mass_1_center', 'inv_mass_counts_intg2']

    # 对inv_mass_1进行积分，得到inv_mass_2的分布
    mass2_df = group_df.groupby('inv_mass_2_center').agg({
        'inv_mass_counts': 'sum'
    }).reset_index()
    mass2_df.columns = ['inv_mass_2_center', 'inv_mass_counts_intg1']

    # 保存积分结果
    integrated_data[(centrality, pair_type)] = {
        'mass1': mass1_df,
        'mass2': mass2_df
    }

# 添加调试信息，显示存储的键
print("integrated_data中存储的键:")
for key in integrated_data.keys():
    print(f"  {key} (centrality type: {type(key[0])}, pair_type type: {type(key[1])})")

# 2. 对积分结果进行拟合
for group_keys, data in integrated_data.items():
    centrality, pair_type = group_keys
    mass1_df = data['mass1']
    mass2_df = data['mass2']

    # 拟合inv_mass_1
    side_mask_1 = (mass1_df['inv_mass_1_center'] < SIGNAL_MIN) | (mass1_df['inv_mass_1_center'] > SIGNAL_MAX)
    x_fit_1 = mass1_df.loc[side_mask_1, 'inv_mass_1_center'].values
    y_fit_1 = mass1_df.loc[side_mask_1, 'inv_mass_counts_intg2'].values

    if len(x_fit_1) < args.poly_order + 3 or np.all(y_fit_1 == 0):
        params1 = [0.0] * (args.poly_order + 3)
    else:
        try:
            # 直接进行指数+多项式拟合
            p0_1 = [np.max(y_fit_1), -5] + [0.0] * (args.poly_order + 1)  # 初始值：A=最大值，k=-5，多项式系数全为0
            popt1, _ = curve_fit(exp_poly_func, x_fit_1, y_fit_1, p0=p0_1, maxfev=args.maxfev)
            params1 = popt1
        except Exception as e:
            print(f"拟合失败 mass1 {group_keys}: {str(e)}")
            params1 = [0.0] * (args.poly_order + 3)

    # 拟合inv_mass_2
    side_mask_2 = (mass2_df['inv_mass_2_center'] < SIGNAL_MIN) | (mass2_df['inv_mass_2_center'] > SIGNAL_MAX)
    x_fit_2 = mass2_df.loc[side_mask_2, 'inv_mass_2_center'].values
    y_fit_2 = mass2_df.loc[side_mask_2, 'inv_mass_counts_intg1'].values

    if len(x_fit_2) < args.poly_order + 3 or np.all(y_fit_2 == 0):
        params2 = [0.0] * (args.poly_order + 3)
    else:
        try:
            # 直接进行指数+多项式拟合
            p0_2 = [np.max(y_fit_2), -5] + [0.0] * (args.poly_order + 1)  # 初始值：A=最大值，k=-5，多项式系数全为0
            popt2, _ = curve_fit(exp_poly_func, x_fit_2, y_fit_2, p0=p0_2, maxfev=args.maxfev)
            params2 = popt2
        except Exception as e:
            print(f"拟合失败 mass2 {group_keys}: {str(e)}")
            params2 = [0.0] * (args.poly_order + 3)

    # 保存拟合参数
    param_dict = {
        'exp_A_mass_1': params1[0],
        'exp_k_mass_1': params1[1],
        **{f'poly{i}_mass_1': params1[i+2] for i in range(args.poly_order + 1)},
        'exp_A_mass_2': params2[0],
        'exp_k_mass_2': params2[1],
        **{f'poly{i}_mass_2': params2[i+2] for i in range(args.poly_order + 1)}
    }
    record = dict(
        centrality=centrality,
        pair_type=pair_type,
        **param_dict
    )
    fit_param_records.append(record)
    fit_param_dict[(centrality, pair_type)] = (params1, params2)

# 3. 计算信号分数
if args.read_frac:
    print('开启 --read-frac，直接从 mass_to_fs.csv 读取 fs_mass_1 和 fs_mass_2')
    frac_df = pd.read_csv('mass_to_fs.csv')
    # 构建查找表：centrality, particle, inv_mass_center -> fs_mass
    frac_lookup = frac_df.set_index(['centrality', 'particle', 'inv_mass_center'])['fs_mass']
    # 需要保证 centrality 类型一致
    def get_fs(cent, particle, mass):
        try:
            return frac_lookup.loc[(cent, particle, mass)]
        except KeyError:
            return np.nan
    fs_mass_1 = np.zeros(len(df))
    fs_mass_2 = np.zeros(len(df))
    for idx, row in df.iterrows():
        cent = row['centrality']
        mass1 = row['inv_mass_1_center']
        mass2 = row['inv_mass_2_center']
        pair_type = row['pair_type']
        # 逻辑判断
        if pair_type == 'LambdaLambda':
            fs_mass_1[idx] = get_fs(cent, 'Lambda', mass1)
            fs_mass_2[idx] = get_fs(cent, 'Lambda', mass2)
        elif pair_type == 'LambdaLambdaBar':
            fs_mass_1[idx] = get_fs(cent, 'Lambda', mass1)
            fs_mass_2[idx] = get_fs(cent, 'LambdaBar', mass2)
        elif pair_type == 'LambdaBarLambda':
            fs_mass_1[idx] = get_fs(cent, 'LambdaBar', mass1)
            fs_mass_2[idx] = get_fs(cent, 'Lambda', mass2)
        elif pair_type == 'LambdaBarLambdaBar':
            fs_mass_1[idx] = get_fs(cent, 'LambdaBar', mass1)
            fs_mass_2[idx] = get_fs(cent, 'LambdaBar', mass2)
        else:
            fs_mass_1[idx] = np.nan
            fs_mass_2[idx] = np.nan
    df['fs_mass_1'] = fs_mass_1
    df['fs_mass_2'] = fs_mass_2
    df.drop(columns=['__rowid__'], inplace=True)
    df.to_csv(output_csv, index=False)
    print(f"已生成 {output_csv}（直接读取fs_mass_1/2）")
    # 后续流程（信号总量统计/画图）可按需要保留或跳过
    exit(0)

fs_mass_1 = np.zeros(len(df))
fs_mass_2 = np.zeros(len(df))

for group_keys, group_df in all_groups:
    centrality, pair_type = group_keys
    # 确保centrality是float类型，pair_type是字符串
    centrality = float(centrality)
    pair_type = str(pair_type)

    group_indices = group_df['__rowid__'].values

    # 获取积分数据
    mass1_df = integrated_data[(centrality, pair_type)]['mass1']
    mass2_df = integrated_data[(centrality, pair_type)]['mass2']

    # 创建查找字典
    mass1_dict = dict(zip(mass1_df['inv_mass_1_center'], mass1_df['inv_mass_counts_intg2']))
    mass2_dict = dict(zip(mass2_df['inv_mass_2_center'], mass2_df['inv_mass_counts_intg1']))

    # 获取拟合参数
    params = fit_param_dict[(centrality, pair_type)]
    params1, params2 = params

    # 计算mass_1的信号分数
    x_all_1 = group_df['inv_mass_1_center'].values
    y_all_1 = np.array([mass1_dict.get(x, 0) for x in x_all_1])  # 使用积分后的计数
    bkg_all_1 = exp_poly_func(x_all_1, *params1)
    signal_1 = y_all_1 - bkg_all_1

    # 根据参数决定是否强制信号为正
    if args.fs_force_non_neg:
        signal_1 = np.maximum(signal_1, 0)

    fs_group_1 = np.zeros_like(signal_1)
    valid_mask_1 = y_all_1 > 0
    fs_group_1[valid_mask_1] = signal_1[valid_mask_1] / y_all_1[valid_mask_1]
    fs_mass_1[group_indices] = fs_group_1

    # 计算mass_2的信号分数
    x_all_2 = group_df['inv_mass_2_center'].values
    y_all_2 = np.array([mass2_dict.get(x, 0) for x in x_all_2])  # 使用积分后的计数
    bkg_all_2 = exp_poly_func(x_all_2, *params2)
    signal_2 = y_all_2 - bkg_all_2

    # 根据参数决定是否强制信号为正
    if args.fs_force_non_neg:
        signal_2 = np.maximum(signal_2, 0)

    fs_group_2 = np.zeros_like(signal_2)
    valid_mask_2 = y_all_2 > 0
    fs_group_2[valid_mask_2] = signal_2[valid_mask_2] / y_all_2[valid_mask_2]
    fs_mass_2[group_indices] = fs_group_2

df['fs_mass_1'] = fs_mass_1
df['fs_mass_2'] = fs_mass_2
df.drop(columns=['__rowid__'], inplace=True)
df.to_csv(output_csv, index=False)
fit_param_df = pd.DataFrame(fit_param_records)
fit_param_df.to_csv(params_csv, index=False)
print(f"已生成 {output_csv}（含信号分数，{'强制大于0' if args.fs_force_non_neg else '允许为负值'}） 和 {params_csv}（拟合参数）")

# 统计信号总量（每个中心度、pair_type、mass1和mass2）
signal_sum_records = []
for (centrality, pair_type), data in integrated_data.items():
    # mass1
    mass1_df = data['mass1']
    params1 = fit_param_dict[(centrality, pair_type)][0]
    in_signal_mask1 = (mass1_df['inv_mass_1_center'] >= SIGNAL_MIN) & (mass1_df['inv_mass_1_center'] <= SIGNAL_MAX)
    y_obs1 = mass1_df.loc[in_signal_mask1, 'inv_mass_counts_intg2'].values
    x1 = mass1_df.loc[in_signal_mask1, 'inv_mass_1_center'].values
    y_bkg1 = exp_poly_func(x1, *params1)
    signal_sum1 = np.sum(y_obs1 - y_bkg1)
    bkg_sum1 = np.sum(y_bkg1)
    signal_bkg_ratio1 = signal_sum1 / bkg_sum1 if bkg_sum1 != 0 else np.nan
    # mass2
    mass2_df = data['mass2']
    params2 = fit_param_dict[(centrality, pair_type)][1]
    in_signal_mask2 = (mass2_df['inv_mass_2_center'] >= SIGNAL_MIN) & (mass2_df['inv_mass_2_center'] <= SIGNAL_MAX)
    y_obs2 = mass2_df.loc[in_signal_mask2, 'inv_mass_counts_intg1'].values
    x2 = mass2_df.loc[in_signal_mask2, 'inv_mass_2_center'].values
    y_bkg2 = exp_poly_func(x2, *params2)
    signal_sum2 = np.sum(y_obs2 - y_bkg2)
    bkg_sum2 = np.sum(y_bkg2)
    signal_bkg_ratio2 = signal_sum2 / bkg_sum2 if bkg_sum2 != 0 else np.nan
    signal_sum_records.append({
        'centrality': centrality,
        'pair_type': pair_type,
        'signal_sum_mass1': signal_sum1,
        'signal_sum_mass2': signal_sum2,
        'signal_bkg_ratio_mass1': signal_bkg_ratio1,
        'signal_bkg_ratio_mass2': signal_bkg_ratio2
    })
signal_sum_df = pd.DataFrame(signal_sum_records)
signal_sum_csv = os.path.join(output_path, f'signal_sum_vs_centrality_{args.task}.csv')
signal_sum_df.to_csv(signal_sum_csv, index=False)
print(f"已生成信号总量文件: {signal_sum_csv}")

# 如果flag打开，画信号总量和信号/背景比值随中心度变化的图
if args.plot_signal_sum:
    for mass_type, ylab in zip(['signal_sum_mass1', 'signal_sum_mass2'], ['Signal Sum (mass1)', 'Signal Sum (mass2)']):
        plt.figure(figsize=(8,6))
        for pair_type in signal_sum_df['pair_type'].unique():
            sub = signal_sum_df[signal_sum_df['pair_type'] == pair_type]
            plt.plot(sub['centrality'], sub[mass_type], 'o-', label=pair_type)
        plt.xlabel('Centrality')
        plt.ylabel(ylab)
        plt.title(f'{ylab} vs Centrality')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        outfig = os.path.join(output_path, f'{mass_type}_vs_centrality_{args.task}.pdf')
        plt.savefig(outfig)
        print(f"已生成信号总量图: {outfig}")
        plt.close()
    # 画信号/背景比值
    for mass_type, ylab in zip(['signal_bkg_ratio_mass1', 'signal_bkg_ratio_mass2'], ['Signal/Bkg Ratio (mass1)', 'Signal/Bkg Ratio (mass2)']):
        plt.figure(figsize=(8,6))
        for pair_type in signal_sum_df['pair_type'].unique():
            sub = signal_sum_df[signal_sum_df['pair_type'] == pair_type]
            plt.plot(sub['centrality'], sub[mass_type], 'o-', label=pair_type)
        plt.xlabel('Centrality')
        plt.ylabel(ylab)
        plt.title(f'{ylab} vs Centrality')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        outfig = os.path.join(output_path, f'{mass_type}_vs_centrality_{args.task}.pdf')
        plt.savefig(outfig)
        print(f"已生成信号/背景比值图: {outfig}")
        plt.close()

def plot_grid_centralities(centrality_list, cendf_list, fit_param_dict, integrated_data, pair_type_key, mass_type):
    n = len(centrality_list)
    ncols = 3
    nrows = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 12), sharex=False, sharey=False)
    axes = axes.flatten()

    for i, (centrality, cendf) in enumerate(zip(centrality_list, cendf_list)):
        # 确保centrality是float类型
        centrality = float(centrality)
        # pair_type_key已经是处理过的格式

        # 使用积分后的数据
        if mass_type == 1:
            mass_df = integrated_data[(centrality, pair_type_key)]['mass1']
            x = mass_df['inv_mass_1_center'].values
            y = mass_df['inv_mass_counts_intg2'].values
            params = fit_param_dict.get((centrality, pair_type_key), ([0]*(args.poly_order+3), [0]*(args.poly_order+3)))
            params = params[0]
        else:
            mass_df = integrated_data[(centrality, pair_type_key)]['mass2']
            x = mass_df['inv_mass_2_center'].values
            y = mass_df['inv_mass_counts_intg1'].values
            params = fit_param_dict.get((centrality, pair_type_key), ([0]*(args.poly_order+3), [0]*(args.poly_order+3)))
            params = params[1]

        ax = axes[i]
        ax.plot(x, y, 'o-', color='red', label='Data')
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

    # Hide unused subplots
    for j in range(i+1, nrows*ncols):
        axes[j].set_visible(False)
    fig.suptitle(f"pair_type={pair_type_key}, mass_{mass_type}", fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig

# 分组输出PDF，每9个centrality一页
print(f"开始生成拟合图像...")
pdf_count = 0

for pair_type, subdf in pdf_groups:
    # 确保pair_type格式与存储时一致
    if isinstance(pair_type, tuple):
        pair_type_key = str(pair_type[0]) if len(pair_type) == 1 else str(pair_type)
    else:
        pair_type_key = str(pair_type)
    
    print(f"处理 pair_type: {pair_type} -> 使用key: {pair_type_key}")
    
    # 生成mass_1的图
    pdfname1 = os.path.join(img_dir, f"fitplots_mass1_{pair_type_key}.pdf")
    centralities = sorted(subdf['centrality'].unique())

    with PdfPages(pdfname1) as pdf:
        # 分页，每9个一组
        for i in range(0, len(centralities), 9):
            centers = centralities[i:i+9]
            cendf_list = [subdf[subdf['centrality'] == c] for c in centers]
            fig = plot_grid_centralities(centers, cendf_list, fit_param_dict, integrated_data, pair_type_key, 1)
            pdf.savefig(fig)
            plt.close(fig)

    # 生成mass_2的图
    pdfname2 = os.path.join(img_dir, f"fitplots_mass2_{pair_type_key}.pdf")
    with PdfPages(pdfname2) as pdf:
        # 分页，每9个一组
        for i in range(0, len(centralities), 9):
            centers = centralities[i:i+9]
            cendf_list = [subdf[subdf['centrality'] == c] for c in centers]
            fig = plot_grid_centralities(centers, cendf_list, fit_param_dict, integrated_data, pair_type_key, 2)
            pdf.savefig(fig)
            plt.close(fig)

    pdf_count += 2
    print(f"已输出 {pdfname1} 和 {pdfname2}")

print(f"所有处理完成！生成了 {pdf_count} 个PDF文件")
print(f"处理模式: {'强制fs非负' if args.fs_force_non_neg else '保留原始fs值'}")
