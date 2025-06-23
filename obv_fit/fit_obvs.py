import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, OptimizeWarning
from matplotlib.backends.backend_pdf import PdfPages
import os
import argparse
import glob
import warnings

# 解析命令行参数
parser = argparse.ArgumentParser(description='拟合不变质量分布的观测量并提取信号和背景贡献')
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
parser.add_argument('--fs-min', type=float, default=None,
                   help='信号分数的最小值过滤（可选）')
parser.add_argument('--fs-max', type=float, default=None,
                   help='信号分数的最大值过滤（可选）')
parser.add_argument('--fit-method', type=str, choices=['lm', 'trf', 'dogbox'], default='trf',
                   help='拟合使用的优化方法: lm (Levenberg-Marquardt), trf (Trust Region Reflective), dogbox (dogleg)')
parser.add_argument('--range-mass', type=float, default=0.0014,
                   help='质量窗口范围半宽度，表示中心质量左右各±该值 (默认: 0.0014 GeV/c²)')
parser.add_argument('--no-fit', action='store_true',
                   help='不进行拟合，直接在质量窗口范围内计算加权平均值')
args = parser.parse_args()

# 设置输出目录
output_dir = args.output_dir if args.output_dir is not None else args.input_dir

# 构建输入输出路径
input_path = os.path.join(args.input_dir, args.dataset, args.particle)
input_csv = os.path.join(input_path, f"fit_inv_mass_{args.task}.csv")

# 检查输入文件是否存在
if not os.path.exists(input_csv):
    print(f"错误：输入文件 {input_csv} 不存在")
    exit(1)

# 创建输出目录
output_path = os.path.join(output_dir, args.dataset, args.particle)
os.makedirs(output_path, exist_ok=True)

# 输出文件路径
output_csv = os.path.join(output_path, f"fit_obvs_{args.task}.csv")
img_dir = os.path.join(output_path, f"fitplots_obvs_{args.task}")
os.makedirs(img_dir, exist_ok=True)

print(f"处理数据集: {args.dataset}, 粒子类型: {args.particle}, 任务: {args.task}")
print(f"输入文件: {input_csv}")
print(f"输出目录: {output_path}")

# 线性模型
def sigbkg_model(fs, a, b):
    return a * fs + b * (1 - fs)

def calculate_weighted_average(inv_mass, inv_mass_center, obv, obv_err):
    """
    计算给定质量窗口内的加权平均值
    
    Args:
        inv_mass: 不变质量
        inv_mass_center: 不变质量中心值
        obv: 观测量
        obv_err: 观测量误差
        
    Returns:
        weighted_mean: 加权平均值
        weighted_error: 加权平均值的误差
        chi2: chi2值（此处为NaN）
        n_points: 窗口内的数据点数
    """
    center_mass = 1.115683  # Lambda粒子的静止质量
    range_mass = args.range_mass  # 窗口半宽度，表示中心质量左右各±该值
    
    # 在指定质量窗口范围内选择数据点
    mask = np.abs(inv_mass - center_mass) < range_mass  # 选择中心质量±range_mass范围内的点
    mask = mask & np.isfinite(obv) & np.isfinite(obv_err) & (obv_err > 0)
    
    if np.count_nonzero(mask) < 1:
        print(f"警告: 质量窗口 {center_mass-range_mass:.6f}-{center_mass+range_mass:.6f} 内没有有效数据点")
        return np.nan, np.nan, np.nan, 0
    
    selected_obv = obv[mask]
    selected_err = obv_err[mask]
    selected_mass = inv_mass_center[mask]
    
    # 计算加权平均值
    weights = 1.0 / (selected_err ** 2)
    weighted_mean = np.sum(selected_obv * weights) / np.sum(weights)
    weighted_error = np.sqrt(1.0 / np.sum(weights))
    
    print(f"在质量窗口 {center_mass-range_mass:.6f}-{center_mass+range_mass:.6f} 内找到 {len(selected_obv)} 个点")
    print(f"加权平均值: {weighted_mean:.6f} ± {weighted_error:.6f}")
    
    return weighted_mean, weighted_error, np.nan, len(selected_obv)

def fit_observable(fs, obv, obv_err, inv_mass=None, inv_mass_center=None):
    # 如果选择不拟合，直接计算加权平均值
    if args.no_fit and inv_mass is not None and inv_mass_center is not None:
        weighted_mean, weighted_error, _, n_points = calculate_weighted_average(
            inv_mass, inv_mass_center, obv, obv_err)
        # 返回信号值作为a，背景值设为NaN
        return weighted_mean, weighted_error, np.nan, np.nan, np.nan, n_points
    
    # 常规拟合流程
    mask = np.isfinite(fs) & np.isfinite(obv) & np.isfinite(obv_err) & (obv_err > 0)

    # 应用fs过滤（如果指定）
    if args.fs_min is not None:
        mask = mask & (fs >= args.fs_min)
    if args.fs_max is not None:
        mask = mask & (fs <= args.fs_max)

    if np.count_nonzero(mask) < 2:
        print(f"警告: 有效数据点不足(<2)，无法进行拟合")
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0

    fs_filtered = fs[mask]
    obv_filtered = obv[mask]
    obv_err_filtered = obv_err[mask]

    # 诊断数据特征
    fs_range = fs_filtered.max() - fs_filtered.min()
    if fs_range < 0.1:
        print(f"警告: fs范围较窄 ({fs_filtered.min():.3f} - {fs_filtered.max():.3f})，可能导致拟合不稳定")

    try:
        # 使用更强健的拟合参数
        popt, pcov = curve_fit(
            sigbkg_model,
            fs_filtered,
            obv_filtered,
            sigma=obv_err_filtered,
            absolute_sigma=True,
            method=args.fit_method,  # 使用命令行指定的方法
            max_nfev=10000  # 增加最大函数评估次数
        )
        a, b = popt

        # 检查协方差矩阵
        if np.any(np.diag(pcov) < 0):
            print(f"警告: 协方差矩阵对角线存在负值，使用替代误差估计")
            # 使用参数值的10%作为替代误差估计
            a_err = abs(0.1 * a) if a != 0 else 0.1
            b_err = abs(0.1 * b) if b != 0 else 0.1
        else:
            try:
                a_err, b_err = np.sqrt(np.diag(pcov))
            except Exception:
                print(f"警告: 无法从协方差矩阵计算参数误差，使用替代误差估计")
                # 使用参数值的10%作为替代误差估计
                a_err = abs(0.1 * a) if a != 0 else 0.1
                b_err = abs(0.1 * b) if b != 0 else 0.1

        residuals = (obv_filtered - sigbkg_model(fs_filtered, *popt)) / obv_err_filtered
        chi2 = np.sum(residuals ** 2)
        ndf = len(obv_filtered) - 2
        chi2_ndf = chi2 / ndf if ndf > 0 else np.nan

        # 检查拟合质量
        if chi2_ndf > 10:
            print(f"警告: 拟合质量较差 (χ²/ndf = {chi2_ndf:.2f})，参数可能不可靠")

        return a, a_err, b, b_err, chi2_ndf, len(obv_filtered)
    except Exception as e:
        print(f"拟合失败: {str(e)}")
        return np.nan, np.nan, np.nan, np.nan, np.nan, len(obv_filtered)

def plot_grid_inv_mass(centrality_list, group_list, a, b, label, sup_title, pdf):
    n = len(centrality_list)
    ncols = 3
    nrows = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 12), sharex=False, sharey=False)
    axes = axes.flatten()
    for i, (centrality, group) in enumerate(zip(centrality_list, group_list)):
        if i >= len(axes):
            break

        inv_mass_center = group['inv_mass_center'].values
        obv = group['obv'].values
        obv_err = group['obv_err'].values
        fs = group['fs_mass'].values

        # 应用fs过滤（如果指定）
        mask = np.ones_like(fs, dtype=bool)
        if args.fs_min is not None:
            mask = mask & (fs >= args.fs_min)
        if args.fs_max is not None:
            mask = mask & (fs <= args.fs_max)

        sort_idx = np.argsort(inv_mass_center[mask])
        x = inv_mass_center[mask][sort_idx]
        y = obv[mask][sort_idx]
        yerr = obv_err[mask][sort_idx]
        fs_sorted = fs[mask][sort_idx]

        # 每个子图都用自己对应的a/b
        if args.no_fit:
            # 使用加权平均值显示一条水平线
            yfit = np.ones_like(fs_sorted) * a[i]
        else:
            # 标准拟合曲线
            yfit = a[i] * fs_sorted + b[i] * (1 - fs_sorted)

        ax = axes[i]
        ax.errorbar(x, y, yerr=yerr, fmt='o', color='blue', label="Data")
        if args.no_fit:
            ax.plot(x, yfit, '-', color='red', label=f"Weighted Avg: {a[i]:.4f}")
        else:
            ax.plot(x, yfit, 's-', color='red', label=f"Fit: a={a[i]:.4f}, b={b[i]:.4f}")
        ax.set_title(f"centrality={centrality}")
        ax.set_xlabel("inv_mass_center (GeV/c$^2$)")
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(alpha=0.3)

    # 隐藏未用子图
    for j in range(i+1, nrows*ncols):
        if j < len(axes):
            axes[j].set_visible(False)

    title = f"{args.dataset}/{args.particle}: {sup_title}"
    if args.no_fit:
        mass_window = f"Mass Window: {1.115683-args.range_mass:.6f}-{1.115683+args.range_mass:.6f} GeV/c²"
        title += f" | {mass_window}"
    elif args.fs_min is not None or args.fs_max is not None:
        fs_range = f"fs∈[{args.fs_min if args.fs_min is not None else '-∞'}, {args.fs_max if args.fs_max is not None else '∞'}]"
        title += f" | {fs_range}"

    fig.suptitle(title, fontsize=16)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)

# 读取数据
print(f"读取输入文件: {input_csv}")
df = pd.read_csv(input_csv)

# 控制警告显示
warnings.filterwarnings('always', category=OptimizeWarning)  # 始终显示优化警告

group_cols = ['centrality', 'diff_type', 'diff_bin', 'pair_type']
out_records = []

# 显示处理模式和过滤信息
if args.no_fit:
    mass_window = f"{1.115683-args.range_mass:.6f}-{1.115683+args.range_mass:.6f} GeV/c²"
    print(f"使用加权平均值模式，质量窗口范围: {mass_window} (中心质量左右各±{args.range_mass} GeV/c²)")
elif args.fs_min is not None or args.fs_max is not None:
    fs_range = f"fs∈[{args.fs_min if args.fs_min is not None else '-∞'}, {args.fs_max if args.fs_max is not None else '∞'}]"
    print(f"使用拟合模式，应用信号分数过滤: {fs_range}")
else:
    print(f"使用拟合模式，不应用信号分数过滤")

# 统计计数
fit_count = 0
fail_count = 0
pdf_count = 0

# 对于每一组(diff_type, diff_bin, pair_type)分开分页输出delta/rawgamma 3x3 pdf
group_main = df.groupby(['diff_type', 'diff_bin', 'pair_type'])
print(f"找到 {len(group_main)} 个分组")

for (diff_type, diff_bin, pair_type), main_group in group_main:
    # 为输出pdf准备所有centrality的子分组
    centrality_list = sorted(main_group['centrality'].unique())
    # 先拟合好每个centrality的delta/rawgamma参数，记录（主表只保留分组参数，不在分页图上重复拟合）
    delta_params = {}
    rawgamma_params = {}

    # 对所有centrality逐个做拟合并记录参数
    print(f"处理 diff_type={diff_type}, diff_bin={diff_bin}, pair_type={pair_type} (centralities={len(centrality_list)})")

    for centrality in centrality_list:
        sub = main_group[main_group['centrality'] == centrality]
        fs = sub['fs_mass'].values
        inv_mass_counts = sub['inv_mass_counts'].values
        delta = sub['delta'].values
        delta_err = sub['delta_err'].values
        rawg = sub['rawgamma'].values
        rawg_err = sub['rawgamma_err'].values

        # 获取不变质量数据（如果使用加权平均方法需要）
        inv_mass = sub['inv_mass'].values if 'inv_mass' in sub.columns else None
        inv_mass_center = sub['inv_mass_center'].values
        
        # 拟合delta或计算加权平均
        a_d, aerr_d, b_d, berr_d, chi2_d, N_d = fit_observable(
            fs, delta, delta_err, inv_mass, inv_mass_center)
            
        # 拟合rawgamma或计算加权平均
        a_g, aerr_g, b_g, berr_g, chi2_g, N_g = fit_observable(
            fs, rawg, rawg_err, inv_mass, inv_mass_center)

        # 计算信号和背景计数
        # 应用fs过滤（如果指定）
        mask = np.ones_like(fs, dtype=bool)
        if args.fs_min is not None:
            mask = mask & (fs >= args.fs_min)
        if args.fs_max is not None:
            mask = mask & (fs <= args.fs_max)

        fs_filtered = fs[mask]
        counts_filtered = inv_mass_counts[mask]

        counts_sig = np.sum(fs_filtered * counts_filtered)
        counts_bkg = np.sum((1 - fs_filtered) * counts_filtered)

        # 记录参数
        out_records.append(dict(
            centrality=centrality,
            diff_type=diff_type,
            diff_bin=diff_bin,
            pair_type=pair_type,
            delta_sig=a_d,
            delta_err_sig=aerr_d,
            delta_bkg=b_d,
            delta_err_bkg=berr_d,
            rawgamma_sig=a_g,
            rawgamma_err_sig=aerr_g,
            rawgamma_bkg=b_g,
            rawgamma_err_bkg=berr_g,
            counts_sig=counts_sig,
            counts_bkg=counts_bkg,
            chi2_ndf_delta=chi2_d,
            chi2_ndf_gamma=chi2_g,
            N_points_delta=N_d,
            N_points_gamma=N_g
        ))

        if not np.isnan(a_d) and not np.isnan(a_g):
            fit_count += 1
        else:
            fail_count += 1

        delta_params[centrality] = (a_d, b_d)
        rawgamma_params[centrality] = (a_g, b_g)

    # --- 分页输出delta 3x3 pdf ---
    pdfname_delta = os.path.join(img_dir, f"delta_vs_inv_mass_grid_{diff_type}_{diff_bin}_{pair_type}.pdf")
    with PdfPages(pdfname_delta) as pdf:
        for i in range(0, len(centrality_list), 9):
            centers = centrality_list[i:i+9]
            sub_groups = []
            for c in centers:
                sub = main_group[main_group['centrality'] == c].copy()
                # 用于画图统一接口
                sub = sub.rename(columns={'delta': 'obv', 'delta_err': 'obv_err'})
                sub_groups.append(sub)
            # 取每个centrality对应的参数
            a_b_list = [delta_params[c] for c in centers]
            # 按9个一页画
            plot_grid_inv_mass(
                centers, sub_groups,
                a=np.array([ab[0] for ab in a_b_list]),
                b=np.array([ab[1] for ab in a_b_list]),
                label=r"$\delta$",
                sup_title=f"delta | diff_type={diff_type}, diff_bin={diff_bin}, pair_type={pair_type}",
                pdf=pdf
            )
    pdf_count += 1

    # --- 分页输出rawgamma 3x3 pdf ---
    pdfname_gamma = os.path.join(img_dir, f"rawgamma_vs_inv_mass_grid_{diff_type}_{diff_bin}_{pair_type}.pdf")
    with PdfPages(pdfname_gamma) as pdf:
        for i in range(0, len(centrality_list), 9):
            centers = centrality_list[i:i+9]
            sub_groups = []
            for c in centers:
                sub = main_group[main_group['centrality'] == c].copy()
                sub = sub.rename(columns={'rawgamma': 'obv', 'rawgamma_err': 'obv_err'})
                sub_groups.append(sub)
            a_b_list = [rawgamma_params[c] for c in centers]
            plot_grid_inv_mass(
                centers, sub_groups,
                a=np.array([ab[0] for ab in a_b_list]),
                b=np.array([ab[1] for ab in a_b_list]),
                label=r"$\gamma$",
                sup_title=f"rawgamma | diff_type={diff_type}, diff_bin={diff_bin}, pair_type={pair_type}",
                pdf=pdf
            )
    pdf_count += 1
    print(f"已输出PDF: delta和rawgamma图像 - {diff_type}_{diff_bin}_{pair_type}")

# --- 输出csv ---
pd.DataFrame(out_records).to_csv(output_csv, index=False)
print(f"拟合成功: {fit_count}, 拟合失败: {fail_count}")
print(f"生成 {pdf_count} 个PDF文件")
print(f"所有分组拟合完成，结果已写入 {output_csv}，所有图片保存在 {img_dir}/")
