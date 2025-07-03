import pandas as pd
import numpy as np
from scipy.optimize import curve_fit, minimize
import ROOT
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='观测量模板拟合脚本')
parser.add_argument('-i', '--input-dir', type=str, default='./',
                   help='输入目录')
parser.add_argument('-d', '--dataset', type=str, default=None,
                   help='数据集名称（可选，自动推断）')
parser.add_argument('-t', '--task', type=str, default='default',
                   help='任务名称')
parser.add_argument('-o', '--output-dir', type=str, default='./',
                   help='输出目录')
parser.add_argument('--fit-mode', choices=['2Dfit', '1Dfit', 'narrMass', 'all'],
                   default='2Dfit', help='拟合方法选择 (默认: 2Dfit)')
parser.add_argument('--range_mass', type=float, default=0.0014, help='质量窗口范围 (默认: 0.0014)')
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

print(f"输入目录: {input_path}")
print(f"输出目录: {output_path}")
print(f"数据集: {args.dataset}")
print(f"任务: {args.task}")
print(f"拟合模式: {args.fit_mode}")
print(f"质量窗口范围: {args.range_mass}")

# Lambda粒子质量
LAMBDA_MASS = 1.115683

# 输入输出文件
input_csv = os.path.join(input_path, f"fit_inv_mass_{args.task}.csv")
output_csv = os.path.join(output_path, f"fit_obvs_{args.task}.csv")
output_root = os.path.join(output_path, f"fit_obvs_plots_{args.task}.root")

print(f"读取输入文件: {input_csv}")

# 读取数据
df = pd.read_csv(input_csv)
print(f"数据形状: {df.shape}")
print(f"列名: {df.columns.tolist()}")

# 检查必要的列是否存在
required_cols = ['centrality', 'pair_type', 'inv_mass_1_center', 'inv_mass_2_center',
                'delta', 'delta_err', 'rawgamma', 'rawgamma_err', 'fs_mass_1', 'fs_mass_2']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    raise ValueError(f"缺少必要的列: {missing_cols}")

print("所有必要的列都存在")

# 模板拟合函数
def template_fit_function(fs_data, obv_ss, obv_sb, obv_bb):
    """
    模板拟合函数
    fs_data: (fs1, fs2) 的数组，形状为 (N, 2)
    obv_ss, obv_sb, obv_bb: 拟合参数
    """
    fs1, fs2 = fs_data[:, 0], fs_data[:, 1]
    return obv_ss * fs1 * fs2 + obv_sb * (fs1 * (1 - fs2) + (1 - fs1) * fs2) + obv_bb * (1 - fs1) * (1 - fs2)

def fit_observable_2D(df_group, observable_col, error_col):
    """
    2D模板拟合方法 (原始方法)
    """
    # 准备数据
    fs1 = df_group['fs_mass_1'].values
    fs2 = df_group['fs_mass_2'].values
    obv_data = df_group[observable_col].values
    obv_err = df_group[error_col].values

    # 过滤掉无效数据
    valid_mask = (~np.isnan(obv_data)) & (~np.isnan(obv_err)) & (obv_err > 0) & (~np.isnan(fs1)) & (~np.isnan(fs2))

    if np.sum(valid_mask) < 3:
        raise ValueError(f"有效数据点太少: {np.sum(valid_mask)}")

    fs1_valid = fs1[valid_mask]
    fs2_valid = fs2[valid_mask]
    obv_data_valid = obv_data[valid_mask]
    obv_err_valid = obv_err[valid_mask]

    # 构建fs数据数组
    fs_data = np.column_stack([fs1_valid, fs2_valid])

    # 定义拟合函数（适配curve_fit的格式）
    def fit_func(dummy, obv_ss, obv_sb, obv_bb):
        return template_fit_function(fs_data, obv_ss, obv_sb, obv_bb)

    # 初始猜测值
    initial_guess = [np.mean(obv_data_valid), 0.0, 0.0]

    try:
        # 使用curve_fit进行拟合
        popt, pcov = curve_fit(fit_func, np.arange(len(obv_data_valid)), obv_data_valid,
                              sigma=obv_err_valid, p0=initial_guess, absolute_sigma=True)

        # 计算参数误差
        param_errors = np.sqrt(np.diag(pcov))

        # 计算拟合结果
        obv_fit = template_fit_function(fs_data, *popt)

        # 计算卡方
        chi2 = np.sum(((obv_data_valid - obv_fit) / obv_err_valid) ** 2)
        ndf = len(obv_data_valid) - 3

        return {
            'success': True,
            'ss': popt[0], 'ss_err': param_errors[0],
            'sb': popt[1], 'sb_err': param_errors[1],
            'bb': popt[2], 'bb_err': param_errors[2],
            'chi2': chi2, 'ndf': ndf,
            'fit_data': obv_fit,
            'original_data': obv_data_valid,
            'fs_data': fs_data,
            'inv_mass_1': df_group['inv_mass_1_center'].values[valid_mask],
            'inv_mass_2': df_group['inv_mass_2_center'].values[valid_mask],
            'valid_mask': valid_mask
        }

    except Exception as e:
        raise RuntimeError(f"2D拟合失败: {str(e)}")

def fit_observable_1D(df_group, observable_col, error_col, range_mass):
    """
    1D拟合方法：卡紧mass_2，使用 obv_data = obv_s * fs + obv_b * (1-fs)
    """
    # 卡紧mass_2在Lambda质量附近
    mass2_mask = np.abs(df_group['inv_mass_2_center'] - LAMBDA_MASS) <= range_mass
    df_filtered = df_group[mass2_mask]

    if len(df_filtered) < 3:
        raise ValueError(f"卡紧mass_2后数据点太少: {len(df_filtered)}")

    # 准备数据
    fs1 = df_filtered['fs_mass_1'].values
    obv_data = df_filtered[observable_col].values
    obv_err = df_filtered[error_col].values

    # 过滤掉无效数据
    valid_mask = (~np.isnan(obv_data)) & (~np.isnan(obv_err)) & (obv_err > 0) & (~np.isnan(fs1))

    if np.sum(valid_mask) < 2:
        raise ValueError(f"有效数据点太少: {np.sum(valid_mask)}")

    fs1_valid = fs1[valid_mask]
    obv_data_valid = obv_data[valid_mask]
    obv_err_valid = obv_err[valid_mask]

    # 定义1D拟合函数
    def fit_func_1d(fs, obv_s, obv_b):
        return obv_s * fs + obv_b * (1 - fs)

    # 初始猜测值
    initial_guess = [np.mean(obv_data_valid), 0.0]

    try:
        # 使用curve_fit进行拟合，只使用误差权重
        popt, pcov = curve_fit(fit_func_1d, fs1_valid, obv_data_valid,
                              sigma=obv_err_valid, p0=initial_guess, absolute_sigma=True)

        # 计算参数误差
        param_errors = np.sqrt(np.diag(pcov))

        # 计算拟合结果
        obv_fit = fit_func_1d(fs1_valid, *popt)

        # 计算卡方
        chi2 = np.sum(((obv_data_valid - obv_fit) / obv_err_valid) ** 2)
        ndf = len(obv_data_valid) - 2

        return {
            'success': True,
            'ss': popt[0], 'ss_err': param_errors[0],
            'sb': popt[1], 'sb_err': param_errors[1],
            'bb': 0.0, 'bb_err': 0.0,  # 1D方法中没有bb项
            'chi2': chi2, 'ndf': ndf,
            'fit_data': obv_fit,
            'original_data': obv_data_valid,
            'fs_data': fs1_valid,
            'inv_mass_1': df_filtered['inv_mass_1_center'].values[valid_mask],
            'inv_mass_2': df_filtered['inv_mass_2_center'].values[valid_mask],
            'valid_mask': valid_mask,
            'method': '1D'
        }

    except Exception as e:
        raise RuntimeError(f"1D拟合失败: {str(e)}")

def fit_observable_narrow(df_group, observable_col, error_col, range_mass):
    """
    窄质量窗口方法：取mass_1和mass_2都在Lambda质量附近的点的加权平均
    """
    # 卡紧mass_1和mass_2都在Lambda质量附近
    mass1_mask = np.abs(df_group['inv_mass_1_center'] - LAMBDA_MASS) <= range_mass
    mass2_mask = np.abs(df_group['inv_mass_2_center'] - LAMBDA_MASS) <= range_mass
    combined_mask = mass1_mask & mass2_mask

    df_filtered = df_group[combined_mask]

    if len(df_filtered) == 0:
        raise ValueError("窄质量窗口内没有数据点")

    # 准备数据
    obv_data = df_filtered[observable_col].values
    obv_err = df_filtered[error_col].values
    inv_mass_counts = df_filtered['inv_mass_counts'].values

    # 过滤掉无效数据
    valid_mask = (~np.isnan(obv_data)) & (~np.isnan(obv_err)) & (obv_err > 0) & (~np.isnan(inv_mass_counts)) & (inv_mass_counts > 0)

    if np.sum(valid_mask) == 0:
        raise ValueError("窄质量窗口内没有有效数据点")

    obv_data_valid = obv_data[valid_mask]
    obv_err_valid = obv_err[valid_mask]
    inv_mass_counts_valid = inv_mass_counts[valid_mask]

    # 计算加权平均，使用inv_mass_counts作为权重因子
    weights = inv_mass_counts_valid / (obv_err_valid ** 2)
    weighted_mean = np.sum(obv_data_valid * weights) / np.sum(weights)
    # 正确的加权平均误差公式
    weighted_err = 1.0 / np.sqrt(np.sum(inv_mass_counts_valid / (obv_err_valid ** 2)))

    return {
        'success': True,
        'ss': weighted_mean, 'ss_err': weighted_err,
        'sb': 0.0, 'sb_err': 0.0,  # 窄质量方法中没有sb和bb项
        'bb': 0.0, 'bb_err': 0.0,
        'chi2': 0.0, 'ndf': len(obv_data_valid) - 1,
        'fit_data': np.full_like(obv_data_valid, weighted_mean),
        'original_data': obv_data_valid,
        'fs_data': None,
        'inv_mass_1': df_filtered['inv_mass_1_center'].values[valid_mask],
        'inv_mass_2': df_filtered['inv_mass_2_center'].values[valid_mask],
        'valid_mask': valid_mask,
        'method': 'narrow',
        'n_points': len(obv_data_valid)
    }

def fit_observable(df_group, observable_col, error_col, mode, range_mass):
    """
    根据模式选择拟合方法
    """
    if mode == '2Dfit':
        return fit_observable_2D(df_group, observable_col, error_col)
    elif mode == '1Dfit':
        return fit_observable_1D(df_group, observable_col, error_col, range_mass)
    elif mode == 'narrMass':
        return fit_observable_narrow(df_group, observable_col, error_col, range_mass)
    else:
        raise ValueError(f"未知的拟合模式: {mode}")

# 按centrality和pair_type分组
group_cols = ['centrality', 'pair_type']
all_groups = df.groupby(group_cols)

# 如果是all模式，需要对每种方法都运行
if args.fit_mode == 'all':
    modes_to_run = ['2Dfit', '1Dfit', 'narrMass']
else:
    modes_to_run = [args.fit_mode]

for current_mode in modes_to_run:
    print(f"\n=== 运行模式: {current_mode} ===")

    # 重新设置输出文件名
    if args.fit_mode == 'all':
        current_output_csv = os.path.join(output_path, f"fit_obvs_{current_mode}_{args.task}.csv")
        current_output_root = os.path.join(output_path, f"fit_obvs_plots_{current_mode}_{args.task}.root")
    else:
        current_output_csv = output_csv
        current_output_root = output_root

    # 新增：提前创建ratio root文件（每个模式一个）
    if args.fit_mode == 'all':
        current_ratio_root = os.path.join(output_path, f"fit_obvs_plots_ratio_{current_mode}_{args.task}.root")
    else:
        current_ratio_root = os.path.join(output_path, f"fit_obvs_plots_ratio_{args.task}.root")
    ratio_root_file = ROOT.TFile(current_ratio_root, "RECREATE")

    # 创建ROOT文件用于保存图形
    current_root_file = ROOT.TFile(current_output_root, "RECREATE")
    print(f"创建ROOT文件: {current_output_root}")

    current_fit_results = []

    for group_keys, group_df in all_groups:
        centrality, pair_type = group_keys
        print(f"\n处理组合: centrality={centrality}, pair_type={pair_type}")
        print(f"数据点数量: {len(group_df)}")

        try:
            # 拟合delta
            print("拟合delta...")
            delta_result = fit_observable(group_df, 'delta', 'delta_err', current_mode, args.range_mass)
            print(f"delta拟合成功: ss={delta_result['ss']:.6f}±{delta_result['ss_err']:.6f}")

            # 拟合rawgamma
            print("拟合rawgamma...")
            gamma_result = fit_observable(group_df, 'rawgamma', 'rawgamma_err', current_mode, args.range_mass)
            print(f"rawgamma拟合成功: ss={gamma_result['ss']:.6f}±{gamma_result['ss_err']:.6f}")

            # 保存拟合结果
            result_dict = {
                'centrality': centrality,
                'pair_type': pair_type,
                'inv_mass_counts': group_df['inv_mass_counts'].sum(),
                'delta_ss': delta_result['ss'], 'delta_ss_err': delta_result['ss_err'],
                'delta_sb': delta_result['sb'], 'delta_sb_err': delta_result['sb_err'],
                'delta_bb': delta_result['bb'], 'delta_bb_err': delta_result['bb_err'],
                'rawgamma_ss': gamma_result['ss'], 'rawgamma_ss_err': gamma_result['ss_err'],
                'rawgamma_sb': gamma_result['sb'], 'rawgamma_sb_err': gamma_result['sb_err'],
                'rawgamma_bb': gamma_result['bb'], 'rawgamma_bb_err': gamma_result['bb_err'],
                'delta_chi2': delta_result['chi2'], 'delta_ndf': delta_result['ndf'],
                'gamma_chi2': gamma_result['chi2'], 'gamma_ndf': gamma_result['ndf'],
                'fit_mode': current_mode
            }

            # 为窄质量方法添加额外信息
            if current_mode == 'narrMass':
                result_dict['delta_n_points'] = delta_result.get('n_points', 0)
                result_dict['gamma_n_points'] = gamma_result.get('n_points', 0)

            current_fit_results.append(result_dict)

            # 创建可视化图形
            canvas_name = f"canvas_{pair_type}_cent_{centrality}_{current_mode}"
            canvas = ROOT.TCanvas(canvas_name, f"{pair_type}, centrality={centrality}, mode={current_mode}", 600, 600)
            canvas.Divide(2, 2)

            # 获取质量范围
            inv_mass_1_vals = delta_result['inv_mass_1']
            inv_mass_2_vals = delta_result['inv_mass_2']

            mass1_min, mass1_max = inv_mass_1_vals.min(), inv_mass_1_vals.max()
            mass2_min, mass2_max = inv_mass_2_vals.min(), inv_mass_2_vals.max()

            # 创建直方图
            nbins = 20

            # Delta原始数据
            canvas.cd(1)
            hist_delta_orig = ROOT.TH2D(f"delta_orig_{pair_type}_{centrality}_{current_mode}",
                                       f"Delta Original;inv_mass_1;inv_mass_2",
                                       nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)

            # Delta拟合结果
            canvas.cd(2)
            hist_delta_fit = ROOT.TH2D(f"delta_fit_{pair_type}_{centrality}_{current_mode}",
                                      f"Delta Fitted;inv_mass_1;inv_mass_2",
                                      nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)

            # Rawgamma原始数据
            canvas.cd(3)
            hist_gamma_orig = ROOT.TH2D(f"gamma_orig_{pair_type}_{centrality}_{current_mode}",
                                       f"Rawgamma Original;inv_mass_1;inv_mass_2",
                                       nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)

            # Rawgamma拟合结果
            canvas.cd(4)
            hist_gamma_fit = ROOT.TH2D(f"gamma_fit_{pair_type}_{centrality}_{current_mode}",
                                      f"Rawgamma Fitted;inv_mass_1;inv_mass_2",
                                      nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)

            # 填充直方图
            for i in range(len(inv_mass_1_vals)):
                hist_delta_orig.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], delta_result['original_data'][i])
                hist_delta_fit.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], delta_result['fit_data'][i])
                hist_gamma_orig.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], gamma_result['original_data'][i])
                hist_gamma_fit.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], gamma_result['fit_data'][i])

            # 绘制直方图，先画fit前的，获取z范围，再画fit后的
            canvas.cd(1)
            hist_delta_orig.Draw("COLZ")
            zmin_delta = hist_delta_orig.GetMinimum()
            zmax_delta = hist_delta_orig.GetMaximum()

            canvas.cd(2)
            hist_delta_fit.SetMinimum(zmin_delta)
            hist_delta_fit.SetMaximum(zmax_delta)
            hist_delta_fit.Draw("COLZ")

            canvas.cd(3)
            hist_gamma_orig.Draw("COLZ")
            zmin_gamma = hist_gamma_orig.GetMinimum()
            zmax_gamma = hist_gamma_orig.GetMaximum()

            canvas.cd(4)
            hist_gamma_fit.SetMinimum(zmin_gamma)
            hist_gamma_fit.SetMaximum(zmax_gamma)
            hist_gamma_fit.Draw("COLZ")

            # 保存到ROOT文件
            current_root_file.cd()
            canvas.Write(canvas.GetName())
            # 保存主图为pdf
            plot_fit_dir = os.path.join(output_path, "plot_fit")
            os.makedirs(plot_fit_dir, exist_ok=True)
            pdf_path = os.path.join(plot_fit_dir, f"{canvas.GetName()}.pdf")
            canvas.SaveAs(pdf_path)

            # 新增：ratio图（delta和gamma）
            # delta ratio
            hist_delta_ratio = ROOT.TH2D(f"delta_ratio_{pair_type}_{centrality}_{current_mode}",
                                         f"Delta Fit/Original Ratio;inv_mass_1;inv_mass_2",
                                         nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)
            for i in range(len(inv_mass_1_vals)):
                orig = delta_result['original_data'][i]
                fit = delta_result['fit_data'][i]
                ratio = fit / orig if orig != 0 else 0
                hist_delta_ratio.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], ratio)
            canvas_delta_ratio = ROOT.TCanvas(f"canvas_delta_ratio_{pair_type}_cent_{centrality}_{current_mode}",
                                             f"Delta Ratio, {pair_type}, centrality={centrality}, mode={current_mode}", 600, 600)
            hist_delta_ratio.Draw("COLZ")
            ratio_root_file.cd()
            canvas_delta_ratio.Write(canvas_delta_ratio.GetName())
            # 保存delta ratio为pdf
            pdf_path_delta = os.path.join(plot_fit_dir, f"{canvas_delta_ratio.GetName()}.pdf")
            canvas_delta_ratio.SaveAs(pdf_path_delta)

            # gamma ratio
            hist_gamma_ratio = ROOT.TH2D(f"gamma_ratio_{pair_type}_{centrality}_{current_mode}",
                                         f"Rawgamma Fit/Original Ratio;inv_mass_1;inv_mass_2",
                                         nbins, mass1_min, mass1_max, nbins, mass2_min, mass2_max)
            for i in range(len(inv_mass_1_vals)):
                orig = gamma_result['original_data'][i]
                fit = gamma_result['fit_data'][i]
                ratio = fit / orig if orig != 0 else 0
                hist_gamma_ratio.Fill(inv_mass_1_vals[i], inv_mass_2_vals[i], ratio)
            canvas_gamma_ratio = ROOT.TCanvas(f"canvas_gamma_ratio_{pair_type}_cent_{centrality}_{current_mode}",
                                             f"Rawgamma Ratio, {pair_type}, centrality={centrality}, mode={current_mode}", 600, 600)
            hist_gamma_ratio.Draw("COLZ")
            ratio_root_file.cd()
            canvas_gamma_ratio.Write(canvas_gamma_ratio.GetName())
            # 保存gamma ratio为pdf
            pdf_path_gamma = os.path.join(plot_fit_dir, f"{canvas_gamma_ratio.GetName()}.pdf")
            canvas_gamma_ratio.SaveAs(pdf_path_gamma)

            # 清理内存
            canvas.Delete()
            hist_delta_orig.Delete()
            hist_delta_fit.Delete()
            hist_gamma_orig.Delete()
            hist_gamma_fit.Delete()

            print(f"已保存图形: {canvas_name}")

        except Exception as e:
            print(f"处理失败: {str(e)}")
            continue

    # 关闭ROOT文件
    current_root_file.Close()
    # 新增：关闭ratio root文件
    ratio_root_file.Close()

    # 保存拟合结果到CSV
    if current_fit_results:
        current_fit_df = pd.DataFrame(current_fit_results)
        current_fit_df.to_csv(current_output_csv, index=False)
        print(f"\n{current_mode}模式拟合完成！结果已保存到: {current_output_csv}")
        print(f"成功拟合的组合数量: {len(current_fit_results)}")
        print(f"可视化图形已保存到: {current_output_root}")

        print(f"\n{current_mode}模式拟合结果预览:")
        print(current_fit_df[['centrality', 'pair_type', 'delta_ss', 'delta_ss_err', 'rawgamma_ss', 'rawgamma_ss_err']].head())
    else:
        print(f"{current_mode}模式没有成功的拟合结果！")

print(f"\n所有拟合模式完成！")
