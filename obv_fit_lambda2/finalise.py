import pandas as pd
import numpy as np
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='汇总Lambda粒子对观测量结果并应用分辨率校正')
parser.add_argument('-i', '--input-dir', type=str, default='./',
                   help='输入目录')
parser.add_argument('-d', '--dataset', type=str, default=None,
                   help='数据集名称（可选，自动推断）')
parser.add_argument('-t', '--task', type=str, default='default',
                   help='任务名称')
parser.add_argument('-o', '--output-dir', type=str, default='./',
                   help='输出目录')
parser.add_argument('-r', '--reso-file', type=str, default="./resolutions.csv",
                   help='分辨率文件路径')
parser.add_argument('--period', type=str, default=None,
                   help='期间标识（如18q，18r），可选，自动推断')
parser.add_argument('--plane', type=str, default="TPC",
                   help='平面标识，默认为TPC')

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

# 如果没有指定period，从数据集推断
if args.period is None:
    if '18q' in args.dataset:
        args.period = '18q'
    elif '18r' in args.dataset:
        args.period = '18r'
    else:
        args.period = '18q'  # 默认值

# 构建输入输出路径
input_path = os.path.join(args.input_dir, args.dataset, particle_type)
output_path = os.path.join(args.output_dir, args.dataset, particle_type)
os.makedirs(output_path, exist_ok=True)

# 构建文件路径
input_csv = os.path.join(input_path, f"fit_obvs_{args.task}.csv")
output_csv = os.path.join(output_path, f"finalise_{args.task}.csv")

print(f"输入目录: {input_path}")
print(f"输出目录: {output_path}")
print(f"数据集: {args.dataset}")
print(f"任务: {args.task}")
print(f"期间: {args.period}")

# 检查输入文件是否存在
if not os.path.exists(input_csv):
    print(f"错误：输入文件 {input_csv} 不存在")
    exit(1)

# 检查分辨率文件是否存在
if not os.path.exists(args.reso_file):
    print(f"错误：分辨率文件 {args.reso_file} 不存在")
    exit(1)

# 读取数据
print(f"读取输入文件: {input_csv}")
main = pd.read_csv(input_csv)

print(f"读取分辨率文件: {args.reso_file}")
reso = pd.read_csv(args.reso_file)

# 只选指定period/plane的resolution
print(f"筛选分辨率数据: period={args.period}, plane={args.plane}")
reso_sel = reso[(reso['period']==args.period) & (reso['plane']==args.plane)][['centrality', 'resolution']]
reso_sel = reso_sel.drop_duplicates('centrality')

if len(reso_sel) == 0:
    print(f"警告：在分辨率文件中找不到匹配的数据 (period={args.period}, plane={args.plane})")

# Lambda粒子对的SS/OS分组
# SS (Same Sign): Lambda-Lambda, LambdaBar-LambdaBar
# OS (Opposite Sign): Lambda-LambdaBar, LambdaBar-Lambda
SS_types = ["LambdaLambda", "LambdaBarLambdaBar"]
OS_types = ["LambdaLambdaBar", "LambdaBarLambda"]

print(f"SS类型 (Same Sign): {SS_types}")
print(f"OS类型 (Opposite Sign): {OS_types}")

# 标准误差加权平均函数
def error_weighted_mean(vals, errs):
    """
    使用标准的误差加权平均公式:
    μ̂ = (μ₁/σ₁² + μ₂/σ₂²) / (1/σ₁² + 1/σ₂²)
    σ_μ̂ = (1/σ₁² + 1/σ₂²)^(-1/2)
    """
    if len(vals) == 0:
        return float('nan'), float('nan')

    valid_mask = ~pd.isna(vals) & ~pd.isna(errs) & (errs > 0)
    if not any(valid_mask):
        return float('nan'), float('nan')

    vals_valid = vals[valid_mask]
    errs_valid = errs[valid_mask]

    # 权重 = 1/σ²
    weights = 1.0 / (errs_valid ** 2)

    # 加权平均
    weighted_mean = np.sum(weights * vals_valid) / np.sum(weights)

    # 加权平均的误差
    weighted_err = 1.0 / np.sqrt(np.sum(weights))

    return weighted_mean, weighted_err

# 汇总新表
records = []
success_count = 0
failed_reso = 0

for centrality, group in main.groupby('centrality'):
    # 选SS、OS
    SS = group[group['pair_type'].isin(SS_types)]
    OS = group[group['pair_type'].isin(OS_types)]

    # 检查是否有匹配的数据
    if len(SS) == 0:
        print(f"警告：在 centrality={centrality} 中找不到SS类型数据")
    if len(OS) == 0:
        print(f"警告：在 centrality={centrality} 中找不到OS类型数据")

    # SS误差加权平均
    delta_ss, delta_err_ss = error_weighted_mean(
        SS['delta_ss'].values if len(SS) > 0 else np.array([]),
        SS['delta_ss_err'].values if len(SS) > 0 else np.array([])
    )
    gamma_ss, gamma_err_ss = error_weighted_mean(
        SS['rawgamma_ss'].values if len(SS) > 0 else np.array([]),
        SS['rawgamma_ss_err'].values if len(SS) > 0 else np.array([])
    )

    # OS误差加权平均
    delta_os, delta_err_os = error_weighted_mean(
        OS['delta_ss'].values if len(OS) > 0 else np.array([]),
        OS['delta_ss_err'].values if len(OS) > 0 else np.array([])
    )
    gamma_os, gamma_err_os = error_weighted_mean(
        OS['rawgamma_ss'].values if len(OS) > 0 else np.array([]),
        OS['rawgamma_ss_err'].values if len(OS) > 0 else np.array([])
    )

    # Del = OS - SS (异号减同号)
    delta_del = delta_os - delta_ss
    delta_err_del = np.sqrt(delta_err_os**2 + delta_err_ss**2) if not pd.isna(delta_err_os) and not pd.isna(delta_err_ss) else float('nan')
    gamma_del = gamma_os - gamma_ss
    gamma_err_del = np.sqrt(gamma_err_os**2 + gamma_err_ss**2) if not pd.isna(gamma_err_os) and not pd.isna(gamma_err_ss) else float('nan')

    # 查找resolution（要求centrality完全匹配）
    try:
        resolution = float(reso_sel[reso_sel['centrality'] == centrality]['resolution'].iloc[0])
        success_count += 1
    except Exception as e:
        print(f"警告：找不到centrality={centrality}的分辨率值: {str(e)}")
        resolution = float('nan')
        failed_reso += 1

    # 归一化gamma（除以分辨率）
    gamma_ss_n = gamma_ss / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')
    gamma_err_ss_n = gamma_err_ss / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')
    gamma_os_n = gamma_os / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')
    gamma_err_os_n = gamma_err_os / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')
    gamma_del_n = gamma_del / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')
    gamma_err_del_n = gamma_err_del / resolution if resolution and not pd.isna(resolution) and resolution != 0 else float('nan')

    # 输出三类: SS, OS, Del
    for ptype, d, de, g, ge in [
        ("SS", delta_ss, delta_err_ss, gamma_ss_n, gamma_err_ss_n),
        ("OS", delta_os, delta_err_os, gamma_os_n, gamma_err_os_n),
        ("Del", delta_del, delta_err_del, gamma_del_n, gamma_err_del_n)
    ]:
        records.append(dict(
            centrality=centrality,
            pair_type=ptype,
            delta=d,
            delta_err=de,
            gamma=g,
            gamma_err=ge,
            resolution=resolution
        ))

# 输出结果
result_df = pd.DataFrame(records)
result_df.to_csv(output_csv, index=False)

print(f"分辨率匹配成功: {success_count}, 失败: {failed_reso}")
print(f"已输出: {output_csv}")
print(f"处理完成！共处理了 {len(result_df)} 条记录")
print(f"使用标准误差加权平均方法")

# 显示结果统计
print("\n结果统计:")
for ptype in ["SS", "OS", "Del"]:
    subset = result_df[result_df['pair_type'] == ptype]
    print(f"{ptype}: {len(subset)} 个中心度点")
    if len(subset) > 0:
        print(f"  Delta平均值: {subset['delta'].mean():.6f}")
        print(f"  Gamma平均值: {subset['gamma'].mean():.6f}")
