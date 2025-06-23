import pandas as pd
import numpy as np
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='汇总观测量结果并应用分辨率校正')
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
parser.add_argument('-r', '--reso-file', type=str, default="./resolutions.csv",
                   help='分辨率文件路径')

parser.add_argument('--plane', type=str, default="TPC",
                   help='平面标识，默认为TPC')
args = parser.parse_args()

# 设置输出目录
output_dir = args.output_dir if args.output_dir is not None else args.input_dir

# 从数据集中提取期间标识
if "18q" in args.dataset:
    period = "18q"
elif "18r" in args.dataset:
    period = "18r"
else:
    print(f"警告：无法从数据集 {args.dataset} 中提取期间标识，使用'18q'作为默认值")
    period = "18q"

print(f"使用期间标识: {period}, 平面标识: {args.plane}")

# 构建输入输出路径
input_path = os.path.join(args.input_dir, args.dataset, args.particle)
input_csv = os.path.join(input_path, f"fit_obvs_{args.task}.csv")

# 检查输入文件是否存在
if not os.path.exists(input_csv):
    print(f"错误：输入文件 {input_csv} 不存在")
    exit(1)

# 检查分辨率文件是否存在
if not os.path.exists(args.reso_file):
    print(f"错误：分辨率文件 {args.reso_file} 不存在")
    exit(1)

# 创建输出目录
output_path = os.path.join(output_dir, args.dataset, args.particle)
os.makedirs(output_path, exist_ok=True)

# 输出文件路径
output_csv = os.path.join(output_path, f"finalise_{args.task}.csv")

# 读取数据
print(f"读取输入文件: {input_csv}")
main = pd.read_csv(input_csv)

print(f"读取分辨率文件: {args.reso_file}")
reso = pd.read_csv(args.reso_file)

# 只选指定period/plane的resolution
print(f"筛选分辨率数据: period={period}, plane={args.plane}")
reso_sel = reso[(reso['period']==period) & (reso['plane']==args.plane)][['centrality', 'resolution']]
reso_sel = reso_sel.drop_duplicates('centrality')

if len(reso_sel) == 0:
    print(f"警告：在分辨率文件中找不到匹配的数据 (period={period}, plane={args.plane})")

# 根据粒子类型设置pair_type映射
if args.particle == "Proton":
    SS_types = ["LambdaProton", "LambdaBarProtonBar"]
    OS_types = ["LambdaProtonBar", "LambdaBarProton"]
elif args.particle == "Pion":
    SS_types = ["LambdaPionPos", "LambdaBarPionNeg"]
    OS_types = ["LambdaPionNeg", "LambdaBarPionPos"]
elif args.particle == "Hadron":
    SS_types = ["LambdaHadronPos", "LambdaBarHadronNeg"]
    OS_types = ["LambdaHadronNeg", "LambdaBarHadronPos"]
else:
    print(f"警告：未知的粒子类型 {args.particle}，使用默认的SS/OS分组")
    SS_types = ["Lambda_0", "LambdaBar_1"]
    OS_types = ["Lambda_1", "LambdaBar_0"]

print(f"SS类型: {SS_types}")
print(f"OS类型: {OS_types}")

# 定义加权平均函数
def weighted_mean_and_err(vals, errs, ws):
    if len(vals)==0:
        return float('nan'), float('nan')
    mean = np.sum(ws*vals)/np.sum(ws) if np.sum(ws)>0 else float('nan')
    err = np.sqrt(np.sum(ws**2 * errs**2)/np.sum(ws)**2) if np.sum(ws)>0 else float('nan')
    return mean, err

# 汇总新表
records = []
group_cols = ['centrality', 'diff_type', 'diff_bin']
success_count = 0
failed_reso = 0

for keys, group in main.groupby(group_cols):
    centrality, diff_type, diff_bin = keys

    # 选SS、OS
    SS = group[group['pair_type'].isin(SS_types)]
    OS = group[group['pair_type'].isin(OS_types)]

    # 检查是否有匹配的数据
    if len(SS) == 0:
        print(f"警告：在 centrality={centrality}, diff_type={diff_type}, diff_bin={diff_bin} 中找不到SS类型数据")
    if len(OS) == 0:
        print(f"警告：在 centrality={centrality}, diff_type={diff_type}, diff_bin={diff_bin} 中找不到OS类型数据")

    # 权重
    SS_w = SS['counts_sig'].values
    OS_w = OS['counts_sig'].values

    # SS
    delta_ss, delta_err_ss = weighted_mean_and_err(SS['delta_sig'].values, SS['delta_err_sig'].values, SS_w)
    gamma_ss, gamma_err_ss = weighted_mean_and_err(SS['rawgamma_sig'].values, SS['rawgamma_err_sig'].values, SS_w)
    # OS
    delta_os, delta_err_os = weighted_mean_and_err(OS['delta_sig'].values, OS['delta_err_sig'].values, OS_w)
    gamma_os, gamma_err_os = weighted_mean_and_err(OS['rawgamma_sig'].values, OS['rawgamma_err_sig'].values, OS_w)

    # Del = OS - SS
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

    # 归一化
    gamma_ss_n = gamma_ss / resolution if resolution and not pd.isna(resolution) else float('nan')
    gamma_err_ss_n = gamma_err_ss / resolution if resolution and not pd.isna(resolution) else float('nan')
    gamma_os_n = gamma_os / resolution if resolution and not pd.isna(resolution) else float('nan')
    gamma_err_os_n = gamma_err_os / resolution if resolution and not pd.isna(resolution) else float('nan')
    gamma_del_n = gamma_del / resolution if resolution and not pd.isna(resolution) else float('nan')
    gamma_err_del_n = gamma_err_del / resolution if resolution and not pd.isna(resolution) else float('nan')

    # 输出三类
    for ptype, d, de, g, ge in [
        ("SS", delta_ss, delta_err_ss, gamma_ss_n, gamma_err_ss_n),
        ("OS", delta_os, delta_err_os, gamma_os_n, gamma_err_os_n),
        ("Del", delta_del, delta_err_del, gamma_del_n, gamma_err_del_n)
    ]:
        records.append(dict(
            centrality=centrality,
            diff_type=diff_type,
            diff_bin=diff_bin,
            pair_type=ptype,
            delta=d,
            delta_err=de,
            gamma=g,
            gamma_err=ge
        ))

# 输出
result_df = pd.DataFrame(records)
result_df.to_csv(output_csv, index=False)
print(f"分辨率匹配成功: {success_count}, 失败: {failed_reso}")
print(f"已输出: {output_csv}")
print(f"处理完成！数据集: {args.dataset}, 粒子类型: {args.particle}, 任务: {args.task}")
