import ROOT
import pandas as pd
import re
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='从ROOT文件中提取数据并扁平化为CSV')
parser.add_argument('-i', '--input-root', type=str,
                    default="./AnalysisResults_CVE2025_18q_TPC_Proton.root",
                    help='输入ROOT文件路径')
parser.add_argument('-t', '--task', type=str, default="default",
                    help='分析任务名称，用于定位TDirectoryFile和TList')
parser.add_argument('-o', '--output-dir', type=str, default="./",
                    help='输出目录，将根据数据集和粒子类型自动创建子目录')
parser.add_argument('-m', '--mode', type=str, default="eff_cali",
                    choices=["raw_mass", "eff_cali"],
                    help='处理模式：raw_mass-使用默认处理方式，eff_cali-使用delta3的bin entries作为计数')
args = parser.parse_args()

# 从文件名提取数据集和粒子类型
def extract_info_from_filename(filename):
    # 提取数据集(18q或18r)
    dataset_match = re.search(r'18[qr]', filename)
    dataset = dataset_match.group(0) if dataset_match else "unknown"

    # 提取粒子类型
    particle_types = ["Proton", "Pion", "Hadron"]
    particle_type = "unknown"
    for pt in particle_types:
        if pt in filename:
            particle_type = pt
            break

    return f"LHC{dataset}", particle_type

# 根据粒子类型设置PAIR_TYPE_MAP
def get_pair_type_map(particle_type):
    if particle_type == "Proton":
        return {
            "0": "LambdaProton",
            "1": "LambdaProtonBar",
            "2": "LambdaBarProton",
            "3": "LambdaBarProtonBar"
        }
    elif particle_type == "Pion":
        return {
            "0": "LambdaPionPos",
            "1": "LambdaPionNeg",
            "2": "LambdaBarPionPos",
            "3": "LambdaBarPionNeg"
        }
    elif particle_type == "Hadron":
        return {
            "0": "LambdaHadronPos",
            "1": "LambdaHadronNeg",
            "2": "LambdaBarHadronPos",
            "3": "LambdaBarHadronNeg"
        }
    else:
        print(f"警告：未知的粒子类型 '{particle_type}'，使用通用映射")
        return {
            "0": "Lambda_0",
            "1": "Lambda_1",
            "2": "LambdaBar_0",
            "3": "LambdaBar_1"
        }

# 从文件名提取信息
input_filename = os.path.basename(args.input_root)
dataset, particle_type = extract_info_from_filename(input_filename)
print(f"检测到数据集: {dataset}, 粒子类型: {particle_type}")

# 设置PAIR_TYPE_MAP
PAIR_TYPE_MAP = get_pair_type_map(particle_type)

# 创建输出目录
output_dir = os.path.join(args.output_dir, dataset, particle_type)
os.makedirs(output_dir, exist_ok=True)
output_csv = os.path.join(output_dir, f"flatten_data_{args.task}.csv")

def get_diff_type(name):
    m = re.search(r'Mass(SPt|DEta|Intg)', name)
    return m.group(1) if m else None

def get_pair_type(name):
    m = re.search(r'_(\d)$', name)
    return PAIR_TYPE_MAP.get(m.group(1)) if m else None

# 打开文件
print(f"打开ROOT文件: {args.input_root}")
f = ROOT.TFile.Open(args.input_root)
if not f or f.IsZombie():
    raise RuntimeError(f"无法打开ROOT文件：{args.input_root}")

# 定位到task/ListResults_task
task_dir = f.Get(args.task)
if not task_dir:
    raise RuntimeError(f"未找到 TDirectoryFile '{args.task}'")

list_name = f"ListResults_{args.task}"
lst = task_dir.Get(list_name)
if not lst:
    raise RuntimeError(f"未找到 TList '{list_name}'")

print(f"成功找到 {args.task}/{list_name}")

# 枚举TList中的对象，按照名字分组
obj_map = {}
for i in range(lst.GetSize()):
    obj = lst.At(i)
    name = obj.GetName()
    diff_type = get_diff_type(name)
    pair_type = get_pair_type(name)
    if diff_type and pair_type:
        key = (diff_type, pair_type)
        if name.startswith("fHist3") or name.startswith("fHist2"):
            obj_map.setdefault(key, {})["hist"] = obj
        elif name.startswith("fProfile3DDelta"):
            obj_map.setdefault(key, {})["delta"] = obj
        elif name.startswith("fProfile3DGamma"):
            obj_map.setdefault(key, {})["gamma"] = obj
        elif name.startswith("fProfile3DDiffDelta"):
            # only set if 'delta' not already set
            obj_map.setdefault(key, {}).setdefault("delta", obj)
        elif name.startswith("fProfile3DDiffGamma"):
            obj_map.setdefault(key, {}).setdefault("gamma", obj)

# 过滤掉不完整的分组
valid_combos = {k: v for k, v in obj_map.items() if set(v) == {"hist", "delta", "gamma"}}
print(f"找到 {len(valid_combos)} 个有效的分组")

records = []

for (diff_type, pair_type), group in valid_combos.items():
    h3 = group["hist"]
    delta3 = group["delta"]
    gamma3 = group["gamma"]
    nx = h3.GetNbinsX()
    ny = h3.GetNbinsY()
    nz = h3.GetNbinsZ()

    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            for iz in range(1, nz+1):
                centrality = h3.GetXaxis().GetBinCenter(ix)
                diff_bin = h3.GetYaxis().GetBinCenter(iy)
                inv_mass_center = h3.GetZaxis().GetBinCenter(iz)

                # 根据不同模式设置inv_mass_counts
                if args.mode == "eff_cali":
                    inv_mass_counts = delta3.GetBinEntries(delta3.GetBin(ix, iy, iz))
                elif args.mode == "raw_mass":
                    inv_mass_counts = h3.GetBinContent(ix, iy, iz)
                else:
                    raise ValueError(f"Invalid mode: {args.mode}")

                delta = delta3.GetBinContent(ix, iy, iz)
                delta_err = delta3.GetBinError(ix, iy, iz)
                rawgamma = gamma3.GetBinContent(ix, iy, iz)
                rawgamma_err = gamma3.GetBinError(ix, iy, iz)

                record = dict(
                    centrality=centrality,
                    diff_type=diff_type,
                    diff_bin=diff_bin,
                    pair_type=pair_type,
                    inv_mass_center=inv_mass_center,
                    inv_mass_counts=inv_mass_counts,
                    delta=delta,
                    delta_err=delta_err,
                    rawgamma=rawgamma,
                    rawgamma_err=rawgamma_err
                )
                records.append(record)

df = pd.DataFrame(records)
df.to_csv(output_csv, index=False)
print(f"数据已写入 {output_csv}")
print(f"数据集: {dataset}, 粒子类型: {particle_type}, 任务: {args.task}, 模式: {args.mode}")
