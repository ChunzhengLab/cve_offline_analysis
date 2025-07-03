import ROOT
import pandas as pd
import re
import argparse
import os

# 解析命令行参数
parser = argparse.ArgumentParser(description='从ROOT文件提取Lambda粒子对数据')
parser.add_argument('-i', '--input-root', type=str, required=True,
                   help='输入ROOT文件路径')
parser.add_argument('-t', '--task', type=str, default='default',
                   help='任务名称')
parser.add_argument('-o', '--output-dir', type=str, default='./',
                   help='输出目录')
parser.add_argument('-d', '--dataset', type=str, default=None,
                   help='数据集名称（可选，自动推断）')
parser.add_argument('-m', '--mode', type=str, default="eff_cali",
                   choices=["raw_mass", "eff_cali"],
                   help='处理模式：raw_mass-使用原始质量，eff_cali-使用delta3的bin entries作为计数')
parser.add_argument('--merge4060', action='store_true',
                   help='将中心度45和55的结果合并为50')

args = parser.parse_args()

# 从文件名自动推断数据集
if args.dataset is None:
    filename = os.path.basename(args.input_root)
    if '18q' in filename:
        args.dataset = 'LHC18q'
    elif '18r' in filename:
        args.dataset = 'LHC18r'
    else:
        args.dataset = 'unknown'

# 创建输出目录结构
particle_type = 'Lambda'  # 固定为Lambda
output_path = os.path.join(args.output_dir, args.dataset, particle_type)
os.makedirs(output_path, exist_ok=True)

# 设置输出文件路径
output_csv = os.path.join(output_path, f"flatten_data_{args.task}.csv")
output_root = os.path.join(output_path, f"flatten_data_{args.task}.root")

print(f"输入文件: {args.input_root}")
print(f"任务: {args.task}")
print(f"数据集: {args.dataset}")
print(f"输出目录: {output_path}")

PAIR_TYPE_MAP = {
    "0": "LambdaLambda",
    "1": "LambdaLambdaBar",
    "2": "LambdaBarLambda",
    "3": "LambdaBarLambdaBar"
}

def get_diff_type(name):
    if "MassMass" in name:
        return "MassMass"
    return None

def get_pair_type(name):
    m = re.search(r'_(\d)$', name)
    return PAIR_TYPE_MAP.get(m.group(1)) if m else None

def process_3d_histogram(h3d, p3d_delta, p3d_gamma, output_file, pair_type):
    records = []
    nx = h3d.GetNbinsX()
    ny = h3d.GetNbinsY()
    nz = h3d.GetNbinsZ()
    print(f"处理直方图: {h3d.GetName()}")
    print(f"X轴范围: {h3d.GetXaxis().GetXmin()} - {h3d.GetXaxis().GetXmax()}, bins: {nx}")
    print(f"Y轴范围: {h3d.GetYaxis().GetXmin()} - {h3d.GetYaxis().GetXmax()}, bins: {ny}")
    print(f"Z轴范围: {h3d.GetZaxis().GetXmin()} - {h3d.GetZaxis().GetXmax()}, bins: {nz}")

    # 创建输出目录
    output_file.cd()

    for ix in range(1, nx+1):
        centrality = h3d.GetXaxis().GetBinCenter(ix)
        print(f"\n处理centrality bin {ix}: {centrality}")

        # 为当前centrality创建2D直方图
        h2d_name = f"h2d_{pair_type}_cent_{centrality}"
        h2d = ROOT.TH2D(h2d_name, h2d_name, ny, h3d.GetYaxis().GetXmin(), h3d.GetYaxis().GetXmax(),
                        nz, h3d.GetZaxis().GetXmin(), h3d.GetZaxis().GetXmax())

        p2d_delta_name = f"p2d_delta_{pair_type}_cent_{centrality}"
        p2d_delta = ROOT.TH2D(p2d_delta_name, p2d_delta_name, ny, p3d_delta.GetYaxis().GetXmin(), p3d_delta.GetYaxis().GetXmax(),
                             nz, p3d_delta.GetZaxis().GetXmin(), p3d_delta.GetZaxis().GetXmax())

        p2d_gamma_name = f"p2d_gamma_{pair_type}_cent_{centrality}"
        p2d_gamma = ROOT.TH2D(p2d_gamma_name, p2d_gamma_name, ny, p3d_gamma.GetYaxis().GetXmin(), p3d_gamma.GetYaxis().GetXmax(),
                             nz, p3d_gamma.GetZaxis().GetXmin(), p3d_gamma.GetZaxis().GetXmax())

        # 直接从3D直方图提取数据填充到2D直方图
        for iy in range(1, ny+1):
            for iz in range(1, nz+1):
                # 获取中心值
                inv_mass_1_center = h3d.GetYaxis().GetBinCenter(iy)
                inv_mass_2_center = h3d.GetZaxis().GetBinCenter(iz)

                # 根据不同模式设置inv_mass_counts
                if args.mode == "eff_cali":
                    inv_mass_counts = p3d_delta.GetBinEntries(p3d_delta.GetBin(ix, iy, iz))
                elif args.mode == "raw_mass":
                    inv_mass_counts = h3d.GetBinContent(ix, iy, iz)
                else:
                    raise ValueError(f"Invalid mode: {args.mode}")

                h2d.SetBinContent(iy, iz, inv_mass_counts)

                # 从Profile3D获取delta和gamma值
                delta = p3d_delta.GetBinContent(ix, iy, iz)
                delta_err = p3d_delta.GetBinError(ix, iy, iz)
                p2d_delta.SetBinContent(iy, iz, delta)
                p2d_delta.SetBinError(iy, iz, delta_err)

                rawgamma = p3d_gamma.GetBinContent(ix, iy, iz)
                rawgamma_err = p3d_gamma.GetBinError(ix, iy, iz)
                p2d_gamma.SetBinContent(iy, iz, rawgamma)
                p2d_gamma.SetBinError(iy, iz, rawgamma_err)

                # 只存储非零计数的记录
                if inv_mass_counts > 0:
                    records.append(dict(
                        centrality=centrality,
                        inv_mass_1_center=inv_mass_1_center,
                        inv_mass_2_center=inv_mass_2_center,
                        inv_mass_counts=inv_mass_counts,
                        delta=delta,
                        delta_err=delta_err,
                        rawgamma=rawgamma,
                        rawgamma_err=rawgamma_err,
                        pair_type=pair_type
                    ))

        # 保存直方图
        h2d.Write()
        p2d_delta.Write()
        p2d_gamma.Write()

        # 清理临时对象
        h2d.Delete()
        p2d_delta.Delete()
        p2d_gamma.Delete()

    return records

# 打开文件
f = ROOT.TFile.Open(args.input_root)
if not f or f.IsZombie():
    raise RuntimeError(f"无法打开ROOT文件：{args.input_root}")

print(f"\n成功打开ROOT文件: {args.input_root}")

# 创建输出ROOT文件
out_file = ROOT.TFile(output_root, "RECREATE")
print(f"创建输出ROOT文件: {output_root}")

# 定位到任务目录
task_dir = f.Get(args.task)
if not task_dir:
    raise RuntimeError(f"未找到 TDirectoryFile '{args.task}'")

print(f"找到{args.task}目录")

lst = task_dir.Get(f"ListResults_{args.task}")
if not lst:
    raise RuntimeError(f"未找到 TList 'ListResults_{args.task}'")

print(f"找到ListResults_{args.task}列表")
print(f"列表中的对象数量: {lst.GetSize()}")

# 枚举TList中的对象，按照名字分组
obj_map = {}
for i in range(lst.GetSize()):
    obj = lst.At(i)
    name = obj.GetName()
    print(f"\n检查对象: {name}")
    diff_type = get_diff_type(name)
    pair_type = get_pair_type(name)
    print(f"diff_type: {diff_type}, pair_type: {pair_type}")

    if diff_type and pair_type:
        print(f"找到对象: {name}, diff_type={diff_type}, pair_type={pair_type}")
        if name.startswith("fHist3"):
            obj_map.setdefault((diff_type, pair_type), {})["hist"] = obj
            print("添加为hist对象")
        elif name.startswith("fProfile3DDelta"):
            obj_map.setdefault((diff_type, pair_type), {})["delta"] = obj
            print("添加为delta对象")
        elif name.startswith("fProfile3DGamma"):
            obj_map.setdefault((diff_type, pair_type), {})["gamma"] = obj
            print("添加为gamma对象")

print(f"\n找到的有效组合数量: {len(obj_map)}")

# 过滤掉不完整的分组
valid_combos = {k: v for k, v in obj_map.items() if set(v) == {"hist", "delta", "gamma"}}
print(f"完整的分组数量: {len(valid_combos)}")

all_records = []

for (diff_type, pair_type), group in valid_combos.items():
    print(f"\n处理粒子对: {pair_type}, diff_type: {diff_type}")
    h3 = group["hist"]
    delta3 = group["delta"]
    gamma3 = group["gamma"]

    # 处理每个粒子对的数据
    records = process_3d_histogram(h3, delta3, gamma3, out_file, pair_type)
    print(f"处理得到 {len(records)} 条记录")

    # pair_type信息已经在record创建时添加
    all_records.extend(records)

print(f"\n总共收集到 {len(all_records)} 条记录")

def merge_4060_records(records, dataset):
    # 只对18q/18r数据集生效
    if not (('18q' in dataset) or ('18r' in dataset)):
        return records
    # 权重
    if '18r' in dataset:
        weight = 1883. / 445
    elif '18q' in dataset:
        weight = 1613. / 638
    else:
        return records
    # 按(inv_mass_1_center, inv_mass_2_center, pair_type)分组
    from collections import defaultdict
    group45 = defaultdict(list)
    group55 = defaultdict(list)
    others = []
    for rec in records:
        cent = rec['centrality']
        key = (rec['inv_mass_1_center'], rec['inv_mass_2_center'], rec['pair_type'])
        if abs(cent - 45) < 1e-3:
            group45[key].append(rec)
        elif abs(cent - 55) < 1e-3:
            group55[key].append(rec)
        else:
            others.append(rec)
    # 合并
    merged = []
    for key in set(group45.keys()).union(group55.keys()):
        rec45s = group45.get(key, [])
        rec55s = group55.get(key, [])
        # 只合并完全匹配的bin
        if rec45s and rec55s:
            rec45 = rec45s[0]
            rec55 = rec55s[0]
            c45 = rec45['inv_mass_counts']
            c55 = rec55['inv_mass_counts'] * weight
            total = c45 + c55
            # delta
            delta = (rec45['delta'] * c45 + rec55['delta'] * c55) / total if total > 0 else 0
            delta_err = (rec45['delta_err'] * c45 + rec55['delta_err'] * c55) / total if total > 0 else 0
            # rawgamma
            rawgamma = (rec45['rawgamma'] * c45 + rec55['rawgamma'] * c55) / total if total > 0 else 0
            rawgamma_err = (rec45['rawgamma_err'] * c45 + rec55['rawgamma_err'] * c55) / total if total > 0 else 0
            merged.append({
                'centrality': 50,
                'pair_type': rec45['pair_type'],
                'inv_mass_1_center': rec45['inv_mass_1_center'],
                'inv_mass_2_center': rec45['inv_mass_2_center'],
                'inv_mass_counts': total,
                'delta': delta,
                'delta_err': delta_err,
                'rawgamma': rawgamma,
                'rawgamma_err': rawgamma_err
            })
        elif rec45s:
            # 只有45
            rec45 = rec45s[0]
            merged.append({**rec45, 'centrality': 50})
        elif rec55s:
            # 只有55
            rec55 = rec55s[0]
            c55 = rec55['inv_mass_counts'] * weight
            merged.append({
                'centrality': 50,
                'pair_type': rec55['pair_type'],
                'inv_mass_1_center': rec55['inv_mass_1_center'],
                'inv_mass_2_center': rec55['inv_mass_2_center'],
                'inv_mass_counts': c55,
                'delta': rec55['delta'],
                'delta_err': rec55['delta_err'],
                'rawgamma': rec55['rawgamma'],
                'rawgamma_err': rec55['rawgamma_err']
            })
    return records + merged

if args.merge4060:
    print("\n启用--merge4060功能，正在合并中心度45和55的数据...")
    all_records = merge_4060_records(all_records, args.dataset)
    print(f"合并后记录总数: {len(all_records)}")
    if all_records:
        print("第一条记录的键:", list(all_records[0].keys()))

df = pd.DataFrame(all_records)
print("\nDataFrame的列:", df.columns.tolist())

# 确保列的顺序符合要求
df = df[['centrality', 'pair_type', 'inv_mass_1_center', 'inv_mass_2_center',
         'inv_mass_counts', 'delta', 'delta_err', 'rawgamma', 'rawgamma_err']]
df.to_csv(output_csv, index=False)
print(f"数据已写入 {output_csv}")
print(f"处理模式: {args.mode}")

# 关闭输出ROOT文件
out_file.Close()
print(f"已保存叠加后的直方图到 {output_root}")
