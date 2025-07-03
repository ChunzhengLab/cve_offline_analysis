import ROOT
import pandas as pd
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='从ROOT文件中提取Lambda和LambdaBar的三维直方图数据并扁平化为CSV')
parser.add_argument('-i', '--input-root', type=str, required=True, help='输入ROOT文件路径')
parser.add_argument('-t', '--task', type=str, required=True, help='分析任务名称，用于定位TDirectoryFile和TList')
parser.add_argument('-o', '--output-dir', type=str, default="./", help='输出目录')
parser.add_argument('--merge4060', action='store_true', help='将中心度45和55的结果合并为50')
args = parser.parse_args()

# 判断数据集后缀
input_filename = os.path.basename(args.input_root)
if '18q' in input_filename:
    dataset_suffix = '_18q'
elif '18r' in input_filename:
    dataset_suffix = '_18r'
else:
    dataset_suffix = '_unknowndata'

# 创建输出目录
os.makedirs(args.output_dir, exist_ok=True)
output_csv_base = f"flatten_data_{args.task}{dataset_suffix}"
if args.merge4060:
    output_csv_base += "_merge4060"
output_csv = os.path.join(args.output_dir, output_csv_base + ".csv")

# 打开ROOT文件
print(f"打开ROOT文件: {args.input_root}")
f = ROOT.TFile.Open(args.input_root)
if not f or f.IsZombie():
    raise RuntimeError(f"无法打开ROOT文件：{args.input_root}")

task_dir = f.Get(args.task)
if not task_dir:
    raise RuntimeError(f"未找到 TDirectoryFile '{args.task}'")

list_name = f"ListQA_{args.task}"
lst = task_dir.Get(list_name)
if not lst:
    raise RuntimeError(f"未找到 TList '{list_name}'")

print(f"成功找到 {args.task}/{list_name}")

# 获取两个TH3D对象
hist_names = {
    "fHist3LambdaCentPtMassWeighted": "Lambda",
    "fHist3AntiLambdaCentPtMassWeighted": "LambdaBar"
}

records = []

for hist_name, particle in hist_names.items():
    h3 = lst.FindObject(hist_name)
    if not h3:
        print(f"警告：未找到 {hist_name}")
        continue
    h3 = h3.RebinZ(2) # 60 --> 30
    nx = h3.GetNbinsX()
    ny = h3.GetNbinsY()
    nz = h3.GetNbinsZ()
    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            for iz in range(1, nz+1):
                centrality = h3.GetXaxis().GetBinCenter(ix)
                pT_mean = h3.GetYaxis().GetBinCenter(iy)
                inv_mass_mean = h3.GetZaxis().GetBinCenter(iz)
                inv_mass_count = h3.GetBinContent(ix, iy, iz)
                record = dict(
                    centrality=centrality,
                    particle=particle,
                    pT_mean=pT_mean,
                    inv_mass_mean=inv_mass_mean,
                    inv_mass_count=inv_mass_count
                )
                records.append(record)

def merge_4060_records(records, dataset):
    if not (('18q' in dataset) or ('18r' in dataset)):
        return records
    if '18r' in dataset:
        weight = 1883. / 445
    elif '18q' in dataset:
        weight = 1613. / 638
    else:
        return records
    from collections import defaultdict
    group45 = defaultdict(list)
    group55 = defaultdict(list)
    key_set_45 = set()
    key_set_55 = set()
    for rec in records:
        cent = rec['centrality']
        key = (rec['particle'], rec['pT_mean'], rec['inv_mass_mean'])
        if abs(cent - 45) < 1e-3:
            group45[key].append(rec)
            key_set_45.add(key)
        elif abs(cent - 55) < 1e-3:
            group55[key].append(rec)
            key_set_55.add(key)
    print(f"\n[DEBUG] 45 keys (原始): {len(key_set_45)}")
    for k in list(key_set_45)[:5]:
        print(f"  45 key: {k}")
    print(f"[DEBUG] 55 keys (原始): {len(key_set_55)}")
    for k in list(key_set_55)[:5]:
        print(f"  55 key: {k}")
    key_set_45_round = set((k[0], round(k[1],2), round(k[2],2)) for k in key_set_45)
    key_set_55_round = set((k[0], round(k[1],2), round(k[2],2)) for k in key_set_55)
    print(f"[DEBUG] 45 keys (round2): {len(key_set_45_round)}")
    for k in list(key_set_45_round)[:5]:
        print(f"  45 key (r2): {k}")
    print(f"[DEBUG] 55 keys (round2): {len(key_set_55_round)}")
    for k in list(key_set_55_round)[:5]:
        print(f"  55 key (r2): {k}")
    print(f"[DEBUG] 45&55交集 (原始): {len(key_set_45 & key_set_55)}")
    print(f"[DEBUG] 45&55交集 (round2): {len(key_set_45_round & key_set_55_round)}")
    merged = []
    for key in set(group45.keys()).union(group55.keys()):
        rec45s = group45.get(key, [])
        rec55s = group55.get(key, [])
        if rec45s and rec55s:
            rec45 = rec45s[0]
            rec55 = rec55s[0]
            c45 = rec45['inv_mass_count']
            c55 = rec55['inv_mass_count'] * weight
            total = c45 + c55
            print(f"[DEBUG] 合并key: {key}, c45={c45}, c55={c55}, total={total}")
            merged.append({
                'centrality': 50,
                'particle': rec45['particle'],
                'pT_mean': rec45['pT_mean'],
                'inv_mass_mean': rec45['inv_mass_mean'],
                'inv_mass_count': total
            })
        elif rec45s:
            rec45 = rec45s[0]
            print(f"[DEBUG] 只有45: {key}, c45={rec45['inv_mass_count']}")
            merged.append({**rec45, 'centrality': 50})
        elif rec55s:
            rec55 = rec55s[0]
            c55 = rec55['inv_mass_count'] * weight
            print(f"[DEBUG] 只有55: {key}, c55={c55}")
            merged.append({
                'centrality': 50,
                'particle': rec55['particle'],
                'pT_mean': rec55['pT_mean'],
                'inv_mass_mean': rec55['inv_mass_mean'],
                'inv_mass_count': c55
            })
    return records + merged

# 输出为CSV
if records:
    df = pd.DataFrame(records)
    df = df[df['pT_mean'] >= 1]  # 只保留pT_mean >= 1的数据
    if args.merge4060:
        print("\n启用--merge4060功能，正在合并中心度45和55的数据...")
        df = pd.DataFrame(merge_4060_records(df.to_dict('records'), dataset_suffix))
        print(f"合并后记录总数: {len(df)}")
    df.to_csv(output_csv, index=False)
    print(f"数据已写入 {output_csv}")
else:
    print("未找到任何有效数据，未生成CSV文件。")
