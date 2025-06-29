import ROOT
import pandas as pd
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='从ROOT文件中提取Lambda和LambdaBar的三维直方图数据并扁平化为CSV')
parser.add_argument('-i', '--input-root', type=str, required=True, help='输入ROOT文件路径')
parser.add_argument('-t', '--task', type=str, required=True, help='分析任务名称，用于定位TDirectoryFile和TList')
parser.add_argument('-o', '--output-dir', type=str, default="./", help='输出目录')
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
output_csv = os.path.join(args.output_dir, f"flatten_data_{args.task}{dataset_suffix}.csv")

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

# 输出为CSV
if records:
    df = pd.DataFrame(records)
    df = df[df['pT_mean'] >= 1]  # 只保留pT_mean >= 1的数据
    df.to_csv(output_csv, index=False)
    print(f"数据已写入 {output_csv}")
else:
    print("未找到任何有效数据，未生成CSV文件。")
