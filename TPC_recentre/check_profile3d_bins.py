#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_profile3d_bins.py

读取 ROOT 文件里 fListTPC (TList) 中的所有 TProfile3D，
统计每个 3D bin 的 Entries 与 Content，
并把它们各自做成分布直方图写出。
"""

import sys
import ROOT

ROOT.gROOT.SetBatch(False)  # 改成 True 可完全后台运行

# ---------------------- 解析命令行 ----------------------
if len(sys.argv) < 2:
    print("Usage: python check_profile3d_bins.py  input.root  [output.root]")
    sys.exit(1)

in_path  = sys.argv[1]
out_path = sys.argv[2] if len(sys.argv) > 2 else "bin_distributions.root"

# ---------------------- 打开输入文件 ---------------------
fin = ROOT.TFile.Open(in_path, "READ")
if not fin or fin.IsZombie():
    raise RuntimeError(f"❌ Cannot open {in_path}")

fListTPC = fin.Get("fListTPC")
if not fListTPC or not fListTPC.InheritsFrom("TList"):
    raise RuntimeError("❌ Cannot find TList 'fListTPC' in the file")

print(f"✅ Loaded fListTPC with {fListTPC.GetSize()} objects")

# ---------------------- 先收集所有数值 -------------------
entries_values  = []   # 所有 bin 的 Entries
content_values  = []   # 所有 bin 的 Content

for obj in fListTPC:
    if not obj.InheritsFrom("TProfile3D"):
        print(f"⚠️  Skip non-TProfile3D object: {obj.GetName()}")
        continue

    nbx, nby, nbz = obj.GetNbinsX(), obj.GetNbinsY(), obj.GetNbinsZ()

    for ix in range(1, nbx + 1):
        for iy in range(1, nby + 1):
            for iz in range(1, nbz + 1):
                bin_idx  = obj.GetBin(ix, iy, iz)
                entries  = obj.GetBinEntries(bin_idx)
                content  = obj.GetBinContent(bin_idx)

                entries_values.append(entries)
                content_values.append(content)

print(f"🎯 Total bins scanned: {len(entries_values)}")

# ---------------------- 创建并填充直方图 -----------------
# 使用可自动扩展的 TH1F
hEntries = ROOT.TH1F("hEntries",
                     "Distribution of TProfile3D bin entries;Entries per bin;Number of bins",
                     1000, 0, 1e5)
# hEntries.SetDirectory(0)

hContent = ROOT.TH1F("hContent",
                     "Distribution of TProfile3D bin contents;Content;Number of bins",
                     1000, 0, 0.05)

h2EntriesContent = ROOT.TH2F("hEntriesContent",
                     "Distribution of TProfile3D bin entries and contents;Entries per bin;Content",
                     1000, 0, 1e5, 1000, 0, 0.05)

# print(entries_values)
for v in entries_values:
    print(v)
    hEntries.Fill(v)
for v in content_values:
    hContent.Fill(v)
for v1, v2 in zip(entries_values, content_values):
    h2EntriesContent.Fill(v1, v2)

# ---------------------- 写入输出文件 ---------------------
fout = ROOT.TFile.Open(out_path, "RECREATE")
hEntries.Write()
hContent.Write()
h2EntriesContent.Write()
fout.Close()
fin.Close()

print(f"✅ Histograms written to {out_path}")

# # ---------------------- 可视化（非 batch 模式下） --------
# if not ROOT.gROOT.IsBatch():
#     c1 = ROOT.TCanvas("cEntries", "Entries distribution", 800, 600)
#     hEntries.Draw("hist")
#     c2 = ROOT.TCanvas("cContent", "Content distribution", 800, 600)
#     hContent.Draw("hist")
#     input("🖱️  Press <Enter> to quit...")   # 保持画布
