#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成效率图和权重图，并保存为ROOT文件。
"""

import ROOT
import math
from array import array

ROOT.gROOT.SetBatch(True)  # 不弹出Canvas
ROOT.TH1.AddDirectory(False)

# ------------------------------------------------------------
# 1. 自定义 pT rebin 边界
# ------------------------------------------------------------
pt_bins_dict = {
    "proton"     : [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "antiproton" : [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "poshadron"  : [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "neghadron"  : [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "pospion_tpc_region" : [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
    "negpion_tpc_region" : [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
    "pospion_tof_region" : [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "negpion_tof_region" : [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.6, 5.0],
    "lambda"     : [0.5, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.45, 1.6, 1.75, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 4.0, 4.4, 4.9, 5.5, 6.3, 7.1, 8.0, 9.0, 10.0],
    "antilambda" : [0.5, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.45, 1.6, 1.75, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6, 4.0, 4.4, 4.9, 5.5, 6.3, 7.1, 8.0, 9.0, 10.0]
}

# 中心度区间定义
cent_ranges = ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%"]

def rebin(h, edges):
    """根据给定的边界值对直方图进行重新分bin"""
    new_edges = array('d', edges)
    nb = len(edges) - 1
    return h.Rebin(nb, h.GetName() + "_rebin", new_edges)

# ------------------------------------------------------------
# 主程序开始
# ------------------------------------------------------------

# 定义要处理的输入文件
input_files = {
    "18q": "../data/mc_18q.root",
    "18r": "../data/mc_18r.root"
}

# 粒子列表
particles = [
    "poshadron", "neghadron",
    "pospion_tpc_region", "negpion_tpc_region",
    "pospion_tof_region", "negpion_tof_region",
    "proton", "antiproton", "lambda", "antilambda"
]

# 创建输出文件
out_eff = ROOT.TFile("eff_pt_cent.root", "RECREATE")
out_nue = ROOT.TFile("eff_pt_calib_cent.root", "RECREATE")

# 创建TList容器
eff_list = ROOT.TList()
eff_list.SetName("fListEFF")
eff_list.SetOwner(False)

nue_list = ROOT.TList()
nue_list.SetName("fListNUE")
nue_list.SetOwner(False)

# 处理两个输入文件
for period, input_file_path in input_files.items():
    # 打开输入文件
    input_file = ROOT.TFile.Open(input_file_path)
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open input file {input_file_path}")
        continue

    # 打印ROOT文件内容以便调试
    print(f"Input file contents for {period}:")
    input_file.ls()

    # 获取ResultsList_default
    reslist = input_file.Get("default/ResultsList_default")
    if not reslist:
        print(f"Error: Could not find default/ResultsList_default in the input file {input_file_path}")
        input_file.Close()
        continue

    # 主循环：处理每种粒子
    for particle in particles:
        edges = pt_bins_dict[particle]

        # 从ResultsList_default获取直方图
        h_mc_name = f"h2_pt_mc_{particle}"
        h_rc_name = f"h2_pt_rc_real_{particle}"

        # Fix syntax error - Python uses 'or' instead of '||'
        if particle == "pospion_tpc_region" or particle == "pospion_tof_region":
            h_mc_name = "h2_pt_mc_pospion"
            h_rc_name = "h2_pt_rc_real_pospion"
        if particle == "negpion_tpc_region" or particle == "negpion_tof_region":
            h_mc_name = "h2_pt_mc_negpion"
            h_rc_name = "h2_pt_rc_real_negpion"

        h_mc = reslist.FindObject(h_mc_name)
        h_rc = reslist.FindObject(h_rc_name)

        print(f"Looking for histograms in {period}: {h_mc_name} and {h_rc_name}")

        # 检查直方图是否存在
        if not h_mc or not h_rc:
            print(f"[Error] skip {particle} in {period}: missing histograms")
            continue

        # 检查是否为TH2类型
        if not isinstance(h_mc, ROOT.TH2) or not isinstance(h_rc, ROOT.TH2):
            print(f"[Error] skip {particle} in {period}: histograms are not TH2D")
            continue

        # 处理每个中心度bin
        for cent_bin in range(1, 8):  # 7个中心度bin，从1开始计数
            # 从TH2D提取特定中心度bin的投影
            h_mc_proj = h_mc.ProjectionY(f"{particle}_{period}_mc_cent{cent_bin-1}_proj", cent_bin, cent_bin)
            h_rc_proj = h_rc.ProjectionY(f"{particle}_{period}_rc_cent{cent_bin-1}_proj", cent_bin, cent_bin)

            print(f"Projection for {particle}, {period}, cent{cent_bin-1}: MC entries={h_mc_proj.GetEntries()}, RC entries={h_rc_proj.GetEntries()}")

            # 如果投影为空，跳过
            if h_mc_proj.GetEntries() == 0 or h_rc_proj.GetEntries() == 0:
                print(f"[Error] skip {particle} {period} cent{cent_bin-1}: empty projection")
                continue

            # 重新分bin
            h_mc_reb = rebin(h_mc_proj, edges)
            h_rc_reb = rebin(h_rc_proj, edges)

            # 准备TGraphErrors数据
            nb = h_mc_reb.GetNbinsX()
            x  = array('d')
            y  = array('d')
            ex = array('d')
            ey = array('d')
            y_w = array('d')
            ey_w = array('d')

            bad = False
            for ib in range(1, nb+1):
                passed = h_rc_reb.GetBinContent(ib)
                total  = h_mc_reb.GetBinContent(ib)

                # 检查数据有效性
                if passed > total:
                    print(f"[WARN] bin {ib} in {particle} {period} cent{cent_bin-1}: passed ({passed}) > total ({total})")
                    bad = True
                    break

                cen  = h_mc_reb.GetBinCenter(ib)
                wid  = h_mc_reb.GetBinWidth(ib)/2.0

                # 计算效率和误差
                eff  = passed/total if total>0 else 0.
                err  = math.sqrt(eff*(1-eff)/total) if 0<eff<1 and total>0 else 0.

                # 添加数据点
                x.append(cen)
                ex.append(wid)
                y.append(eff)
                ey.append(err)

                # 计算权重和误差
                w     = 1./eff if eff>0 else 0.
                err_w = err/(eff**2) if eff>0 else 0.
                y_w.append(w)
                ey_w.append(err_w)

            if bad:
                print(f"[WARN] skip {particle} {period} cent{cent_bin-1}: passed>total")
                continue

            # 创建效率图
            g_eff = ROOT.TGraphErrors(nb, x, y, ex, ey)
            g_eff_name = f"eff_pt_{particle}_{period}_cent{cent_bin-1}"
            g_eff.SetName(g_eff_name)
            g_eff.SetTitle(f"Efficiency {particle} {period} {cent_ranges[cent_bin-1]};p_T;#varepsilon")

            # 创建权重图
            g_w = ROOT.TGraphErrors(nb, x, y_w, ex, ey_w)
            g_w_name = f"nue_pt_{particle}_{period}_cent{cent_bin-1}"
            g_w.SetName(g_w_name)
            g_w.SetTitle(f"1/eff {particle} {period} {cent_ranges[cent_bin-1]};p_T;weight")

            # 添加到输出列表
            eff_list.Add(g_eff)
            nue_list.Add(g_w)

            print(f" --> added {g_eff.GetName()} to eff_pt_cent.root and {g_w.GetName()} to eff_pt_calib_cent.root")

    # 关闭当前输入文件
    input_file.Close()

# 将图形写入ROOT文件
out_eff.cd()
eff_list.Write("fListEFF", ROOT.TObject.kSingleKey)
out_eff.Close()

out_nue.cd()
nue_list.Write("fListNUE", ROOT.TObject.kSingleKey)
out_nue.Close()

print("\nEfficiency graphs saved to eff_pt_cent.root (key: fListEFF)")
print("Weight graphs saved to eff_pt_calib_cent.root (key: fListNUE)")
