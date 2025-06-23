#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
from array import array
import os

# 设置ROOT样式
ROOT.gROOT.SetBatch(True)  # 不显示画布
ROOT.gStyle.SetOptStat(0)  # 不显示统计框
ROOT.gStyle.SetPadTickX(1)  # 显示X轴刻度
ROOT.gStyle.SetPadTickY(1)  # 显示Y轴刻度
ROOT.gStyle.SetLegendBorderSize(0)  # 图例无边框

# 创建输出目录
output_dir = "plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 读取输入文件
eff_file = ROOT.TFile.Open("eff_pt_cent.root")
nue_file = ROOT.TFile.Open("eff_pt_calib_cent.root")
pt_file = ROOT.TFile.Open("./pt_dist.root")  # 添加pT分布文件

cent_ranges = ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%"]
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+1, ROOT.kViolet]

# 获取列表
eff_list = eff_file.Get("fListEFF")
nue_list = nue_file.Get("fListNUE")

datasets = ["18q", "18r"]
particles = [
    "poshadron", "neghadron",
    # "pospion_tpc_region", "negpion_tpc_region",
    # "pospion_tof_region", "negpion_tof_region",
    "proton", "antiproton", "lambda", "antilambda"
]

# 粒子显示名称映射
particle_display_names = {
    "poshadron": "h^{+}",
    "neghadron": "h^{-}",
    # "pospion_tpc_region": "#pi^{+} (TPC)",
    # "negpion_tpc_region": "#pi^{-} (TPC)",
    # "pospion_tof_region": "#pi^{+} (TOF)",
    # "negpion_tof_region": "#pi^{-} (TOF)",
    "proton": "p",
    "antiproton": "#bar{p}",
    "lambda": "#Lambda",
    "antilambda": "#bar{#Lambda}"
}

# 粒子标记样式
particle_markers = {
    "poshadron": 20,
    "neghadron": 20,
    # "pospion_tpc_region": 22,
    # "negpion_tpc_region": 23,
    # "pospion_tof_region": 24,
    # "negpion_tof_region": 25,
    "proton": 20,
    "antiproton": 20,
    "lambda": 20,
    "antilambda": 20
}

# 粒子颜色
particle_colors = {
    "poshadron": ROOT.kRed,
    "neghadron": ROOT.kRed+2,
    # "pospion_tpc_region": ROOT.kBlue,
    # "negpion_tpc_region": ROOT.kBlue+2,
    # "pospion_tof_region": ROOT.kGreen+2,
    # "negpion_tof_region": ROOT.kGreen+4,
    "proton": ROOT.kMagenta,
    "antiproton": ROOT.kMagenta+2,
    "lambda": ROOT.kOrange,
    "antilambda": ROOT.kOrange+2
}

# 定义获取图形的函数
def get_graph(list_obj, particle, dataset, cent):
    """从TList中获取指定的TGraphErrors"""
    graph_name = f"eff_pt_{particle}_{dataset}_cent{cent}"
    graph = list_obj.FindObject(graph_name)
    if not graph:
        print(f"Warning: Graph {graph_name} not found")
        return None
    return graph

def get_graph_nue(list_obj, particle, dataset, cent):
    """从TList中获取指定的 NUE TGraphErrors"""
    graph_name = f"nue_pt_{particle}_{dataset}_cent{cent}"
    graph = list_obj.FindObject(graph_name)
    if not graph:
        print(f"Warning: NUE Graph {graph_name} not found")
        return None
    return graph

# 定义样式设置函数
def set_graph_style(graph, color, marker_style=20, marker_size=1.2, line_width=2):
    """设置图形样式"""
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(marker_size)
    graph.SetLineWidth(line_width)

# 创建图例的函数
def create_legend(x1, y1, x2, y2, header=None):
    """创建图例"""
    legend = ROOT.TLegend(x1, y1, x2, y2)
    if header:
        legend.SetHeader(header)
    legend.SetTextSize(0.035)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    return legend

# 创建比值图的函数
def create_ratio_graph(g_num, g_den):
    """创建比值图 g_num/g_den"""
    if not g_num or not g_den:
        return None

    # 获取点数
    n_num = g_num.GetN()
    n_den = g_den.GetN()

    # 创建新的图形
    g_ratio = ROOT.TGraphErrors()
    g_ratio.SetName(f"{g_num.GetName()}_ratio")

    # 计算比值
    point_idx = 0
    for i in range(n_num):
        x_num = array('d', [0.0])
        y_num = array('d', [0.0])
        g_num.GetPoint(i, x_num, y_num)
        e_num = g_num.GetErrorY(i)

        for j in range(n_den):
            x_den = array('d', [0.0])
            y_den = array('d', [0.0])
            g_den.GetPoint(j, x_den, y_den)
            e_den = g_den.GetErrorY(j)

            if abs(x_num[0] - x_den[0]) < 1e-6 and y_den[0] != 0:
                ratio = y_num[0] / y_den[0]
                err_ratio = ratio * math.sqrt((e_num / y_num[0]) ** 2 + (e_den / y_den[0]) ** 2) if y_num[0] != 0 else 0
                g_ratio.SetPoint(point_idx, x_num[0], ratio)
                g_ratio.SetPointError(point_idx, 0, err_ratio)
                point_idx += 1
                break

                g_ratio.SetPoint(point_idx, x_num, ratio)
                g_ratio.SetPointError(point_idx, 0, err_ratio)
                point_idx += 1
                break

    return g_ratio

# 1. 画18q和18r在cent3下所有粒子的效率对比图
def plot_all_particles_cent3():
    """为18q和18r在cent3的所有粒子创建效率对比图"""
    cent_idx = 3  # 对应cent3

    canvas = ROOT.TCanvas("c_cent3_comparison", "Cent3 Comparison", 1200, 600)
    canvas.Divide(2, 1)  # 左右两个子图

    legends = []
    frames = []

    for i, dataset in enumerate(datasets):
        canvas.cd(i+1)

        # 创建Frame
        frame = ROOT.gPad.DrawFrame(0, 0, 5, 1.3, f"Efficiency {dataset};p_{{T}} (GeV/c);#varepsilon")
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetTitleSize(0.045)
        frames.append(frame)

        # 创建图例
        legend = create_legend(0.5, 0.65, 0.85, 0.85)
        legends.append(legend)

        # 画所有粒子
        for particle in particles:
            graph = get_graph(eff_list, particle, dataset, cent_idx)
            if graph:
                set_graph_style(graph, particle_colors[particle], particle_markers[particle])
                graph.Draw("P SAME")
                legend.AddEntry(graph, particle_display_names[particle], "P")

        legend.Draw()

    # 保存图像
    canvas.SaveAs(f"{output_dir}/all_particles_cent3_comparison.pdf")

import math

# 2. 为每个dataset和粒子创建所有中心度效率图及其与cent0的比值
def plot_centrality_comparison():
    for particle in particles:
        x_max = 10.0 if particle in ["lambda", "antilambda"] else 5.0

        canvas = ROOT.TCanvas(f"c_{particle}_split", f"{particle} split 18q/18r", 1800, 800)
        canvas.Divide(2, 1)

        all_ratio_graphs = []  # 关键1：在大循环外面存活

        # 为每个pad分别画18q/18r
        for idx, dataset in enumerate(datasets):  # idx=0是18q，idx=1是18r
            canvas.cd(idx + 1)

            pad1 = ROOT.TPad(f"pad1_{dataset}", "pad1", 0.0, 0.3, 1.0, 1.0)
            pad1.SetBottomMargin(0.0)
            pad1.SetLeftMargin(0.13)
            pad1.SetRightMargin(0.06)
            pad1.Draw()
            pad1.cd()

            frame1 = pad1.DrawFrame(0, 0, x_max, 1.05,
                                    f"Efficiency {particle_display_names[particle]} {dataset};p_{{T}} (GeV/c);#varepsilon")
            frame1.GetXaxis().SetTitleSize(0.045)
            frame1.GetYaxis().SetTitleSize(0.045)
            legend = create_legend(0.65, 0.45, 0.95, 0.85)

            for cent_idx in range(len(cent_ranges)):
                graph = get_graph(eff_list, particle, dataset, cent_idx)
                if graph:
                    set_graph_style(graph, colors[cent_idx], 20)
                    graph.Draw("P SAME")
                    legend.AddEntry(graph, cent_ranges[cent_idx], "P")
            legend.Draw()
            pad1.Update()

            # 下pad（比值）
            canvas.cd(idx + 1)
            pad2 = ROOT.TPad(f"pad2_{dataset}", "pad2", 0.0, 0.0, 1.0, 0.3)
            pad2.SetTopMargin(0.0)
            pad2.SetBottomMargin(0.28)
            pad2.SetLeftMargin(0.13)
            pad2.SetRightMargin(0.06)
            pad2.Draw()
            pad2.cd()

            frame2 = pad2.DrawFrame(0, 0.8, x_max, 3,
                                    ";p_{{T}} (GeV/c);Ratio to Cent 0-10%")
            frame2.GetXaxis().SetTitleSize(0.1)
            frame2.GetYaxis().SetTitleSize(0.1)
            frame2.GetXaxis().SetLabelSize(0.09)
            frame2.GetYaxis().SetLabelSize(0.09)
            unity = ROOT.TLine(0, 1, x_max, 1)
            unity.SetLineStyle(2)
            unity.SetLineWidth(2)
            unity.SetLineColor(ROOT.kBlack)
            unity.Draw("SAME")

            g_ref = get_graph(eff_list, particle, dataset, 0)
            ratio_graphs = []
            for cent_idx in range(len(cent_ranges)):
                graph = get_graph(eff_list, particle, dataset, cent_idx)
                if graph:
                    ratio_graph = create_ratio_graph(graph, g_ref)
                    if ratio_graph:
                        set_graph_style(ratio_graph, colors[cent_idx], 20)
                        ratio_graphs.append(ratio_graph)
                        all_ratio_graphs.append(ratio_graph)  # 关键2：引用给大list，绝不会被销毁！

            # 画全部
            for ratio_graph in ratio_graphs:
                ratio_graph.Draw("P SAME")
            pad2.Update()

        canvas.Update()
        canvas.SaveAs(f"{output_dir}/{particle}_split_18q_18r_comparison.pdf")

# 3 画18q和18r在cent3下所有粒子的NUE Weight对比图
def plot_all_particles_cent3_nue():
    """为18q和18r在cent3的所有粒子画 NUEWeight 对比图"""
    cent_idx = 3  # 对应cent3

    canvas = ROOT.TCanvas("c_cent3_nue_comparison", "Cent3 NUE Comparison", 1200, 600)
    canvas.Divide(2, 1)  # 左右两个子图

    legends = []
    frames = []

    for i, dataset in enumerate(datasets):
        canvas.cd(i+1)

        # 创建Frame
        frame = ROOT.gPad.DrawFrame(0, 1, 5, 1000, f"NUE Weight {dataset};p_{{T}} (GeV/c);NUE")
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetTitleSize(0.045)
        frames.append(frame)
        ROOT.gPad.SetLogy()

        # 创建图例
        legend = create_legend(0.5, 0.65, 0.85, 0.85)
        legends.append(legend)

        # 画所有粒子
        for particle in particles:
            graph = get_graph_nue(nue_list, particle, dataset, cent_idx)
            if graph:
                set_graph_style(graph, particle_colors[particle], particle_markers[particle])
                graph.Draw("L SAME")
                legend.AddEntry(graph, particle_display_names[particle], "L")

        legend.Draw()

    # 保存图像
    canvas.SaveAs(f"{output_dir}/all_particles_cent3_nue_comparison.pdf")

# 新增：绘制pT分布的函数
def plot_pt_distributions():
    """绘制所有粒子的pT分布"""
    canvas = ROOT.TCanvas("c_pt_dist", "pT Distributions", 1800, 800)
    canvas.Divide(2, 1)

    for idx, dataset in enumerate(datasets):
        canvas.cd(idx + 1)
        
        # 创建Frame
        frame = ROOT.gPad.DrawFrame(0, 0, 5, 1.1, f"pT Distribution {dataset};p_{{T}} (GeV/c);Counts(a.u.)")
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetTitleSize(0.045)
        
        # 创建图例
        legend = create_legend(0.65, 0.65, 0.95, 0.95)
        
        # 画所有粒子
        for particle in particles:
            hist_name = f"h_pt_{particle}_{dataset}_cent3_eff_col"
            h = pt_file.Get(hist_name)
            if h:
                # 计算积分
                integral = h.Integral()
                if integral > 0:
                    # Scale by Nxbins/Integral
                    h.Scale(0.15 * h.GetNbinsX() / integral)
                    
                    # 设置样式
                    h.SetLineColor(particle_colors[particle])
                    h.SetMarkerColor(particle_colors[particle])
                    h.SetMarkerStyle(particle_markers[particle])
                    h.SetMarkerSize(1.2)
                    h.SetLineWidth(2)
                    
                    # 绘制
                    h.Draw("HIST SAME")
                    legend.AddEntry(h, particle_display_names[particle], "L")
        
        legend.Draw()
        # canvas.cd(idx + 1).SetLogy()
    
    canvas.SaveAs(f"{output_dir}/pt_distributions.pdf")

# 执行绘图
plot_all_particles_cent3()
plot_centrality_comparison()
plot_all_particles_cent3_nue()
plot_pt_distributions()  # 添加新的pT分布绘图

print("Plots have been saved to the 'plots' directory.")
