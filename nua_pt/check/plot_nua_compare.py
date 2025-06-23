#!/usr/bin/env python3
import ROOT
import os
import argparse
import re
from collections import defaultdict
import numpy as np

# 设置ROOT画图样式
def set_root_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadGridX(True)
    ROOT.gStyle.SetPadGridY(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "xyz")
    ROOT.gStyle.SetTitleFont(42, "xyz")
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)

# 粒子类型名称映射
particle_names = {
    "poshadron": "Positive Hadron",
    "neghadron": "Negative Hadron",
    "proton": "Proton",
    "antiproton": "Anti-Proton",
    "lambda": "Lambda",
    "antilambda": "Anti-Lambda"
}

def extract_info_from_name(hist_name):
    """从直方图名称中提取信息"""
    # 检查是否包含pt范围
    pt_match = re.search(r"_pt(\d+\.\d+)_(\d+\.\d+)", hist_name)
    if pt_match:
        pt_low = float(pt_match.group(1))
        pt_high = float(pt_match.group(2))
        pt_info = f"{pt_low:.1f} < p_{{T}} < {pt_high:.1f} GeV/c"
    else:
        pt_info = "All p_{T}"

    # 提取粒子类型和状态
    if "Proton" in hist_name and "Anti" not in hist_name:
        particle = "proton"
    elif "AntiProton" in hist_name:
        particle = "antiproton"
    elif "Lambda" in hist_name and "Anti" not in hist_name:
        particle = "lambda"
    elif "AntiLambda" in hist_name:
        particle = "antilambda"
    else:
        particle = "poshadron"

    state = "before" if "bfNUA" in hist_name else "after"

    return {"particle": particle, "state": state, "pt_info": pt_info}

def plot_comparison(input_file, output_dir="plots"):
    """读取ROOT文件并创建比较图"""
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打开ROOT文件
    root_file = ROOT.TFile.Open(input_file, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Cannot open {input_file}")
        return

    # 获取所有直方图的键
    all_keys = [key.GetName() for key in root_file.GetListOfKeys()]

    # 按粒子类型和pt分组
    histogram_groups = defaultdict(dict)

    for key in all_keys:
        # 获取直方图
        hist = root_file.Get(key)
        if not hist:
            continue

        # 从名称中提取信息
        info = extract_info_from_name(key)
        particle = info["particle"]
        state = info["state"]
        pt_info = info["pt_info"]

        # 组织直方图
        group_key = f"{particle}_{pt_info}"
        histogram_groups[group_key][state] = hist

    # 设置ROOT样式
    set_root_style()

    # 为每组创建一个画布并绘制比较图
    for group_key, hists in histogram_groups.items():
        # 检查是否同时有before和after状态的直方图
        if "before" not in hists or "after" not in hists:
            continue

        # 获取直方图
        hist_before = hists["before"]
        hist_after = hists["after"]

        # 提取粒子名称和pt信息
        particle, pt_info = group_key.split("_", 1)
        particle_name = particle_names.get(particle, particle)

        # 创建画布
        canvas_name = f"canvas_{group_key}"
        canvas = ROOT.TCanvas(canvas_name, f"{particle_name} {pt_info}", 800, 600)

        # 添加比率面板
        canvas.cd()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
        pad1.SetBottomMargin(0.02)
        pad1.Draw()
        pad1.cd()

        # 获取x轴和y轴范围
        xmin = hist_before.GetXaxis().GetXmin()
        xmax = hist_before.GetXaxis().GetXmax()
        max_val = max(hist_before.GetMaximum(), hist_after.GetMaximum()) * 1.1

        # 使用DrawFrame创建框架
        frame1 = pad1.DrawFrame(xmin, 0, xmax, max_val,
                              f"{particle_name} #phi Distribution - {pt_info}")
        frame1.GetXaxis().SetLabelSize(0)  # 隐藏上面板x轴标签
        frame1.GetYaxis().SetTitle("Normalized Counts")
        frame1.GetYaxis().SetTitleSize(0.05)
        frame1.GetYaxis().SetTitleOffset(1.2)

        # 设置直方图样式
        hist_before.SetLineColor(ROOT.kBlue)
        hist_before.SetMarkerColor(ROOT.kBlue)
        hist_before.SetMarkerStyle(20)
        hist_before.SetMarkerSize(1.0)
        hist_before.SetLineWidth(2)
        hist_before.SetTitle("")  # 清除标题以避免重复

        hist_after.SetLineColor(ROOT.kRed)
        hist_after.SetMarkerColor(ROOT.kRed)
        hist_after.SetMarkerStyle(21)
        hist_after.SetMarkerSize(0.8)
        hist_after.SetLineWidth(2)
        hist_after.SetTitle("")  # 清除标题以避免重复

        # 绘制直方图
        hist_before.Draw("HIST SAME")
        hist_after.Draw("HIST SAME")

        # 添加图例
        legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(hist_before, "Before NUA", "l")
        legend.AddEntry(hist_after, "After NUA", "l")
        legend.Draw()

        # 创建比率直方图
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        pad2.cd()

        ratio_hist = hist_after.Clone(f"ratio_{group_key}")
        ratio_hist.Divide(hist_before)

        # 设置比率范围
        min_ratio = max(0.7, ratio_hist.GetMinimum() * 0.9)
        max_ratio = min(1.3, ratio_hist.GetMaximum() * 1.1)

        # 使用DrawFrame创建框架
        frame2 = pad2.DrawFrame(xmin, min_ratio, xmax, max_ratio, "")
        frame2.GetYaxis().SetTitle("After/Before")
        frame2.GetYaxis().SetTitleSize(0.12)
        frame2.GetYaxis().SetTitleOffset(0.4)
        frame2.GetYaxis().SetLabelSize(0.10)
        frame2.GetXaxis().SetLabelSize(0.10)
        frame2.GetXaxis().SetTitleSize(0.12)
        frame2.GetXaxis().SetTitle("#phi (radians)")
        frame2.GetXaxis().SetTitleOffset(1.0)

        ratio_hist.SetLineColor(ROOT.kBlack)
        ratio_hist.SetMarkerColor(ROOT.kBlack)
        ratio_hist.SetMarkerStyle(20)
        ratio_hist.SetMarkerSize(0.6)
        ratio_hist.SetTitle("")

        # 绘制参考线
        line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
        line.SetLineColor(ROOT.kGray+2)
        line.SetLineStyle(2)
        line.Draw()

        # 绘制比率直方图
        ratio_hist.Draw("HIST SAME")

        # 保存图像
        canvas.SaveAs(f"{output_dir}/{particle}_{pt_info.replace(' < ', '_').replace(' > ', '_').replace(' ', '_').replace('p_{T}', 'pt').replace('/c', '')}.pdf")
        print(f"Created plot: {output_dir}/{particle}_{pt_info.replace(' < ', '_').replace(' > ', '_').replace(' ', '_')}.pdf")

    # 关闭文件
    root_file.Close()

def plot_particle_summary(input_file, output_dir="plots"):
    """为每种粒子创建一个总结图（所有pt区间）"""
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打开ROOT文件
    root_file = ROOT.TFile.Open(input_file, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Cannot open {input_file}")
        return

    # 按粒子类型收集所有全pt投影
    particle_histograms = defaultdict(dict)

    # 获取所有直方图的键
    all_keys = [key.GetName() for key in root_file.GetListOfKeys()]

    # 筛选全pt投影直方图（带有_projY但不带_pt的）
    for key in all_keys:
        if "_projY" in key and not "_pt" in key:
            # 获取直方图
            hist = root_file.Get(key)
            if not hist:
                continue

            # 从名称中提取信息
            info = extract_info_from_name(key)
            particle = info["particle"]
            state = info["state"]

            # 保存直方图
            particle_histograms[particle][state] = hist

    # 设置ROOT样式
    set_root_style()

    # 为每种粒子创建一个画布
    for particle, hists in particle_histograms.items():
        # 检查是否同时有before和after状态的直方图
        if "before" not in hists or "after" not in hists:
            continue

        # 获取直方图
        hist_before = hists["before"]
        hist_after = hists["after"]

        # 获取粒子名称
        particle_name = particle_names.get(particle, particle)

        # 创建画布
        canvas_name = f"canvas_{particle}_summary"
        canvas = ROOT.TCanvas(canvas_name, f"{particle_name} Summary", 800, 600)

        # 添加比率面板
        canvas.cd()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
        pad1.SetBottomMargin(0.02)
        pad1.Draw()
        pad1.cd()

        # 获取x轴和y轴范围
        xmin = hist_before.GetXaxis().GetXmin()
        xmax = hist_before.GetXaxis().GetXmax()
        max_val = max(hist_before.GetMaximum(), hist_after.GetMaximum()) * 1.1

        # 使用DrawFrame创建框架
        frame1 = pad1.DrawFrame(xmin, 0, xmax, max_val,
                              f"{particle_name} #phi Distribution - All p_{{T}}")
        frame1.GetXaxis().SetLabelSize(0)  # 隐藏上面板x轴标签
        frame1.GetYaxis().SetTitle("Normalized Counts")
        frame1.GetYaxis().SetTitleSize(0.05)
        frame1.GetYaxis().SetTitleOffset(1.2)

        # 设置直方图样式
        hist_before.SetLineColor(ROOT.kBlue)
        hist_before.SetMarkerColor(ROOT.kBlue)
        hist_before.SetMarkerStyle(20)
        hist_before.SetMarkerSize(1.0)
        hist_before.SetLineWidth(2)
        hist_before.SetTitle("")  # 清除标题以避免重复

        hist_after.SetLineColor(ROOT.kRed)
        hist_after.SetMarkerColor(ROOT.kRed)
        hist_after.SetMarkerStyle(21)
        hist_after.SetMarkerSize(0.8)
        hist_after.SetLineWidth(2)
        hist_after.SetTitle("")  # 清除标题以避免重复

        # 绘制直方图
        hist_before.Draw("HIST SAME")
        hist_after.Draw("HIST SAME")

        # 添加图例
        legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.85)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.AddEntry(hist_before, "Before NUA", "l")
        legend.AddEntry(hist_after, "After NUA", "l")
        legend.Draw()

        # 创建比率直方图
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        pad2.cd()

        ratio_hist = hist_after.Clone(f"ratio_{particle}_summary")
        ratio_hist.Divide(hist_before)

        # 设置比率范围
        min_ratio = max(0.3, ratio_hist.GetMinimum() * 0.9)
        max_ratio = min(3, ratio_hist.GetMaximum() * 1.1)

        # 使用DrawFrame创建框架
        frame2 = pad2.DrawFrame(xmin, min_ratio, xmax, max_ratio, "")
        frame2.GetYaxis().SetTitle("After/Before")
        frame2.GetYaxis().SetTitleSize(0.12)
        frame2.GetYaxis().SetTitleOffset(0.4)
        frame2.GetYaxis().SetLabelSize(0.10)
        frame2.GetXaxis().SetLabelSize(0.10)
        frame2.GetXaxis().SetTitleSize(0.12)
        frame2.GetXaxis().SetTitle("#phi (radians)")
        frame2.GetXaxis().SetTitleOffset(1.0)

        ratio_hist.SetLineColor(ROOT.kBlack)
        ratio_hist.SetMarkerColor(ROOT.kBlack)
        ratio_hist.SetMarkerStyle(20)
        ratio_hist.SetMarkerSize(0.6)
        ratio_hist.SetTitle("")

        # 绘制参考线
        line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
        line.SetLineColor(ROOT.kGray+2)
        line.SetLineStyle(2)
        line.Draw()

        # 绘制比率直方图
        ratio_hist.Draw("HIST SAME")

        # 保存图像
        canvas.SaveAs(f"{output_dir}/{particle}_summary.pdf")
        print(f"Created summary plot: {output_dir}/{particle}_summary.pdf")

    # 关闭文件
    root_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='绘制NUA校正前后对比图')
    parser.add_argument('input_file', help='输入ROOT文件路径 (nua_compare.root)')
    parser.add_argument('--output-dir', default="plots", help='输出图像目录 (默认: plots)')

    args = parser.parse_args()

    # 绘制所有比较图
    plot_comparison(args.input_file, args.output_dir)

    # 绘制粒子总结图
    plot_particle_summary(args.input_file, args.output_dir)

    print(f"所有图形已保存到 {args.output_dir} 目录")
