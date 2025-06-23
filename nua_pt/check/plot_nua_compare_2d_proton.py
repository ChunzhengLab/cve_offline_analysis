#!/usr/bin/env python3
import ROOT
import os
import argparse
import numpy as np
from array import array

# 定义pt_ranges，与nua_compare.py保持一致
pt_ranges = {
    "poshadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "neghadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "proton":     (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "antiproton": (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "lambda":     (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0),
    "antilambda": (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0)
}

# 定义直方图映射，与nua_compare.py保持一致
histogram_info = {
    "fHistPhi_bfNUA": {"particle": "poshadron", "state": "before"},
    "fHistPhi_afNUA": {"particle": "poshadron", "state": "after"},
    "fHistProtonPhi_bfNUA": {"particle": "proton", "state": "before"},
    "fHistProtonPhi_afNUA": {"particle": "proton", "state": "after"},
    "fHistAntiProtonPhi_bfNUA": {"particle": "antiproton", "state": "before"},
    "fHistAntiProtonPhi_afNUA": {"particle": "antiproton", "state": "after"},
    "fHistLambdaPhi_bfNUA": {"particle": "lambda", "state": "before"},
    "fHistLambdaPhi_afNUA": {"particle": "lambda", "state": "after"},
    "fHistAntiLambdaPhi_bfNUA": {"particle": "antilambda", "state": "before"},
    "fHistAntiLambdaPhi_afNUA": {"particle": "antilambda", "state": "after"}
}

# 设置ROOT画图样式
def set_root_style():
    # 设置ROOT样式
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "xyz")
    ROOT.gStyle.SetTitleFont(42, "xyz")
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)  # 增大右边距以适应色标
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.08)
    # 设置调色板
    ROOT.gStyle.SetPalette(ROOT.kBird)  # 使用较好的彩虹色调色板
    ROOT.gStyle.SetNumberContours(100)  # 增加颜色层次

def process_and_plot_2d(input_file, output_dir="plots"):
    """处理质子和反质子的TH2D直方图，归一化每个pT slice并绘制比较图

    Args:
        input_file: 输入ROOT文件路径
        output_dir: 输出图像目录
    """
    # 启用ROOT批处理模式以提高性能并避免弹出窗口
    ROOT.gROOT.SetBatch(True)
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打开ROOT文件
    root_file = ROOT.TFile.Open(input_file, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Cannot open {input_file}")
        return

    # 获取结果列表
    lst = root_file.Get("default/ListQA_default")
    if not lst:
        print(f"Cannot find ListQA_default, trying alternatives...")
        lst = root_file.Get("default/ResultsList_default")
    if not lst:
        print(f"Cannot find any list, exiting")
        return

    # 要处理的粒子类型
    particles_to_process = ["proton", "antiproton"]

    # 存储找到的直方图
    histograms = {}

    # 查找直方图
    for hist_name, info in histogram_info.items():
        # 只处理质子和反质子
        if info["particle"] not in particles_to_process:
            continue

        for i in range(lst.GetEntries()):
            obj = lst.At(i)
            if not obj:
                continue

            obj_name = obj.GetName()
            if obj_name == hist_name:
                histograms[hist_name] = obj
                print(f"Found histogram: {hist_name}")
                break

        if hist_name not in histograms:
            print(f"Could not find histogram {hist_name}")

    # 设置ROOT样式
    set_root_style()

    # 处理质子和反质子
    for particle_type in particles_to_process:
        # 获取对应的直方图名称
        hist_before_name = None
        hist_after_name = None

        for hist_name, info in histogram_info.items():
            if info["particle"] == particle_type:
                if info["state"] == "before":
                    hist_before_name = hist_name
                elif info["state"] == "after":
                    hist_after_name = hist_name

        print(f"Processing {particle_type}: before={hist_before_name}, after={hist_after_name}")

        # 检查直方图是否存在
        if hist_before_name not in histograms or hist_after_name not in histograms:
            print(f"Skipping {particle_type} due to missing histograms")
            continue

        hist_before = histograms[hist_before_name]
        hist_after = histograms[hist_after_name]

        # 创建归一化的2D直方图
        norm_before = normalize_2d_hist(hist_before)
        if not norm_before:
            print(f"Error: Failed to normalize histogram for {particle_type} (before)")
            continue

        norm_after = normalize_2d_hist(hist_after)
        if not norm_after:
            print(f"Error: Failed to normalize histogram for {particle_type} (after)")
            continue

        # 创建输出ROOT文件
        output_root_path = f"{output_dir}/{particle_type}_nua_compare.root"
        output_root = ROOT.TFile(output_root_path, "RECREATE")
        print(f"Creating output ROOT file: {output_root_path}")

        # 获取粒子显示名称
        particle_display_name = particle_type.capitalize()
        if particle_type == "antiproton":
            particle_display_name = "Anti-Proton"

        # 将直方图保存到ROOT文件中
        output_root.cd()

        # 保存原始直方图
        hist_before.Write(f"{particle_type}_before")
        hist_after.Write(f"{particle_type}_after")

        # 保存归一化直方图
        norm_before.Write(f"{particle_type}_normalized_before")
        norm_after.Write(f"{particle_type}_normalized_after")

        # 直接访问和保存归一化过程中创建的平均值直方图
        # 不依赖Python的属性管理，而是在ROOT文件中重新创建
        mean_hist_before = ROOT.TH1D(f"{particle_type}_mean_values_before",
                                    f"{particle_display_name} Phi Bin Mean Values - Before",
                                    hist_before.GetNbinsY(),
                                    hist_before.GetYaxis().GetXmin(),
                                    hist_before.GetYaxis().GetXmax())

        mean_hist_after = ROOT.TH1D(f"{particle_type}_mean_values_after",
                                   f"{particle_display_name} Phi Bin Mean Values - After",
                                   hist_after.GetNbinsY(),
                                   hist_after.GetYaxis().GetXmin(),
                                   hist_after.GetYaxis().GetXmax())

        # 直接从归一化过程中重新计算并填充平均值
        # 这样避免依赖Python的属性存储
        calculate_phi_means(hist_before, mean_hist_before)
        calculate_phi_means(hist_after, mean_hist_after)

        # 保存平均值直方图
        mean_hist_before.Write()
        mean_hist_after.Write()

        # 创建对比画布
        canvas = ROOT.TCanvas(f"canvas_{particle_type}", f"{particle_display_name} NUA Comparison", 1200, 500)
        canvas.Divide(2, 1)  # 左右分割

        # 绘制校正前的直方图
        canvas.cd(1)
        norm_before.SetTitle(f"{particle_display_name} #phi vs p_{{T}} - Before NUA")
        norm_before.GetXaxis().SetTitle("p_{T} (GeV/c)")
        norm_before.GetYaxis().SetTitle("#phi (radians)")
        norm_before.GetZaxis().SetTitle("Normalized Counts")
        norm_after.GetZaxis().SetRangeUser(0, 2)
        norm_before.Draw("COLZ")

        # 绘制校正后的直方图
        canvas.cd(2)
        norm_after.SetTitle(f"{particle_display_name} #phi vs p_{{T}} - After NUA")
        norm_after.GetXaxis().SetTitle("p_{T} (GeV/c)")
        norm_after.GetYaxis().SetTitle("#phi (radians)")
        norm_after.GetZaxis().SetTitle("Normalized Counts")
        norm_after.GetZaxis().SetRangeUser(0, 2)
        norm_after.Draw("COLZ")

        # 保存画布
        canvas.SaveAs(f"{output_dir}/{particle_type}_2d_nua_comparison.pdf")
        print(f"Created 2D plot: {output_dir}/{particle_type}_2d_nua_comparison.pdf")

        # 将画布也保存到ROOT文件中
        canvas.Write(f"{particle_type}_canvas")

    # 关闭ROOT文件
    output_root.Close()
    root_file.Close()

    print(f"All plots and ROOT files saved to {output_dir}")

def normalize_2d_hist(hist2d):
    """对2D直方图的每个pt slice进行归一化处理，使每个pt段内的phi分布归一化到1附近

    对每个pt段(x轴)，计算该pt段中所有phi角度的平均值，然后用该pt段内各phi bin的值除以这个平均值。
    这样处理后，每个pt段内的平均值接近1，可以突显phi角度分布在每个pt段上的不均匀性。

    返回一个新的归一化后的2D直方图。
    """
    # 检查输入直方图是否有效
    if not hist2d or not isinstance(hist2d, ROOT.TH2):
        print(f"Error: Invalid histogram provided")
        return None

    print(f"Normalizing histogram: {hist2d.GetName()}")

    # 获取直方图中的粒子类型
    particle_type = None
    for hist_name, info in histogram_info.items():
        if hist2d.GetName() == hist_name:
            particle_type = info["particle"]
            break

    if not particle_type:
        print(f"Warning: Could not determine particle type for {hist2d.GetName()}")
        particle_type = "proton"  # 默认值

    # 创建一个与原始直方图完全相同的新直方图
    hist_name = hist2d.GetName() + "_normalized"
    normalized_hist = hist2d.Clone(hist_name)
    normalized_hist.Reset()  # 清空所有bin的内容
    normalized_hist.SetDirectory(0)  # 防止被ROOT内存管理回收

    # 获取直方图维度norm_before.GetZaxis().SetRangeUser(0, 2)
    nbins_x = hist2d.GetNbinsX()
    nbins_y = hist2d.GetNbinsY()

    print(f"  Dimensions: {nbins_x} x {nbins_y} bins")

    # 对每个pt slice (x轴)单独进行归一化
    for i in range(1, nbins_x + 1):  # 遍历pt bin (x轴)
        # 计算该pt段内所有phi角度的平均值
        total_content = 0
        count = 0

        for j in range(1, nbins_y + 1):  # 遍历该pt段内的所有phi bin
            content = hist2d.GetBinContent(i, j)
            if content > 0:
                total_content += content
                count += 1

        # 计算该pt段的平均值
        mean_val = total_content / count if count > 0 else 1.0

        # 对该pt段内的所有phi bin进行归一化
        for j in range(1, nbins_y + 1):
            original_content = hist2d.GetBinContent(i, j)

            # 只处理非零值
            if original_content > 0 and mean_val > 0:
                normalized_value = original_content / mean_val
            else:
                normalized_value = 1.0  # 设置为中性值

            normalized_hist.SetBinContent(i, j, normalized_value)
            normalized_hist.SetBinError(i, j, 0)  # 设置误差为0

    # 额外检查并处理异常值
    for i in range(1, nbins_x + 1):
        for j in range(1, nbins_y + 1):
            bin_content = normalized_hist.GetBinContent(i, j)
            if np.isnan(bin_content) or np.isinf(bin_content) or bin_content <= 0:
                normalized_hist.SetBinContent(i, j, 1.0)  # 设置为中性值

    print(f"  Normalization complete for {hist2d.GetName()}")
    return normalized_hist

def calculate_phi_means(hist2d, mean_hist):
    """计算2D直方图中每个pt slice内phi bin的平均值并填充到1D直方图中"""
    nbins_x = hist2d.GetNbinsX()
    nbins_y = hist2d.GetNbinsY()

    # 对每个pt slice计算平均值
    for i in range(1, nbins_x + 1):
        total_content = 0
        count = 0

        for j in range(1, nbins_y + 1):
            content = hist2d.GetBinContent(i, j)
            if content > 0:
                total_content += content
                count += 1

        # 计算该pt slice的平均值
        mean_val = total_content / count if count > 0 else 0

        # 填充到平均值直方图中 - 注意这里我们使用pt bin索引
        mean_hist.SetBinContent(i, mean_val)
        mean_hist.SetBinError(i, 0)  # 设置误差为0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理并绘制质子和反质子的2D NUA比较图')
    parser.add_argument('input_file', help='输入ROOT文件路径')
    parser.add_argument('--output-dir', default="plots", help='输出图像目录 (默认: plots)')

    args = parser.parse_args()

    process_and_plot_2d(args.input_file, args.output_dir)
