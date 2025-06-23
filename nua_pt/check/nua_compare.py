import ROOT
import os
import array


pt_ranges = {
    "poshadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "neghadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "proton":     (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "antiproton": (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "lambda":     (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0),
    "antilambda": (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0)
}

# Mapping of histogram names and particle types
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


def process_histograms(input_file, rebin=2, output_file_name="nua_compare.root"):
    """从ROOT文件中提取TH2D直方图，投影处理并保存到新文件

    Args:
        input_file (str): 输入ROOT文件路径
        rebin (int): 对投影后的TH1D进行rebin的参数
        output_file_name (str): 输出ROOT文件名
    """

    # 打开ROOT文件
    root_file = ROOT.TFile.Open(input_file)
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

    # 创建输出文件
    output_file = ROOT.TFile(output_file_name, "RECREATE")

    # 处理每个直方图
    for hist_name, info in histogram_info.items():
        # 从列表中查找直方图
        hist2d = None

        # 在列表中查找直方图
        for i in range(lst.GetEntries()):
            obj = lst.At(i)
            if not obj:
                continue

            obj_name = obj.GetName()
            if obj_name == hist_name:
                hist2d = obj
                break

        # 如果找不到直方图，跳过
        if not hist2d:
            print(f"Could not find histogram {hist_name}")
            continue

        # 确认是TH2D类型
        if not hist2d.InheritsFrom("TH2"):
            print(f"{hist_name} is not a TH2D histogram, skipping")
            continue

        particle_type = info["particle"]
        state = info["state"]

        # 创建整体的投影
        proj_y = hist2d.ProjectionY(f"{hist_name}_projY", 1, hist2d.GetNbinsX())

        # 在缩放前先进行rebin
        proj_y.Rebin(rebin)

        # 缩放整体投影
        if proj_y.Integral() > 0:
            scale_factor = proj_y.GetNbinsX() / proj_y.Integral()
            proj_y.Scale(scale_factor)

        # 禁用统计误差
        for bin in range(1, proj_y.GetNbinsX() + 1):
            proj_y.SetBinError(bin, 0)

        # 保存整体投影
        output_file.cd()
        proj_y.Write(f"{hist_name}_projY")

        # 对每个pt区间进行处理
        pt_bins = pt_ranges[particle_type]
        for i in range(len(pt_bins) - 1):
            pt_low = pt_bins[i]
            pt_high = pt_bins[i+1]

            # 找到对应的bin范围
            bin_low = hist2d.GetXaxis().FindBin(pt_low)
            bin_high = hist2d.GetXaxis().FindBin(pt_high) - 1

            # 这个pt区间的名称
            pt_slice_name = f"{hist_name}_pt{pt_low:.1f}_{pt_high:.1f}"

            # 创建这个pt区间的phi分布
            pt_slice = hist2d.ProjectionY(pt_slice_name, bin_low, bin_high)

            # 在缩放前先进行rebin
            pt_slice.Rebin(rebin)

            # 缩放这个pt区间的投影
            if pt_slice.Integral() > 0:
                pt_scale_factor = pt_slice.GetNbinsX() / pt_slice.Integral()
                pt_slice.Scale(pt_scale_factor)

            # 禁用统计误差
            for bin in range(1, pt_slice.GetNbinsX() + 1):
                pt_slice.SetBinError(bin, 0)

            # 保存这个pt区间的投影
            output_file.cd()
            pt_slice.Write(pt_slice_name)

            print(f"Processed {pt_slice_name}")

    # 关闭输出文件
    output_file.Close()

    print(f"All histograms processed and saved to nua_compare.root")

if __name__ == "__main__":
    input_file = "AnalysisResults_CVE2025_18q_TPC_Proton.root"
    process_histograms(input_file)
