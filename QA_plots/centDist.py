#!/usr/bin/env python3
import ROOT
import os
import sys

def create_pileup_plots(input_dir):
    # 获取目录下所有的root文件
    root_files = [f for f in os.listdir(input_dir) if f.endswith('.root')]
    
    for root_file in root_files:
        file_path = os.path.join(input_dir, root_file)
        print(f"Processing file: {file_path}")
        
        # 打开ROOT文件
        f = ROOT.TFile(file_path, "READ")
        
        # 获取TList
        qa_list = f.Get("default/ListQA_default")
        if not qa_list:
            print(f"Warning: Could not find ListQA_default in {file_path}")
            continue
            
        # 获取所需的直方图
        h_v0m_spd1_bf = qa_list.FindObject("fHist2CentQA_V0M_SPD1_BfCut")
        h_v0m_spd1_af = qa_list.FindObject("fHist2CentQA_V0M_SPD1_AfCut")
        h_mult_cent_bf = qa_list.FindObject("fHist2MultCentQA_BfCut")
        h_mult_cent_af = qa_list.FindObject("fHist2MultCentQA_AfCut")
        
        # 检查所有直方图是否存在
        histograms = {
            "V0M vs SPD1 Before": h_v0m_spd1_bf,
            "V0M vs SPD1 After": h_v0m_spd1_af,
            "Multiplicity vs Centrality Before": h_mult_cent_bf,
            "Multiplicity vs Centrality After": h_mult_cent_af
        }
        
        missing_hists = [name for name, hist in histograms.items() if not hist]
        if missing_hists:
            print(f"Warning: The following histograms are missing in {file_path}:")
            for name in missing_hists:
                print(f"  - {name}")
            continue
        
        # 设置直方图属性
        for hist in [h_v0m_spd1_bf, h_v0m_spd1_af]:
            hist.SetTitle("")
            hist.GetXaxis().SetTitle("Centrality V0M (%)")
            hist.GetYaxis().SetTitle("Centrality SPD1 (%)")
            hist.SetStats(0)  # 关闭统计框
            # 设置坐标轴标题偏移
            hist.GetXaxis().SetTitleOffset(1.2)
            hist.GetYaxis().SetTitleOffset(1.4)
            # 设置边距
            hist.GetXaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetXaxis().SetTitleSize(0.04)
            hist.GetYaxis().SetTitleSize(0.04)
            
        for hist in [h_mult_cent_bf, h_mult_cent_af]:
            hist.SetTitle("")
            hist.GetXaxis().SetTitle("Centrality V0M (%)")
            hist.GetYaxis().SetTitle("Multiplicity (FB32)")
            hist.SetStats(0)  # 关闭统计框
            # 设置坐标轴标题偏移
            hist.GetXaxis().SetTitleOffset(1.2)
            hist.GetYaxis().SetTitleOffset(1.4)
            # 设置边距
            hist.GetXaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetXaxis().SetTitleSize(0.04)
            hist.GetYaxis().SetTitleSize(0.04)
            # 设置Y轴为科学计数法
            hist.GetYaxis().SetNdivisions(505)
            # 设置科学计数法显示
            ROOT.TGaxis.SetMaxDigits(3)
        
        # 创建第一个画布 - V0M vs SPD1
        c1 = ROOT.TCanvas("c1", "V0M vs SPD1 Comparison", 1200, 600)
        c1.Divide(2, 1)
        
        # 设置画布边距
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.15)
        c1.SetBottomMargin(0.15)
        c1.SetTopMargin(0.15)
        
        # 绘制V0M vs SPD1图
        c1.cd(1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetTopMargin(0.15)
        h_v0m_spd1_bf.SetTitle("Before pile up rejection")
        h_v0m_spd1_bf.Draw("colz")
        ROOT.gPad.SetLogz()
        
        c1.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetTopMargin(0.15)
        h_v0m_spd1_af.SetTitle("After pile up rejection")
        h_v0m_spd1_af.Draw("colz")
        ROOT.gPad.SetLogz()
        
        # 保存第一个画布
        output_name1 = f"{os.path.splitext(root_file)[0]}_v0m_spd1.pdf"
        c1.SaveAs(output_name1)
        print(f"Created: {output_name1}")
        
        # 创建第二个画布 - Multiplicity vs Centrality
        c2 = ROOT.TCanvas("c2", "Multiplicity vs Centrality Comparison", 1200, 600)
        c2.Divide(2, 1)
        
        # 设置画布边距
        c2.SetLeftMargin(0.15)
        c2.SetRightMargin(0.15)
        c2.SetBottomMargin(0.15)
        c2.SetTopMargin(0.15)
        
        # 绘制Multiplicity vs Centrality图
        c2.cd(1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetTopMargin(0.15)
        h_mult_cent_bf.SetTitle("Before pile up rejection")
        h_mult_cent_bf.Draw("colz")
        ROOT.gPad.SetLogz()
        
        c2.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetTopMargin(0.15)
        h_mult_cent_af.SetTitle("After pile up rejection")
        h_mult_cent_af.Draw("colz")
        ROOT.gPad.SetLogz()
        
        # 保存第二个画布
        output_name2 = f"{os.path.splitext(root_file)[0]}_mult_cent.pdf"
        c2.SaveAs(output_name2)
        print(f"Created: {output_name2}")
        
        # 关闭文件
        f.Close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python pileUp.py <input_directory>")
        sys.exit(1)
        
    input_dir = sys.argv[1]
    create_pileup_plots(input_dir)
