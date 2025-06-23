#!/usr/bin/env python3
import ROOT
import os

def compare_centrality():
    # 打开文件
    f18q = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18q_TPC_Proton.root", "READ")
    f18r = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18r_TPC_Proton.root", "READ")
    
    # 获取TList
    qa_list_18q = f18q.Get("default/ListQA_default")
    qa_list_18r = f18r.Get("default/ListQA_default")
    
    # 获取直方图
    h_cent_18q = qa_list_18q.FindObject("fHistCentAfCut")
    h_cent_18r = qa_list_18r.FindObject("fHistCentAfCut")
    
    # 创建画布
    c = ROOT.TCanvas("c", "Centrality Comparison", 800, 600)
    
    # 使用DrawFrame创建框架
    c.DrawFrame(0, 0, 70, 2e6, "Distribution of centrality;Centrality(%);Event counts")
    
    # 设置直方图属性
    h_cent_18q.SetLineColor(ROOT.kBlue)
    h_cent_18q.SetLineWidth(2)
    h_cent_18r.SetLineColor(ROOT.kRed)
    h_cent_18r.SetLineWidth(2)
    
    # 绘制直方图
    h_cent_18q.Draw("same")
    h_cent_18r.Draw("same")
    
    # 添加图例
    legend = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
    legend.SetBorderSize(0)  # 移除边框
    legend.AddEntry(h_cent_18q, "18q", "l")
    legend.AddEntry(h_cent_18r, "18r", "l")
    legend.Draw()
    
    # 保存图片
    c.SaveAs("centrality_comparison.pdf")
    print("Created: centrality_comparison.pdf")
    
    # 关闭文件
    f18q.Close()
    f18r.Close()

def compare_vz():
    # 打开文件
    f18q = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18q_TPC_Proton.root", "READ")
    f18r = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18r_TPC_Proton.root", "READ")
    
    # 获取TList
    qa_list_18q = f18q.Get("default/ListQA_default")
    qa_list_18r = f18r.Get("default/ListQA_default")
    
    # 获取直方图
    h_vz_18q = qa_list_18q.FindObject("fHistVzAfCut")
    h_vz_18r = qa_list_18r.FindObject("fHistVzAfCut")
    
    # 创建画布
    c = ROOT.TCanvas("c_vz", "Vz Comparison", 800, 600)
    
    # 使用DrawFrame创建框架
    c.DrawFrame(-20, 0, 20, 2e6, "Distribution of Vz;Vz (cm);Event counts")
    
    # 设置直方图属性
    h_vz_18q.SetLineColor(ROOT.kBlue)
    h_vz_18q.SetLineWidth(2)
    h_vz_18r.SetLineColor(ROOT.kRed)
    h_vz_18r.SetLineWidth(2)
    
    # 绘制直方图
    h_vz_18q.Draw("same")
    h_vz_18r.Draw("same")
    
    # 添加图例
    legend = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
    legend.SetBorderSize(0)  # 移除边框
    legend.AddEntry(h_vz_18q, "18q", "l")
    legend.AddEntry(h_vz_18r, "18r", "l")
    legend.Draw()
    
    # 保存图片
    c.SaveAs("vz_comparison.pdf")
    print("Created: vz_comparison.pdf")
    
    # 关闭文件
    f18q.Close()
    f18r.Close()

def apply_efficiency_correction():
    # 打开文件
    f18q = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18q_TPC_Proton.root", "READ")
    f18r = ROOT.TFile("/Users/wangchunzheng/works/Experiments/cve/train_output/NUE/AnalysisResults_CVE2025_18r_TPC_Proton.root", "READ")
    f_nua = ROOT.TFile("../nua_pt/eff_pt_calib_cent.root", "READ")
    
    # 获取TList
    qa_list_18q = f18q.Get("default/ListQA_default")
    qa_list_18r = f18r.Get("default/ListQA_default")
    nua_list = f_nua.Get("fListNUENUA")
    
    # 获取效率修正的TGraph
    eff_graph = nua_list.FindObject("nue_pt_proton_18q_cent0~6")
    
    # 获取pt直方图
    h_pt_18q = qa_list_18q.FindObject("fHistPtAfCut")
    h_pt_18r = qa_list_18r.FindObject("fHistPtAfCut")
    
    # 创建新的直方图用于存储修正后的结果
    h_pt_18q_eff = h_pt_18q.Clone("h_pt_18q_eff_col")
    h_pt_18r_eff = h_pt_18r.Clone("h_pt_18r_eff_col")
    
    # 对每个bin进行效率修正
    for i in range(1, h_pt_18q.GetNbinsX() + 1):
        bin_center = h_pt_18q.GetBinCenter(i)
        efficiency = eff_graph.Eval(bin_center)
        
        # 应用效率修正
        h_pt_18q_eff.SetBinContent(i, h_pt_18q.GetBinContent(i) * efficiency)
        h_pt_18r_eff.SetBinContent(i, h_pt_18r.GetBinContent(i) * efficiency)
    
    # 将修正后的直方图写入文件
    f18q.cd()
    h_pt_18q_eff.Write()
    
    f18r.cd()
    h_pt_18r_eff.Write()
    
    print("Efficiency correction applied and saved to files")
    
    # 关闭文件
    f18q.Close()
    f18r.Close()
    f_nua.Close()

if __name__ == "__main__":
    compare_centrality()
    compare_vz()
    apply_efficiency_correction() 