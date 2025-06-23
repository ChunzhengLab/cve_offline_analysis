#!/usr/bin/env python3
import ROOT
from ROOT import TMath, TGraphErrors
from array import array
import csv
import os

def get_tgraph(file, name):
    obj = file.Get(name)
    if not obj:
        print(f"未找到 {name}")
        print("当前文件对象有：")
        for key in file.GetListOfKeys():
            print(key.GetName(), key.GetClassName())
        exit(1)
    # 确保是 TGraphErrors 类型
    if obj.ClassName() != "TGraphErrors":
        print(f"{name} 类型错误: 实际类型是 {obj.ClassName()}")
        exit(1)
    obj_noerr = ROOT.TGraph(obj)
    return obj_noerr

def calculate_ratio(h_total, h_component, pt_min=0):
    """计算比例和误差"""
    bin_min = h_total.FindBin(pt_min)
    nbins = h_total.GetNbinsX()
    
    total_integral = h_total.Integral(bin_min, nbins)
    component_integral = h_component.Integral(bin_min, nbins)
    
    total_err = array('d', [0])
    component_err = array('d', [0])
    h_total.IntegralAndError(bin_min, nbins, total_err)
    h_component.IntegralAndError(bin_min, nbins, component_err)
    
    ratio = component_integral / total_integral if total_integral > 0 else 0
    ratio_err = ratio * TMath.Sqrt((component_err[0]/component_integral)**2 + (total_err[0]/total_integral)**2) if total_integral > 0 else 0
    
    return ratio, ratio_err

def process_centrality(cent):
    """处理特定中心度的数据"""
    # 关闭全局统计框
    ROOT.gStyle.SetOptStat(0)
    
    # --- 打开输入文件 ---
    f_pt = ROOT.TFile.Open("./pt_dist.root")
    f_frac = ROOT.TFile.Open("./src_fraction.root")

    # --- 读取pT分布 ---
    h_pt_proton_18q = f_pt.Get(f"h_pt_proton_18q_cent{cent}_eff_col")
    h_pt_proton_18r = f_pt.Get(f"h_pt_proton_18r_cent{cent}_eff_col")
    h_pt_antiproton_18q = f_pt.Get(f"h_pt_antiproton_18q_cent{cent}_eff_col")
    h_pt_antiproton_18r = f_pt.Get(f"h_pt_antiproton_18r_cent{cent}_eff_col")

    if not all([h_pt_proton_18q, h_pt_proton_18r, h_pt_antiproton_18q, h_pt_antiproton_18r]):
        print(f"未找到中心度 {cent} 所需的pT分布")
        return

    # --- 获取 TGraphErrors 分数 ---
    # 所有中心度都使用cent1的分数
    g_frac_p_primary_18q = get_tgraph(f_frac, "g_primary_proton_18q_cent1")
    g_frac_p_secondary_18q = get_tgraph(f_frac, "g_secondary_proton_18q_cent1")
    g_frac_p_primary_18r = get_tgraph(f_frac, "g_primary_proton_18r_cent1")
    g_frac_p_secondary_18r = get_tgraph(f_frac, "g_secondary_proton_18r_cent1")
    
    g_frac_ap_primary_18q = get_tgraph(f_frac, "g_primary_antiproton_18q_cent1")
    g_frac_ap_secondary_18q = get_tgraph(f_frac, "g_secondary_antiproton_18q_cent1")
    g_frac_ap_primary_18r = get_tgraph(f_frac, "g_primary_antiproton_18r_cent1")
    g_frac_ap_secondary_18r = get_tgraph(f_frac, "g_secondary_antiproton_18r_cent1")

    # --- 画fraction图 ---
    c3 = ROOT.TCanvas(f"c3_cent{cent}", f"Fraction Distribution (Centrality {cent})", 1200, 600)
    c3.Divide(2, 1)

    # 创建图例
    leg_frac = ROOT.TLegend(0.3, 0.4, 0.8, 0.6)
    leg_frac.SetBorderSize(0)

    # 绘制18q的fraction
    c3.cd(1)
    g_frac_p_primary_18q.SetLineColor(ROOT.kRed)
    g_frac_p_secondary_18q.SetLineColor(ROOT.kRed)
    g_frac_ap_primary_18q.SetLineColor(ROOT.kBlue)
    g_frac_ap_secondary_18q.SetLineColor(ROOT.kBlue)
    
    g_frac_p_primary_18q.SetMarkerColor(ROOT.kRed)
    g_frac_p_secondary_18q.SetMarkerColor(ROOT.kRed)
    g_frac_ap_primary_18q.SetMarkerColor(ROOT.kBlue)
    g_frac_ap_secondary_18q.SetMarkerColor(ROOT.kBlue)

    g_frac_p_primary_18q.SetMarkerStyle(20)
    g_frac_p_secondary_18q.SetMarkerStyle(24)
    g_frac_ap_primary_18q.SetMarkerStyle(21)
    g_frac_ap_secondary_18q.SetMarkerStyle(25)
    
    g_frac_p_primary_18q.Draw("ALP")
    g_frac_p_secondary_18q.Draw("LP SAME")
    g_frac_ap_primary_18q.Draw("LP SAME")
    g_frac_ap_secondary_18q.Draw("LP SAME")
    
    leg_frac.AddEntry(g_frac_p_primary_18q, "Proton Primary", "lp")
    leg_frac.AddEntry(g_frac_p_secondary_18q, "Proton Secondary", "lp")
    leg_frac.AddEntry(g_frac_ap_primary_18q, "Antiproton Primary", "lp")
    leg_frac.AddEntry(g_frac_ap_secondary_18q, "Antiproton Secondary", "lp")
    leg_frac.Draw()
    
    g_frac_p_primary_18q.SetTitle(f"Proton source fraction (Centrality {cent});pT (GeV/c);Fraction")
    g_frac_p_primary_18q.GetYaxis().SetRangeUser(0, 1)
    c3.cd(1).SetGrid()

    # 绘制18r的fraction
    c3.cd(2)
    g_frac_p_primary_18r.SetLineColor(ROOT.kRed)
    g_frac_p_secondary_18r.SetLineColor(ROOT.kRed)
    g_frac_ap_primary_18r.SetLineColor(ROOT.kBlue)
    g_frac_ap_secondary_18r.SetLineColor(ROOT.kBlue)

    g_frac_p_primary_18r.SetMarkerColor(ROOT.kRed)
    g_frac_p_secondary_18r.SetMarkerColor(ROOT.kRed)
    g_frac_ap_primary_18r.SetMarkerColor(ROOT.kBlue)
    g_frac_ap_secondary_18r.SetMarkerColor(ROOT.kBlue)

    g_frac_p_primary_18r.SetMarkerStyle(20)
    g_frac_p_secondary_18r.SetMarkerStyle(24)
    g_frac_ap_primary_18r.SetMarkerStyle(21)
    g_frac_ap_secondary_18r.SetMarkerStyle(25)
    
    g_frac_p_primary_18r.Draw("ALP")
    g_frac_p_secondary_18r.Draw("LP SAME")
    g_frac_ap_primary_18r.Draw("LP SAME")
    g_frac_ap_secondary_18r.Draw("LP SAME")
    
    g_frac_p_primary_18r.SetTitle(f"LHC18r Fraction (Centrality {cent});pT (GeV/c);Fraction")
    g_frac_p_primary_18r.GetYaxis().SetRangeUser(0, 1)
    c3.cd(2).SetGrid()

    c3.SaveAs(f"fraction_distribution_cent{cent}.pdf")

    # --- 新的直方图 ---
    nbins = h_pt_proton_18q.GetNbinsX()
    xlow = h_pt_proton_18q.GetXaxis().GetXmin()
    xup = h_pt_proton_18q.GetXaxis().GetXmax()

    # 创建新的直方图
    hists = {
        'proton_18q': {
            'primary': ROOT.TH1D(f"h_pt_primary_proton_18q_cent{cent}", f"Primary pT Proton 18q (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup),
            'secondary': ROOT.TH1D(f"h_pt_secondary_proton_18q_cent{cent}", f"Secondary pT Proton 18q (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup)
        },
        'proton_18r': {
            'primary': ROOT.TH1D(f"h_pt_primary_proton_18r_cent{cent}", f"Primary pT Proton 18r (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup),
            'secondary': ROOT.TH1D(f"h_pt_secondary_proton_18r_cent{cent}", f"Secondary pT Proton 18r (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup)
        },
        'antiproton_18q': {
            'primary': ROOT.TH1D(f"h_pt_primary_antiproton_18q_cent{cent}", f"Primary pT Antiproton 18q (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup),
            'secondary': ROOT.TH1D(f"h_pt_secondary_antiproton_18q_cent{cent}", f"Secondary pT Antiproton 18q (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup)
        },
        'antiproton_18r': {
            'primary': ROOT.TH1D(f"h_pt_primary_antiproton_18r_cent{cent}", f"Primary pT Antiproton 18r (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup),
            'secondary': ROOT.TH1D(f"h_pt_secondary_antiproton_18r_cent{cent}", f"Secondary pT Antiproton 18r (Centrality {cent});pT (GeV/c);Counts", nbins, xlow, xup)
        }
    }

    # --- 填充新的分布 ---
    for i in range(1, nbins+1):
        pt = h_pt_proton_18q.GetBinCenter(i)
        
        # 处理质子和反质子的18q和18r数据
        for period in ['18q', '18r']:
            # 质子
            h_pt = h_pt_proton_18q if period == '18q' else h_pt_proton_18r
            g_primary = g_frac_p_primary_18q if period == '18q' else g_frac_p_primary_18r
            g_secondary = g_frac_p_secondary_18q if period == '18q' else g_frac_p_secondary_18r
            
            cnt = h_pt.GetBinContent(i)
            err = h_pt.GetBinError(i)
            f_primary = g_primary.Eval(pt)
            f_secondary = g_secondary.Eval(pt)
            e_primary = g_primary.GetErrorY(i-1)
            e_secondary = g_secondary.GetErrorY(i-1)
            
            val_primary = cnt * f_primary
            err_primary = TMath.Sqrt((err*f_primary)**2 + (cnt*e_primary)**2)
            val_secondary = cnt * f_secondary
            err_secondary = TMath.Sqrt((err*f_secondary)**2 + (cnt*e_secondary)**2)
            
            hists[f'proton_{period}']['primary'].SetBinContent(i, val_primary)
            hists[f'proton_{period}']['primary'].SetBinError(i, err_primary)
            hists[f'proton_{period}']['secondary'].SetBinContent(i, val_secondary)
            hists[f'proton_{period}']['secondary'].SetBinError(i, err_secondary)
            
            # 反质子
            h_pt = h_pt_antiproton_18q if period == '18q' else h_pt_antiproton_18r
            g_primary = g_frac_ap_primary_18q if period == '18q' else g_frac_ap_primary_18r
            g_secondary = g_frac_ap_secondary_18q if period == '18q' else g_frac_ap_secondary_18r
            
            cnt = h_pt.GetBinContent(i)
            err = h_pt.GetBinError(i)
            f_primary = g_primary.Eval(pt)
            f_secondary = g_secondary.Eval(pt)
            e_primary = g_primary.GetErrorY(i-1)
            e_secondary = g_secondary.GetErrorY(i-1)
            
            val_primary = cnt * f_primary
            err_primary = TMath.Sqrt((err*f_primary)**2 + (cnt*e_primary)**2)
            val_secondary = cnt * f_secondary
            err_secondary = TMath.Sqrt((err*f_secondary)**2 + (cnt*e_secondary)**2)
            
            hists[f'antiproton_{period}']['primary'].SetBinContent(i, val_primary)
            hists[f'antiproton_{period}']['primary'].SetBinError(i, err_primary)
            hists[f'antiproton_{period}']['secondary'].SetBinContent(i, val_secondary)
            hists[f'antiproton_{period}']['secondary'].SetBinError(i, err_secondary)

    # --- 画图 ---
    c1 = ROOT.TCanvas(f"c1_cent{cent}", f"pT Distribution (Centrality {cent})", 1600, 1200)
    c1.Divide(2, 2)

    # 创建图例
    legs = [ROOT.TLegend(0.6, 0.7, 0.9, 0.9) for _ in range(4)]

    # 绘制四个子图
    pads = [
        (1, f"Proton 18q (Centrality {cent})", h_pt_proton_18q, hists['proton_18q']),
        (2, f"Proton 18r (Centrality {cent})", h_pt_proton_18r, hists['proton_18r']),
        (3, f"Antiproton 18q (Centrality {cent})", h_pt_antiproton_18q, hists['antiproton_18q']),
        (4, f"Antiproton 18r (Centrality {cent})", h_pt_antiproton_18r, hists['antiproton_18r'])
    ]

    for pad_num, title, h_total, h_components in pads:
        c1.cd(pad_num)

        h_total.SetTitle("")

        h_total.GetXaxis().SetRangeUser(0.7,5)
        h_components['primary'].GetXaxis().SetRangeUser(0.7,5)
        h_components['secondary'].GetXaxis().SetRangeUser(0.7,5)

        h_total.SetLineColor(ROOT.kBlack)

        h_components['primary'].SetLineColor(ROOT.kRed)
        h_components['secondary'].SetLineColor(ROOT.kBlue)
        
        h_total.SetMarkerStyle(20)
        h_components['primary'].SetMarkerStyle(24)
        h_components['secondary'].SetMarkerStyle(25)
        
        h_total.Draw("HIST")
        h_components['primary'].Draw("HIST SAME")
        h_components['secondary'].Draw("HIST SAME")
        
        leg = legs[pad_num-1]
        leg.AddEntry(h_total, "Total", "l")
        leg.AddEntry(h_components['primary'], "Primary", "l")
        leg.AddEntry(h_components['secondary'], "Secondary", "l")
        leg.Draw()

    c1.SaveAs(f"pt_distribution_all_cent{cent}.pdf")

    # --- 计算比例 ---
    c2 = ROOT.TCanvas(f"c2_cent{cent}", f"Primary/Total Ratio (Centrality {cent})", 1200, 600)
    c2.Divide(2, 2)

    # 创建总结直方图
    h_summary = ROOT.TH1D(f"h_summary_cent{cent}", f"Primary/Total Ratio (Centrality {cent});Category;Ratio", 4, 0.5, 4.5)
    h_summary_secondary = ROOT.TH1D(f"h_summary_secondary_cent{cent}", f"Secondary/Total Ratio (Centrality {cent});Category;Ratio", 4, 0.5, 4.5)
    labels = [f"Proton 18q (Cent {cent})", f"Proton 18r (Cent {cent})", f"Antiproton 18q (Cent {cent})", f"Antiproton 18r (Cent {cent})"]
    for ib, lab in enumerate(labels, start=1):
        h_summary.GetXaxis().SetBinLabel(ib, lab)
        h_summary_secondary.GetXaxis().SetBinLabel(ib, lab)

    # 计算比例和误差
    ratios_primary = []
    ratios_secondary = []
    ratios_primary_pt07 = []
    ratios_secondary_pt07 = []

    for period in ['18q', '18r']:
        for particle in ['proton', 'antiproton']:
            h_total = h_pt_proton_18q if particle == 'proton' and period == '18q' else \
                     h_pt_proton_18r if particle == 'proton' and period == '18r' else \
                     h_pt_antiproton_18q if particle == 'antiproton' and period == '18q' else \
                     h_pt_antiproton_18r
            
            h_primary = hists[f'{particle}_{period}']['primary']
            h_secondary = hists[f'{particle}_{period}']['secondary']
            
            # 计算全pT范围的比例
            ratio_primary, ratio_primary_err = calculate_ratio(h_total, h_primary)
            ratio_secondary, ratio_secondary_err = calculate_ratio(h_total, h_secondary)
            ratios_primary.append((ratio_primary, ratio_primary_err))
            ratios_secondary.append((ratio_secondary, ratio_secondary_err))
            
            # 计算pT > 0.7的比例
            ratio_primary_pt07, ratio_primary_err_pt07 = calculate_ratio(h_total, h_primary, 0.7)
            ratio_secondary_pt07, ratio_secondary_err_pt07 = calculate_ratio(h_total, h_secondary, 0.7)
            ratios_primary_pt07.append((ratio_primary_pt07, ratio_primary_err_pt07))
            ratios_secondary_pt07.append((ratio_secondary_pt07, ratio_secondary_err_pt07))

    # 填充总结直方图
    for ib, (ratio, ratio_err) in enumerate(ratios_primary, start=1):
        h_summary.SetBinContent(ib, ratio)
        h_summary.SetBinError(ib, ratio_err)
    
    for ib, (ratio, ratio_err) in enumerate(ratios_secondary, start=1):
        h_summary_secondary.SetBinContent(ib, ratio)
        h_summary_secondary.SetBinError(ib, ratio_err)

    # 绘制总结图
    c2.cd(1)
    h_summary.SetMarkerStyle(22)
    h_summary.SetMarkerSize(2)
    h_summary.SetMarkerColor(ROOT.kRed)
    h_summary.SetLineColor(ROOT.kRed)
    h_summary.Draw("E1")
    c2.cd(1).SetGrid()
    h_summary.GetYaxis().SetRangeUser(0, 1)
    h_summary.SetTitle(f"Primary/Total Ratio Summary (Centrality {cent})")
    
    c2.cd(2)
    h_summary_secondary.SetMarkerStyle(22)
    h_summary_secondary.SetMarkerSize(2)
    h_summary_secondary.SetMarkerColor(ROOT.kBlue)
    h_summary_secondary.SetLineColor(ROOT.kBlue)
    h_summary_secondary.Draw("E1")
    c2.cd(2).SetGrid()
    h_summary_secondary.GetYaxis().SetRangeUser(0, 1)
    h_summary_secondary.SetTitle(f"Secondary/Total Ratio Summary (Centrality {cent})")
    
    c2.SaveAs(f"ratio_summary_cent{cent}.pdf")

    # --- 输出结果 ---
    fout = ROOT.TFile(f"../results/output_cent{cent}.root", "RECREATE")
    
    # 保存所有直方图
    for period in ['18q', '18r']:
        for particle in ['proton', 'antiproton']:
            for comp in ['primary', 'secondary']:
                hists[f'{particle}_{period}'][comp].Write()
    
    # 保存总结直方图
    h_summary.Write()
    h_summary_secondary.Write()
    
    fout.Close()

    # 打印结果
    print(f"\n中心度 {cent} 的结果:")
    print("\n全pT范围的比例:")
    print("\nPrimary/Total Ratios:")
    for i, (ratio, ratio_err) in enumerate(ratios_primary):
        print(f"{labels[i]}: {ratio:.3f} ± {ratio_err:.3f}")
    
    print("\nSecondary/Total Ratios:")
    for i, (ratio, ratio_err) in enumerate(ratios_secondary):
        print(f"{labels[i]}: {ratio:.3f} ± {ratio_err:.3f}")
    
    print("\npT > 0.7 GeV/c的比例:")
    print("\nPrimary/Total Ratios (pT > 0.7):")
    for i, (ratio, ratio_err) in enumerate(ratios_primary_pt07):
        print(f"{labels[i]}: {ratio:.3f} ± {ratio_err:.3f}")
    
    print("\nSecondary/Total Ratios (pT > 0.7):")
    for i, (ratio, ratio_err) in enumerate(ratios_secondary_pt07):
        print(f"{labels[i]}: {ratio:.3f} ± {ratio_err:.3f}")

    # 返回结果供CSV输出使用
    return {
        'primary_pt07': ratios_primary_pt07,
        'secondary_pt07': ratios_secondary_pt07
    }

def main():
    # 创建结果目录
    if not os.path.exists("../results"):
        os.makedirs("../results")
    
    # 创建CSV文件
    csv_file = "../results/integration_results.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # 写入表头
        writer.writerow(['Centrality', 'Period', 'Particle', 'Primary/Total', 'Primary/Total Error', 
                        'Secondary/Total', 'Secondary/Total Error'])
        
        # 处理所有中心度
        centralities = [0, 1, 2, 3, 4, 5]
        for cent in centralities:
            # 获取该中心度的结果
            results = process_centrality(cent)
            if results is None:
                continue
                
            # 写入结果
            for period in ['18q', '18r']:
                for particle in ['proton', 'antiproton']:
                    idx = 0 if period == '18q' else 2
                    idx += 0 if particle == 'proton' else 1
                    
                    primary_ratio, primary_err = results['primary_pt07'][idx]
                    secondary_ratio, secondary_err = results['secondary_pt07'][idx]
                    
                    writer.writerow([
                        cent,
                        period,
                        particle,
                        f"{primary_ratio:.4f}",
                        f"{primary_err:.4f}",
                        f"{secondary_ratio:.4f}",
                        f"{secondary_err:.4f}"
                    ])
    
    print(f"\n结果已保存到 {csv_file}")

if __name__ == "__main__":
    main()
