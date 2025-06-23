import ROOT, re, array, os
import numpy as np

def integrate_phi_with_pt(input_file):
    """
    Open a ROOT file, integrate phi distributions weighted by pT distributions,
    and save the results to a new ROOT file.
    """
    # 打开输入文件
    f = ROOT.TFile.Open(input_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {input_file}")

    # 从文件名中提取标签
    tag_match = re.search(r'(18[qr])', input_file)
    tag = tag_match.group(1) if tag_match else "unknown"
    
    # 从make_nua.py中获取的粒子列表和中心度边界
    particle_list = (
        "poshadron","neghadron","proton","antiproton","lambda","antilambda"
    )
    cent_edges = [0,7]  # 7 个 centrality 区间
    
    # 输出文件名
    output_file = f"intg_check_{tag}.root"
    fout = ROOT.TFile(output_file, "RECREATE")
    
    # 中心度区间数量
    n_cent = len(cent_edges) - 1
    
    # 创建输出目录
    intg_dir = fout.mkdir("integrated_phi_distributions")
    intg_dir.cd()
    
    # 创建NUA检查目录（应用NUA权重后的phi分布）
    nua_check_dir = fout.mkdir("nua_corrected_phi")
    nua_check_dir.cd()
    
    # 对每个粒子和每个中心度进行积分
    for p in particle_list:
        print(f"Processing {p} with {n_cent} centrality bins")
        
        for ic in range(n_cent):
            # 获取pT分布
            h_pt = f.Get(f"pt_distributions/pt_dist_{p}_{tag}_cent{ic}")
            if not h_pt:
                print(f"WARN: Missing pT distribution for {p}, centrality {ic}")
                continue
            
            # 获取phi分布
            h_phi_pt = f.Get(f"phi_distributions/phi_dist_pt_{p}_{tag}_cent{ic}")
            if not h_phi_pt:
                print(f"WARN: Missing phi distribution for {p}, centrality {ic}")
                continue
                
            # 获取NUA修正
            h_nua = f.Get(f"nua_corrections/nua_pt_{p}_{tag}_cent{ic}")
            if not h_nua:
                print(f"WARN: Missing NUA correction for {p}, centrality {ic}")
                continue
            
            # 获取phi轴上的bin数
            n_phi = h_phi_pt.GetNbinsY()
            phi_min = h_phi_pt.GetYaxis().GetXmin()
            phi_max = h_phi_pt.GetYaxis().GetXmax()
            
            # 创建积分后的phi分布直方图
            h_phi_intg = ROOT.TH1D(
                f"phi_intg_{p}_{tag}_cent{ic}", 
                f"Integrated #phi Distribution: {p}, Centrality {cent_edges[ic]}-{cent_edges[ic+1]}%",
                n_phi, phi_min, phi_max
            )
            h_phi_intg.SetDirectory(0)
            h_phi_intg.GetXaxis().SetTitle("#phi (rad)")
            h_phi_intg.GetYaxis().SetTitle("Weighted Counts")
            
            # 积分过程：遍历每个phi bin
            for iphi in range(1, n_phi+1):
                weighted_sum = 0.0
                pt_sum = 0.0
                
                # 对每个pT bin进行加权求和
                for ipt in range(1, h_pt.GetNbinsX()+1):
                    pt_val = h_pt.GetBinContent(ipt)
                    phi_val = h_phi_pt.GetBinContent(ipt, iphi)
                    
                    # 加权求和：pT分布 * phi分布
                    weighted_sum += pt_val * phi_val
                    pt_sum += pt_val
                
                # 填入积分后的直方图
                if pt_sum > 0:
                    h_phi_intg.SetBinContent(iphi, weighted_sum / pt_sum)
                else:
                    h_phi_intg.SetBinContent(iphi, 0)
            
            # 归一化：使得平均值为1
            # 这确保了输出的分布是对平均值的偏差
            mean = 0.0
            count = 0
            for i in range(1, n_phi+1):
                content = h_phi_intg.GetBinContent(i)
                if content > 0:
                    mean += content
                    count += 1
            
            if count > 0:
                mean /= count
                if mean > 0:
                    h_phi_intg.Scale(1.0 / mean)
            
            # 保存到输出文件
            intg_dir.cd()
            h_phi_intg.Write()
            
            # 创建应用NUA权重后的phi分布
            h_nua_corrected = ROOT.TH1D(
                f"nua_corrected_phi_{p}_{tag}_cent{ic}", 
                f"NUA Corrected #phi Distribution: {p}, Centrality {cent_edges[ic]}-{cent_edges[ic+1]}%",
                n_phi, phi_min, phi_max
            )
            h_nua_corrected.SetDirectory(0)
            h_nua_corrected.GetXaxis().SetTitle("#phi (rad)")
            h_nua_corrected.GetYaxis().SetTitle("Counts (NUA corrected)")
            
            # 对每个phi bin应用NUA权重
            for iphi in range(1, n_phi+1):
                weighted_sum = 0.0
                pt_sum = 0.0
                
                # 对每个pT bin进行加权求和，并应用NUA修正
                for ipt in range(1, h_pt.GetNbinsX()+1):
                    pt_val = h_pt.GetBinContent(ipt)
                    phi_val = h_phi_pt.GetBinContent(ipt, iphi)
                    nua_val = h_nua.GetBinContent(ipt, iphi)
                    
                    # 加权求和：pT分布 * phi分布 * NUA权重
                    corrected_val = phi_val * nua_val if nua_val > 0 else 0
                    weighted_sum += pt_val * corrected_val
                    pt_sum += pt_val
                
                # 填入NUA修正后的直方图
                if pt_sum > 0:
                    h_nua_corrected.SetBinContent(iphi, weighted_sum / pt_sum)
                else:
                    h_nua_corrected.SetBinContent(iphi, 0)
            
            # 归一化：使得平均值为1
            mean_nua = 0.0
            count_nua = 0
            for i in range(1, n_phi+1):
                content = h_nua_corrected.GetBinContent(i)
                if content > 0:
                    mean_nua += content
                    count_nua += 1
            
            if count_nua > 0:
                mean_nua /= count_nua
                if mean_nua > 0:
                    h_nua_corrected.Scale(1.0 / mean_nua)
            
            # 保存NUA修正后的分布
            nua_check_dir.cd()
            h_nua_corrected.Write()
    
    fout.Close()
    f.Close()
    
    print(f"Done → {output_file}")
    return output_file

if __name__ == "__main__":
    # 默认输入文件
    input_file = "nua_18q.root"
    
    # 检查命令行参数
    import sys
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    
    try:
        output_file = integrate_phi_with_pt(input_file)
        print(f"Successfully integrated phi distributions with pT weights and saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")