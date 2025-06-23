import ROOT, re, array, os

# -------- 配置 --------
particle_list = (
    "poshadron","neghadron","proton","antiproton","lambda","antilambda"
)
cent_edges = [0,7]  # 7 个 centrality 区间
# phi bin 合并因子: 每 phi_rebin_factor 个 bin 合并为一个
phi_rebin_factor = 1

pt_ranges = {
    "poshadron":  (0.2,5.0),
    "neghadron":  (0.2,5.0),
    "proton":     (0.7,1.0),
    "antiproton": (0.7,1.0),
    "lambda":     (1.0,5.0,10.0),
    "antilambda": (1.0,5.0,10.0)
}
#
#
pt_ranges = {
    "poshadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "neghadron":  (0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.3,2.5,2.7,2.9,3.1,3.5,4.0,5.0),
    "proton":     (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "antiproton": (0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.9,2.2,2.6,3.1,3.5,4.0,5.0),
    "lambda":     (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0),
    "antilambda": (1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.7,2.9,3.1,3.5,4.0,5.0,10.0)
}


def get_run_tag(fn):
    m = re.search(r'(18[qr])', fn)
    return m.group(1) if m else "unknown"

def extract_nua2d(input_file, tag):
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {input_file}")
    lst = f.Get("default/ResultsList_default")
    if not lst:
        raise RuntimeError("Cannot find default/ResultsList_default")

    out = {}
    phi_dist = {}
    pt_dist = {}
    for p in particle_list:
        h3 = lst.FindObject(f"h3_pt_phi_{p}")
        if not h3:
            print(f"WARN: missing h3_pt_phi_{p}")
            continue

        # 原始 phi bin 数
        n_phi_orig = h3.GetZaxis().GetNbins()
        if n_phi_orig % phi_rebin_factor != 0:
            print(f"WARN: original phi bins {n_phi_orig} not divisible by rebin factor {phi_rebin_factor}, skipping {p}")
            continue
        # 合并后 phi bin 数
        n_phi = n_phi_orig // phi_rebin_factor

        # 构造新的 phi 边界：按照 rebin factor 取边界
        phi_edges_orig = [h3.GetZaxis().GetBinLowEdge(i+1) for i in range(n_phi_orig)]
        phi_edges_orig.append(h3.GetZaxis().GetBinUpEdge(n_phi_orig))
        phi_edges = array.array('d', [phi_edges_orig[phi_rebin_factor*i] for i in range(n_phi+1)])

        pt_edges = pt_ranges[p]
        n_pt = len(pt_edges) - 1
        n_cent = len(cent_edges) - 1

        h2s = []
        phi_dist_h2s = []
        for ic in range(1, n_cent+1):
            h2 = ROOT.TH2F(f"nua_pt_{p}_{tag}_cent{ic-1}", f"NUA Correction: {p}, Centrality {cent_edges[ic-1]}-{cent_edges[ic]}%",
                           n_pt, array.array('d',pt_edges),
                           n_phi, phi_edges)
            h2.SetDirectory(0)
            h2.GetXaxis().SetTitle("p_{T} (GeV/c)")
            h2.GetYaxis().SetTitle("#phi (rad)")

            # 创建phi分布直方图
            h2_phi_dist = ROOT.TH2F(f"phi_dist_pt_{p}_{tag}_cent{ic-1}", f"Normalized #phi Distribution (value/mean): {p}, Centrality {cent_edges[ic-1]}-{cent_edges[ic]}%",
                                    n_pt, array.array('d',pt_edges),
                                    n_phi, phi_edges)
            h2_phi_dist.SetDirectory(0)
            h2_phi_dist.GetXaxis().SetTitle("p_{T} (GeV/c)")
            h2_phi_dist.GetYaxis().SetTitle("#phi (rad)")

            # 对每个 pt bin 计算 phi 合并后的 nua
            for ip in range(n_pt):
                pt_lo, pt_hi = pt_edges[ip], pt_edges[ip+1]
                biny1 = h3.GetYaxis().FindBin(pt_lo)
                biny2 = h3.GetYaxis().FindBin(pt_hi - 1e-9)

                # 先算原始每个 phi-bin 的总和
                P_orig = [0.0]*n_phi_orig
                for iph in range(1, n_phi_orig+1):
                    s = 0.0
                    for jpt in range(biny1, biny2+1):
                        s += h3.GetBinContent(ic, jpt, iph)
                    P_orig[iph-1] = s

                # 按照 rebin factor 合并 bin → P_grp
                P_grp = []
                for j in range(n_phi):
                    bin_sum = 0
                    for k in range(phi_rebin_factor):
                        bin_sum += P_orig[phi_rebin_factor*j + k]
                    P_grp.append(bin_sum)

                # 计算均值 M_grp
                M_grp = sum(P_grp)/len(P_grp) if P_grp else 0.0

                # 填入 TH2D：nua = M_grp / P_grp[j]
                for j in range(n_phi):
                    val = P_grp[j]
                    nua = M_grp/val if val>0 else 0.0
                    h2.SetBinContent(ip+1, j+1, nua)

                    # 填入归一化的phi分布直方图 (normalized to mean=1, 这样分布均匀时所有bin都是1.0)
                    normalized_val = val/M_grp if M_grp>0 else 0.0
                    h2_phi_dist.SetBinContent(ip+1, j+1, normalized_val)

            h2s.append(h2)
            phi_dist_h2s.append(h2_phi_dist)

        # 创建pt分布直方图
        pt_dist_h1s = []
        for ic in range(1, n_cent+1):
            h1 = ROOT.TH1D(f"pt_dist_{p}_{tag}_cent{ic-1}", f"p_{{T}} Distribution: {p}, Centrality {cent_edges[ic-1]}-{cent_edges[ic]}%",
                          n_pt, array.array('d', pt_edges))
            h1.SetDirectory(0)
            h1.GetXaxis().SetTitle("p_{T} (GeV/c)")
            h1.GetYaxis().SetTitle("dN/dp_{T}")

            # 计算每个pt bin的计数
            for ip in range(n_pt):
                pt_lo, pt_hi = pt_edges[ip], pt_edges[ip+1]
                biny1 = h3.GetYaxis().FindBin(pt_lo)
                biny2 = h3.GetYaxis().FindBin(pt_hi - 1e-9)

                # 对所有phi bin求和得到pT分布
                s = 0.0
                for jpt in range(biny1, biny2+1):
                    for iph in range(1, n_phi_orig+1):
                        s += h3.GetBinContent(ic, jpt, iph)

                # pT bin宽度归一化
                bin_width = pt_hi - pt_lo
                if bin_width > 0:
                    s = s / bin_width

                h1.SetBinContent(ip+1, s)

            pt_dist_h1s.append(h1)

        out[p] = h2s
        phi_dist[p] = phi_dist_h2s
        pt_dist[p] = pt_dist_h1s

    return out, phi_dist, pt_dist

def write_nua2d(out, phi_dist, pt_dist, infile, outfile):
    tag = get_run_tag(infile)
    fout = ROOT.TFile(outfile, "RECREATE")

    # 写入NUA直方图
    nua_dir = fout.mkdir("nua_corrections")
    nua_dir.cd()
    for p, h2s in out.items():
        for ic, h2 in enumerate(h2s):
            h2.Write()

    # 写入phi分布直方图
    phi_dir = fout.mkdir("phi_distributions")
    phi_dir.cd()
    for p, h2s in phi_dist.items():
        for ic, h2 in enumerate(h2s):
            h2.Write()

    # 写入pt分布直方图
    pt_dir = fout.mkdir("pt_distributions")
    pt_dir.cd()
    for p, h1s in pt_dist.items():
        for ic, h1 in enumerate(h1s):
            h1.Write()

    fout.Close()

if __name__=="__main__":
    infile  = "data_18r.root"
    outfile = "nua_18r.root"
    tag = get_run_tag(infile)
    res, phi_dist, pt_dist = extract_nua2d(infile, tag)
    write_nua2d(res, phi_dist, pt_dist, infile, outfile)
    print("Done →", outfile)
    print(f"Generated NUA corrections, phi distributions, and pT distributions for {len(res)} particles")

    infile = "data_18q.root"
    outfile = "nua_18q.root"
    tag = get_run_tag(infile)
    res, phi_dist, pt_dist = extract_nua2d(infile, tag)
    write_nua2d(res, phi_dist, pt_dist, infile, outfile)
    print("Done →", outfile)
    print(f"Generated NUA corrections, phi distributions, and pT distributions for {len(res)} particles")
