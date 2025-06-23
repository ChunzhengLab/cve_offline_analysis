import ROOT
from array import array
import math
import sys

# =================== 统一管理每个粒子的bin边界 ===================
pt_bins_dict = {
    "proton":      [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 5.0],
    "antiproton":  [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 5.0],
    "poshadron":   [0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    "neghadron":   [0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    "pospion":     [0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    "negpion":     [0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    "lambda":      [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0],
    "antilambda":  [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                    2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0],
}

def get_rebin_edges(particle):
    return pt_bins_dict.get(particle, None)

def rebin_hist(hist, new_edges):
    arr = array('d', new_edges)
    name = hist.GetName() + "_rebin"
    nbins = len(new_edges) - 1
    return hist.Rebin(nbins, name, arr)

def find_inconsistent_bins(hnum, hden):
    inconsistent_bins = []
    for ibin in range(1, hnum.GetNbinsX() + 1):
        passed = hnum.GetBinContent(ibin)
        total = hden.GetBinContent(ibin)
        if passed > total + 1e-6:  # allow a tiny numerical tolerance
            inconsistent_bins.append(ibin)
    return inconsistent_bins

def pt_ranges_from_bins(hist, bins):
    edges = [hist.GetXaxis().GetBinLowEdge(i) for i in bins]
    edges2 = [hist.GetXaxis().GetBinUpEdge(bins[0])] if len(bins)==1 else [hist.GetXaxis().GetBinUpEdge(i) for i in bins]
    return [f"[{edges[i-1]:.3f}, {edges2[i-1]:.3f})" for i in range(1, len(bins)+1)]

def create_tgrapherrors(h_rc, h_mc):
    n_bins = h_rc.GetNbinsX()
    x = array('d', [0] * n_bins)
    y = array('d', [0] * n_bins)
    ex = array('d', [0] * n_bins)
    ey = array('d', [0] * n_bins)

    for i in range(1, n_bins + 1):
        bin_center = h_rc.GetBinCenter(i)
        bin_width = h_rc.GetBinWidth(i)

        passed = h_rc.GetBinContent(i)
        total = h_mc.GetBinContent(i)

        # Calculate efficiency and error
        eff = passed / total if total > 0 else 0
        # 使用二项分布误差计算
        err = math.sqrt(eff * (1 - eff) / total) if total > 0 and 0 < eff < 1 else 0

        x[i-1] = bin_center
        y[i-1] = eff
        ex[i-1] = bin_width / 2
        ey[i-1] = err

    graph = ROOT.TGraphErrors(n_bins, x, y, ex, ey)
    return graph

def main():
    ROOT.TH1.AddDirectory(False)  # Prevent ROOT from automatically managing objects
    file_mc_18q = ROOT.TFile.Open("../data/mc_18q.root")
    file_mc_18r = ROOT.TFile.Open("../data/mc_18r.root")
    out = ROOT.TFile("all_eff_pt.root", "RECREATE")

    # Create TList to store histograms
    fListNUE = ROOT.TList()
    fListNUE.SetName("fListNUE")
    fListNUE.SetOwner(False)  # Do not take ownership of objects

    # Store information for later conversion to TH1D
    eff_info = []

    particles = [
        "poshadron", "pospion", "neghadron", "negpion",
        "proton", "antiproton", "lambda", "antilambda"
    ]
    datasets = ["18q", "18r"]
    kinds = ["rc", "rc_real"]

    for dataset, f in zip(datasets, [file_mc_18q, file_mc_18r]):
        dname = "MyTask"
        l = f.Get(f"{dname}/ResultsList")
        for particle in particles:
            for kind in kinds:
                hname_mc = f"h_pt_{particle}_mc"
                hname_rc = f"h_pt_{particle}_{kind}"
                print(f"Trying {particle} {dataset} {kind}...")

                h_mc = l.FindObject(hname_mc)
                h_rc = l.FindObject(hname_rc)
                if not h_mc or not h_rc:
                    print(f"Skip {particle} {dataset} {kind}: missing hist")
                    continue

                # rebin
                custom_bins = get_rebin_edges(particle)
                if custom_bins:
                    h_mc_reb = rebin_hist(h_mc, custom_bins)
                    h_rc_reb = rebin_hist(h_rc, custom_bins)
                else:
                    h_mc_reb = h_mc.Clone()
                    h_rc_reb = h_rc.Clone()

                inconsistent = find_inconsistent_bins(h_rc_reb, h_mc_reb)
                if inconsistent:
                    ptbins = []
                    for ibin in inconsistent:
                        low = h_mc_reb.GetXaxis().GetBinLowEdge(ibin)
                        high = h_mc_reb.GetXaxis().GetBinUpEdge(ibin)
                        ptbins.append(f"[{low:.3f}, {high:.3f})")
                    print(f"    -- 某些bin: passed > total, 比如第 {inconsistent} bin, pT区间: {ptbins}（仅显示前5个）")
                    print(f"    -> Skip {particle} {dataset} {kind}: bin content inconsistent or missing hist")
                    continue

                out.cd()
                # 检查是否为rc类型且不是hadron粒子
                is_hadron = "hadron" in particle  # 检查particle名称中是否包含"hadron"

                if kind == "rc" and not is_hadron:
                    # 对于rc类型且非hadron粒子，使用TGraphErrors
                    graph = create_tgrapherrors(h_rc_reb, h_mc_reb)
                    graph.SetName(f"graph_ratio_pt_{particle}_{dataset}")
                    graph.SetTitle(f"Efficiency for {particle} ({dataset})")
                    fListNUE.Add(graph)
                    print(f"    -> Added graph_ratio_pt_{particle}_{dataset} as TGraphErrors to fListNUE")
                else:  # rc_real或hadron粒子
                    # 对于rc_real或hadron粒子仍然使用TEfficiency
                    eff = ROOT.TEfficiency(h_rc_reb, h_mc_reb)
                    if kind == "rc_real":
                        eff.SetName(f"eff_pt_{particle}_{dataset}")
                    else:
                        eff.SetName(f"eff_pt_{particle}_{dataset}")
                    fListNUE.Add(eff)
                    # Store efficiency name and particle/dataset/kind information for later conversion
                    eff_name = eff.GetName()
                    eff_info.append({
                        "name": eff_name,
                        "particle": particle,
                        "dataset": dataset,
                        "kind": kind
                    })
                    print(f"    -> Added eff_pt_{particle}{'_real' if kind=='rc_real' else ''}_{dataset} as TEfficiency to fListNUE")

    # Write the TList to the output file
    out.cd()
    fListNUE.Write("fListNUE", ROOT.TObject.kSingleKey)
    # Clear the list without deleting objects before closing the file
    fListNUE.SetOwner(False)
    fListNUE.Clear()
    out.Close()
    print("All done! Saved fListNUE to all_eff_pt.root")

    # Convert TEfficiency objects to TH1D and save to a new file
    print("\nConverting TEfficiency objects to TH1D...")

    # Re-read the input files to process the TEfficiency objects
    ROOT.TH1.AddDirectory(False)  # Prevent ROOT from automatically managing objects
    file_mc_18q = ROOT.TFile.Open("mc_18q.root")
    file_mc_18r = ROOT.TFile.Open("mc_18r.root")
    calib_file = ROOT.TFile("eff_pt_calib.root", "RECREATE")

    # Also open the original efficiency file to use as reference
    orig_file = ROOT.TFile("all_eff_pt.root", "READ")

    # Create TList for calibration file
    fListNUE_calib = ROOT.TList()
    fListNUE_calib.SetName("fListNUE")
    fListNUE_calib.SetOwner(False)  # Do not take ownership of objects

    for info in eff_info:
        particle = info["particle"]
        dataset = info["dataset"]
        kind = info["kind"]
        eff_name = info["name"]

        # Get data from the correct source file
        f = file_mc_18q if dataset == "18q" else file_mc_18r
        dname = "MyTask"
        l = f.Get(f"{dname}/ResultsList")

        # Get the original histograms
        hname_mc = f"h_pt_{particle}_mc"
        hname_rc = f"h_pt_{particle}_{kind}"

        h_mc = l.FindObject(hname_mc)
        h_rc = l.FindObject(hname_rc)

        if not h_mc or not h_rc:
            print(f"    -> Skip converting {eff_name}: missing original histograms")
            continue

        # Apply rebinning as in the original process
        custom_bins = get_rebin_edges(particle)
        if custom_bins:
            h_mc_reb = rebin_hist(h_mc, custom_bins)
            h_rc_reb = rebin_hist(h_rc, custom_bins)
        else:
            h_mc_reb = h_mc.Clone()
            h_rc_reb = h_rc.Clone()

        # Get the original TEfficiency from the file as reference
        orig_eff = orig_file.Get(eff_name)

        # Create a new TH1D with the same binning
        h_eff = h_mc_reb.Clone(eff_name + "_hist")
        h_eff.Reset()  # Clear all bin contents
        h_eff.SetTitle(eff_name + " (Histogram)")

        # Fill the histogram with efficiency values
        for i in range(1, h_eff.GetNbinsX() + 1):
            passed = h_rc_reb.GetBinContent(i)
            total = h_mc_reb.GetBinContent(i)

            # Calculate efficiency
            if total > 0:
                eff_value = passed / total
                # Calculate error using binomial distribution
                eff_error = math.sqrt(eff_value * (1 - eff_value) / total) if 0 < eff_value < 1 else 0

                h_eff.SetBinContent(i, eff_value)
                h_eff.SetBinError(i, eff_error)

        # Add the TH1D to the TList
        fListNUE_calib.Add(h_eff)
        print(f"    -> Converted {eff_name} to TH1D and added to fListNUE")

    # Write the TList to the calibration file
    calib_file.cd()
    fListNUE_calib.Write("fListNUE", ROOT.TObject.kSingleKey)

    # Clear the list without deleting objects before closing the file
    fListNUE_calib.SetOwner(False)
    fListNUE_calib.Clear()

    file_mc_18q.Close()
    file_mc_18r.Close()
    orig_file.Close()
    calib_file.Close()
    print("All done! Saved fListNUE to eff_pt_calib.root")

if __name__ == "__main__":
    main()
