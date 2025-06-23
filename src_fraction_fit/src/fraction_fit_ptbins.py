import ROOT
import os
import numpy as np
import array
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter1d

def moving_average(y, window_size=8):
    return np.convolve(y, np.ones(window_size)/window_size, mode='same')

def gaussian_smooth(y, sigma=1):
    return gaussian_filter1d(y, sigma)

def get_bin_centers_and_counts(h, cut_at_last_nonzero=True):
    nbins = h.GetNbinsX()
    centers = np.array([h.GetXaxis().GetBinCenter(i+1) for i in range(nbins)])
    counts = np.array([h.GetBinContent(i+1) for i in range(nbins)])

    if cut_at_last_nonzero:
        nonzero_indices = np.where(counts > 0)[0]
        if len(nonzero_indices) > 1:
            cut_idx = nonzero_indices[-2] + 1
            centers = centers[:cut_idx]
            counts = counts[:cut_idx]
        # 如果全是零或者只有一个非零点就不切割

    return centers, counts

from scipy.optimize import minimize
import numpy as np

def template_fraction_fit_bootstrap(
    data_dca, data_counts,
    mc0_dca, mc0_counts,
    mc1_dca, mc1_counts,
    n_bootstrap=500,
    bin_weights=None
):
    """
    用最小二乘拟合 data = f0*mc0 + f1*mc1 (f0+f1=1)，并用bootstrap估算误差。
    返回 f0, f1, f0_err, f1_err, chi2_val, ndf
    """

    def single_fit(d_counts):
        n_bins = len(d_counts)
        if bin_weights is None:
            weights = np.ones(n_bins)
        else:
            weights = bin_weights
        def chi2(frac):
            f0, f1 = frac
            model = f0*mc0_counts + f1*mc1_counts
            err = np.sqrt(np.maximum(d_counts, 1e-9))
            return np.sum(weights * ((d_counts - model)/err)**2)
        cons = {'type':'eq', 'fun': lambda x: x[0] + x[1] - 1}
        bounds = [(0, 1), (0, 1)]
        x0 = [0.5, 0.5]
        result = minimize(chi2, x0, bounds=bounds, constraints=cons)
        return result.x if result.success else [np.nan, np.nan]

    # 先用真实data拟合一次
    best_frac = single_fit(data_counts)
    f0, f1 = best_frac

    # 计算chi2和ndf
    model = f0 * mc0_counts + f1 * mc1_counts
    err = np.sqrt(np.maximum(data_counts, 1e-9))
    chi2_val = np.sum(((data_counts - model)/err) ** 2)
    ndf = len(data_counts) - 2  # 参数个数为2

    # bootstrap采样，估算误差
    all_frac = []
    n_bins = len(data_counts)
    total_entries = np.sum(data_counts)
    for _ in range(n_bootstrap):
        # 对data histogram做Poisson采样
        resample = np.random.poisson(lam=data_counts * total_entries)
        # 避免全0，保证归一化
        if np.sum(resample) > 0:
            resample = resample / (np.sum(resample) + 1e-20)
        else:
            resample = np.zeros_like(resample)
        frac = single_fit(resample)
        all_frac.append(frac)
    all_frac = np.array(all_frac)
    f0_err = np.nanstd(all_frac[:,0])
    f1_err = np.nanstd(all_frac[:,1])

    print(f"Chi-square value: {chi2_val}")
    print(f"Degrees of freedom: {ndf}")
    print(f"Bootstrap error: f0_err={f0_err}, f1_err={f1_err}")

    return f0, f1, f0_err, f1_err, chi2_val, ndf

def template_fraction_fit_local(
    data_dca, data_counts,
    mc0_dca, mc0_counts,
    mc1_dca, mc1_counts,
    bin_weights=None
):
    n_bins = len(data_counts)
    if bin_weights is None:
        bin_weights = np.ones(n_bins)

    def chi2(frac):
        f0, f1 = frac
        model = f0*mc0_counts + f1*mc1_counts
        err = np.sqrt(np.maximum(data_counts, 1e-9))
        return np.sum(bin_weights * ((data_counts - model)/err)**2)
    cons = {'type':'eq', 'fun': lambda x: x[0] + x[1] - 1}
    bounds = [(0, 1), (0, 1)]
    x0 = [0.5, 0.5]
    # trust-constr 可以返回 hess_inv
    result = minimize(chi2, x0, bounds=bounds, constraints=cons, method="trust-constr")
    f0, f1 = result.x

    # 计算chi2和ndf
    model = f0 * mc0_counts + f1 * mc1_counts
    err = np.sqrt(np.maximum(data_counts, 1e-9))
    chi2_val = np.sum(((data_counts - model)/err) ** 2)
    ndf = len(data_counts) - 2

    # Hessian 估算误差
    if hasattr(result, "hess_inv") and result.hess_inv is not None:
        try:
            cov = result.hess_inv.todense()
            f0_err = np.sqrt(abs(cov[0,0]))
            f1_err = np.sqrt(abs(cov[1,1]))
        except Exception:
            f0_err = f1_err = np.nan
    else:
        f0_err = f1_err = np.nan

    return f0, f1, f0_err, f1_err, chi2_val, ndf

def fraction_fit_ptbins(dataset="18q", particle="proton", output=None):
    ROOT.gStyle.SetOptStat(0)

    data_file = f"../data/data_{dataset}.root"
    mc_file = f"../data/mc_{dataset}.root"
    out_file = output if output else f"fraction_{particle}_{dataset}.root"

    pt_bins = [0.4, 0.6, 0.7, 1, 1.4, 1.9, 2.5, 3.2, 4.0, 5.0]
    pt_ranges = list(zip(pt_bins[:-1], pt_bins[1:]))
    print(f"{pt_ranges}")

    # 创建plots目录
    plots_dir = "../plots"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    fMC = ROOT.TFile.Open(mc_file)
    fout = ROOT.TFile(out_file, "RECREATE")

    fData = ROOT.TFile.Open(data_file)
    dirData = fData.Get("default")
    listData = dirData.Get("ResultsList_default")
    h3_data = listData.FindObject(f"h3_pt_dcaXY_{particle}")
    if h3_data:
        print(f"Found TH3D data {particle}")
    h3_data.SetName(f"h3_data_{particle}")

    dirMC = fMC.Get("default")
    listMC = dirMC.Get("ResultsList_default")
    h3_mc_origin = listMC.FindObject(f"h3_pt_dcaXY_origin_rc_{particle}")
    h3_mc_lambda = listMC.FindObject(f"h3_pt_dcaXY_lambda_rc_{particle}")
    h3_mc_material = listMC.FindObject(f"h3_pt_dcaXY_material_rc_{particle}")
    h3_mc_other = listMC.FindObject(f"h3_pt_dcaXY_other_rc_{particle}")
    if h3_mc_origin and h3_mc_lambda and h3_mc_material and h3_mc_other:
        print(f"Found TH3D mc {particle}")

    h3_mc_primary = h3_mc_origin.Clone("h3_mc_primary")
    h3_mc_secondary = h3_mc_lambda.Clone("h3_mc_secondary")
    h3_mc_secondary.Add(h3_mc_material)
    h3_mc_secondary.Add(h3_mc_other)

    h3_mc_primary = h3_mc_primary.RebinX(7)
    h3_mc_secondary = h3_mc_secondary.RebinX(7)
    h3_mc_primary.Draw()

    for iCentBin in range(1, h3_mc_primary.GetNbinsX()+1):
        # 为每个中心度创建TGraphErrors
        pt_centers = []
        primary_fractions = []
        primary_errors = []
        secondary_fractions = []
        secondary_errors = []

        # 获取实际要处理的pT区间数量
        n_pt_bins = len(pt_ranges)

        # 创建子图的PDF
        fig = plt.figure(figsize=(15, 10))
        
        # 创建主GridSpec，用于控制三组图之间的间距
        outer_gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.6)
        
        # 为每组创建内部的GridSpec，用于控制主图和ratio图
        inner_gs = []
        for i in range(3):
            inner_gs.append(gridspec.GridSpecFromSubplotSpec(2, 3, 
                                                           subplot_spec=outer_gs[i],
                                                           height_ratios=[3, 1],
                                                           hspace=0.0))
        
        # 先进行所有pt区间的拟合
        for i, (pt_low, pt_high) in enumerate(pt_ranges):
            pt_low_bin = h3_mc_primary.GetYaxis().FindBin(pt_low + 1.e-6)
            pt_high_bin = h3_mc_primary.GetYaxis().FindBin(pt_high - 1.e-6)

            h_mc_dcaXY_primary = h3_mc_primary.ProjectionZ(f"h_mc_dcaXY_primary_{iCentBin}_{pt_low}_{pt_high}", iCentBin, iCentBin, pt_low_bin, pt_high_bin)
            h_mc_dcaXY_secondary = h3_mc_secondary.ProjectionZ(f"h_mc_dcaXY_secondary_{iCentBin}_{pt_low}_{pt_high}", iCentBin, iCentBin, pt_low_bin, pt_high_bin)
            h_data_dcaXY = h3_data.ProjectionZ(f"h_data_dcaXY_{iCentBin}_{pt_low}_{pt_high}", iCentBin, iCentBin, pt_low_bin, pt_high_bin)

            h_mc_dcaXY_primary.Scale(1. / h_mc_dcaXY_primary.Integral())
            h_mc_dcaXY_secondary.Scale(1. / h_mc_dcaXY_secondary.Integral())
            h_data_dcaXY.Scale(1. / h_data_dcaXY.Integral())

            h_mc_dcaXY_primary.Rebin(4)
            h_mc_dcaXY_secondary.Rebin(4)
            h_data_dcaXY.Rebin(4)

            print(f"h_mc_dcaXY_primary bins {h_mc_dcaXY_primary.GetXaxis().GetNbins()}")
            print(f"h_mc_dcaXY_secondary bins {h_mc_dcaXY_secondary.GetXaxis().GetNbins()}")
            print(f"h_data_dcaXY bins {h_data_dcaXY.GetXaxis().GetNbins()}")

            # ---------- 自定义chi2最小化并bootstrap误差 -------------
            data_dca, data_counts = get_bin_centers_and_counts(h_data_dcaXY)
            mc0_dca, mc0_counts = get_bin_centers_and_counts(h_mc_dcaXY_primary)
            mc1_dca, mc1_counts = get_bin_centers_and_counts(h_mc_dcaXY_secondary)

            # 确保所有数组长度一致，以较小的为准
            min_len = min(len(data_dca), len(mc0_dca), len(mc1_dca))
            data_dca = data_dca[:min_len]
            data_counts = data_counts[:min_len]
            mc0_dca = mc0_dca[:min_len]
            mc0_counts = mc0_counts[:min_len]
            mc1_dca = mc1_dca[:min_len]
            mc1_counts = mc1_counts[:min_len]

            f0, f1, f0_err, f1_err, chi2_val, ndf = template_fraction_fit_bootstrap(
                data_dca, data_counts,
                mc0_dca, mc0_counts,
                mc1_dca, mc1_counts
            )
            print(f"Centrality Bin {iCentBin}, pt [{pt_low}, {pt_high}]:")
            print(f"Primary fraction: {f0:.6f} ± {f0_err:.6f}")
            print(f"Secondary fraction: {f1:.6f} ± {f1_err:.6f}")
            print(f"Chi2/NDF: {chi2_val:.2f}/{ndf}")

            # 存储拟合结果用于TGraphErrors
            pt_center = (pt_low + pt_high) / 2
            pt_centers.append(pt_center)
            primary_fractions.append(f0)
            primary_errors.append(f0_err)
            secondary_fractions.append(f1)
            secondary_errors.append(f1_err)

            # 只绘制前9个pt区间的图形
            if i < 9:
                row = i // 3
                col = i % 3
                
                # 创建主图和ratio图
                ax_upper = plt.subplot(inner_gs[row][0, col])
                ax_lower = plt.subplot(inner_gs[row][1, col], sharex=ax_upper)
                
                # 归一化
                data_counts_norm = data_counts / np.sum(data_counts)
                mc0_counts_norm = mc0_counts / np.sum(mc0_counts)
                mc1_counts_norm = mc1_counts / np.sum(mc1_counts)
                
                # 创建拟合结果
                fit_counts_norm = f0 * mc0_counts_norm + f1 * mc1_counts_norm
                
                ax_upper.set_title(f"{dataset} {particle}\npT: {pt_low:.1f}-{pt_high:.1f} GeV/c")
                ax_upper.set_yscale('log')
                ax_upper.set_ylim(1e-5, max(np.max(data_counts_norm),
                                            np.max(mc0_counts_norm),
                                            np.max(mc1_counts_norm)) * 1.5)
                
                # 设置x轴范围为当前数据最大值的1.1倍
                current_x_max = np.max(data_dca) * 1.05
                ax_upper.set_xlim(0, current_x_max)
                
                # 绘制数据和模型
                ax_upper.plot(data_dca, data_counts_norm, 'ko', markersize=4, label='Data')
                ax_upper.plot(mc0_dca, mc0_counts_norm, 'r-', linewidth=2, label='Primary')
                ax_upper.plot(mc1_dca, mc1_counts_norm, 'b-', linewidth=2, label='Secondary')
                ax_upper.plot(mc0_dca, fit_counts_norm, 'm-', linewidth=2, label=f'Fit: {f0:.3f}±{f0_err:.3f}')
                
                # 添加图例
                ax_upper.legend(loc='best', fontsize=8)
                ax_upper.set_ylabel('Counts')
                ax_upper.tick_params(axis='x', labelbottom=False)  # 隐藏x轴标签
                
                ax_lower.set_ylim(0.11, 20)
                ax_lower.set_yscale('log')
                
                # 计算比率
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio_primary = np.divide(mc0_counts_norm, data_counts_norm, out=np.ones_like(mc0_counts_norm), where=data_counts_norm != 0)
                    ratio_secondary = np.divide(mc1_counts_norm, data_counts_norm, out=np.ones_like(mc1_counts_norm), where=data_counts_norm != 0)
                    ratio_fit = np.divide(fit_counts_norm, data_counts_norm, out=np.ones_like(fit_counts_norm), where=data_counts_norm != 0)
                
                # 画一条y=1的参考线
                ax_lower.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7)
                ax_lower.plot(mc0_dca, ratio_primary, 'r-', linewidth=2)
                ax_lower.plot(mc1_dca, ratio_secondary, 'b-', linewidth=2)
                ax_lower.plot(mc0_dca, ratio_fit, 'm-', linewidth=2)
                
                ax_lower.set_xlabel('DCA$_{XY}$ (cm)')
                ax_lower.set_ylabel('Ratio')
                ax_lower.tick_params(axis='y', labelsize=8)
                # 移除上下图之间的空白
                plt.setp(ax_upper.get_xticklabels(), visible=False)

        # 保存图形到PDF
        fig.subplots_adjust(
            top=0.95,      # 顶部边距
            bottom=0.05,   # 底部边距
            left=0.1,      # 左侧边距
            right=0.9,     # 右侧边距
            wspace=0.3     # 列间距
        )
        
        pdf_filename = f"../plots/fraction_fit_{particle}_{dataset}_cent{iCentBin}.pdf"
        plt.savefig(pdf_filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"已保存PDF: {pdf_filename}")

        # 创建并保存TGraphErrors
        if len(pt_centers) > 0:
            # 将列表转换为数组
            x = array.array('d', pt_centers)
            y_pri = array.array('d', primary_fractions)
            y_pri_err = array.array('d', primary_errors)
            y_sec = array.array('d', secondary_fractions)
            y_sec_err = array.array('d', secondary_errors)
            # x轴误差设为0
            x_err = array.array('d', [0.0] * len(pt_centers))

            # 创建TGraphErrors
            g_primary = ROOT.TGraphErrors(len(pt_centers), x, y_pri, x_err, y_pri_err)
            g_primary.SetName(f"g_primary_{particle}_{dataset}_cent{iCentBin}")
            g_primary.SetTitle(f"Primary Fraction (Cent. Bin {iCentBin})")
            g_primary.GetXaxis().SetTitle("p_{T} (GeV/c)")
            g_primary.GetYaxis().SetTitle("Primary Fraction")

            g_secondary = ROOT.TGraphErrors(len(pt_centers), x, y_sec, x_err, y_sec_err)
            g_secondary.SetName(f"g_secondary_{particle}_{dataset}_cent{iCentBin}")
            g_secondary.SetTitle(f"Secondary Fraction (Cent. Bin {iCentBin})")
            g_secondary.GetXaxis().SetTitle("p_{T} (GeV/c)")
            g_secondary.GetYaxis().SetTitle("Secondary Fraction")

            # 保存到输出文件
            fout.cd()
            g_primary.Write()
            g_secondary.Write()

            print(f"已保存TGraphErrors: g_primary_cent{iCentBin}, g_secondary_cent{iCentBin}")

    # 关闭文件
    fout.Close()
    fData.Close()
    fMC.Close()

    print(f"完成所有拟合，结果已保存到: {out_file}")

def main():
    fraction_fit_ptbins(dataset = "18q", particle = "proton")
    fraction_fit_ptbins(dataset = "18q", particle = "antiproton")
    fraction_fit_ptbins(dataset = "18r", particle = "proton")
    fraction_fit_ptbins(dataset = "18r", particle = "antiproton")

if __name__ == "__main__":
    main()
