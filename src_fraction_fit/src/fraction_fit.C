#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TObjArray.h>
#include <TFractionFitter.h>
#include <TGraphErrors.h>
#include <vector>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TLine.h>
#include <string>
#include <fstream>

// ---- source结构体统一管理 ----
struct SourceInfo {
    std::string hist2d_name;   // MC模板的TH2D名字
    std::string parname;       // 拟合参数名
    double init_val;
    double step;
    double vmin, vmax;
    int color;
    std::string particle;      // "proton" 或 "antiproton"
    std::string dataset;       // "18q" 或 "18r"
    size_t cent_bin;          // 中心度bin索引
};

// 中心度bin定义
const std::vector<std::pair<double, double>> cent_bins = {
    {0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 60}, {60, 70}
};

// 获取MC模板名称
std::string get_mc_template_name(const std::string& type, const std::string& particle) {
    if (type == "primary") {
        return Form("h3_pt_dcaXY_origin_rc_real_%s", particle.c_str());
    } else if (type == "secondary") {
        return Form("h3_pt_dcaXY_material_rc_real_%s", particle.c_str());
    } else if (type == "material") {
        return Form("h3_pt_dcaXY_material_rc_real_%s", particle.c_str());
    } else if (type == "lambda") {
        return Form("h3_pt_dcaXY_lambda_rc_real_%s", particle.c_str());
    } else if (type == "other") {
        return Form("h3_pt_dcaXY_other_rc_real_%s", particle.c_str());
    }
    return "";
}

// 画图函数也用SourceInfo信息
void DrawFitCompareHistos(TH1D* h_data, const std::vector<TH1D*>& h_mc, const std::vector<SourceInfo>& sources,
                          TH1* h_fit, int ipt, double ptlow, double pthigh)
{
    // 获取中心度信息
    double cent_low = cent_bins[sources[0].cent_bin].first;
    double cent_high = cent_bins[sources[0].cent_bin].second;

    // 创建副本用于归一化，避免修改原始直方图
    TH1D* h_data_norm = (TH1D*)h_data->Clone("h_data_norm");
    std::vector<TH1D*> h_mc_norm;
    for(size_t i = 0; i < h_mc.size(); ++i) {
        h_mc_norm.push_back((TH1D*)h_mc[i]->Clone(Form("h_mc_norm_%zu", i)));
    }
    TH1D* h_fit_norm = nullptr;
    if(h_fit) h_fit_norm = (TH1D*)h_fit->Clone("h_fit_norm");

    // 归一化副本
    if(h_data_norm->Integral() > 0) h_data_norm->Scale(1.0 / h_data_norm->Integral());
    for(auto& h : h_mc_norm) if(h->Integral() > 0) h->Scale(1.0 / h->Integral());
    if(h_fit_norm && h_fit_norm->Integral() > 0) h_fit_norm->Scale(1.0 / h_fit_norm->Integral());

    // 画布和pad
    TCanvas* c1 = new TCanvas(Form("c_fit_cent%zu_ptbin%d", sources[0].cent_bin, ipt),
                             Form("Cent %.0f-%.0f%%, pT %.1f-%.1f", cent_low, cent_high, ptlow, pthigh),
                             800, 800);

    // 创建上下两个pad
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.1);
    pad1->SetLogy();
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    // 上半区：归一化叠加
    pad1->cd();

    h_data_norm->GetXaxis()->SetRangeUser(0, 0.093);
    h_data_norm->GetYaxis()->SetRangeUser(1e-6, 1);
    h_data_norm->SetLineColor(kBlack);
    h_data_norm->SetLineWidth(2);
    h_data_norm->SetMarkerStyle(20);
    h_data_norm->SetMarkerColor(kBlack);
    h_data_norm->SetTitle(Form("%s %s, Cent %.0f-%.0f%%, pT %.1f-%.1f GeV/c",
                         sources[0].dataset.c_str(),
                         sources[0].particle.c_str(),
                         cent_low, cent_high,
                         ptlow, pthigh));
    h_data_norm->GetYaxis()->SetTitle("Normalized Entries");
    h_data_norm->GetYaxis()->SetTitleSize(0.05);
    h_data_norm->GetXaxis()->SetLabelSize(0);  // 隐藏上图的x轴标签
    h_data_norm->Draw("E");

    TLegend* leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data_norm, "Data", "pl");

    for(size_t i=0; i<h_mc_norm.size(); ++i) {
        h_mc_norm[i]->GetXaxis()->SetRangeUser(0, 0.093);
        h_mc_norm[i]->SetLineColor(sources[i].color);
        h_mc_norm[i]->SetLineWidth(2);
        h_mc_norm[i]->Draw("HIST SAME");
        leg->AddEntry(h_mc_norm[i], sources[i].parname.c_str(), "l");
    }

    if(h_fit_norm){
        h_fit_norm->SetLineColor(kMagenta+1);
        h_fit_norm->SetLineWidth(3);
        h_fit_norm->SetLineStyle(2);
        h_fit_norm->Draw("HIST SAME");
        leg->AddEntry(h_fit_norm, "Template Fit", "l");
    }
    leg->Draw();

    // 下半区：ratio plot
    pad2->cd();

    std::vector<TH1D*> h_ratios;
    for(size_t i=0; i<h_mc.size(); ++i) {
        TString rname = Form("ratio_%s_cent%zu_ptbin%d", sources[i].parname.c_str(), sources[0].cent_bin, ipt);
        TH1D* h_ratio = (TH1D*)h_mc[i]->Clone(rname);
        h_ratio->Divide(h_data);
        h_ratio->SetLineColor(sources[i].color);
        h_ratio->SetLineWidth(2);
        h_ratio->SetStats(0);
        h_ratio->SetTitle("");
        h_ratio->GetYaxis()->SetTitle("MC/Data");
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);
        h_ratio->GetYaxis()->SetLabelSize(0.1);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetXaxis()->SetLabelSize(0.1);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetTitleOffset(1.0);
        h_ratio->GetXaxis()->SetTitle("DCA_{xy} (cm)");
        h_ratio->GetXaxis()->SetRangeUser(0, 0.093);
        h_ratio->SetMaximum(3.0);
        h_ratio->SetMinimum(0.0);
        h_ratios.push_back(h_ratio);
    }

    TH1D* h_ratio_fit = nullptr;
    if(h_fit){
        h_ratio_fit = (TH1D*)h_fit->Clone("ratio_fitsum");
        h_ratio_fit->Divide(h_data);
        h_ratio_fit->SetLineColor(kMagenta+1);
        h_ratio_fit->SetLineWidth(3);
        h_ratio_fit->SetLineStyle(2);
        h_ratio_fit->SetStats(0);
        h_ratio_fit->SetTitle("");
        h_ratio_fit->GetXaxis()->SetRangeUser(0, 0.093);
    }

    // 先画第一个
    bool first = true;
    for(auto h : h_ratios){
        if(first) {
            h->Draw("HIST");
            first = false;
        } else {
            h->Draw("HIST SAME");
        }
    }
    if(h_ratio_fit) h_ratio_fit->Draw("HIST SAME");

    // 水平参考线 y=1
    TLine* line = new TLine(0, 1, 0.093, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+1);
    line->Draw();

    c1->SaveAs(Form("fit_compare_%s_%s_cent%zu_ptbin%d.pdf",
                    sources[0].dataset.c_str(),
                    sources[0].particle.c_str(),
                    sources[0].cent_bin,
                    ipt));
    delete c1;
    delete h_data_norm;
    for(auto h : h_mc_norm) delete h;
    if(h_fit_norm) delete h_fit_norm;
    for(auto h : h_ratios) delete h;
    if(h_ratio_fit) delete h_ratio_fit;
}


// 主体
void fraction_fit_ptbins(
    const char* dataset = "18q",      // "18q" 或 "18r"
    const char* particle = "antiproton",   // "proton" 或 "antiproton"
    const char* output = nullptr)      // 如果为nullptr，则自动生成文件名
{
    std::cout << "Starting fraction_fit_ptbins..." << std::endl;
    std::cout << "Parameters: dataset=" << dataset << ", particle=" << particle << std::endl;

    gStyle->SetOptStat(0);
    // 自动生成文件名
    std::string data_file = Form("../data/data_%s.root", dataset);
    std::string mc_file = Form("../data/mc_%s.root", dataset);
    std::string out_file;
    if(output == nullptr) {
        out_file = Form("fraction_%s_%s.root", particle, dataset);
    } else {
        out_file = output;
    }
    std::cout << "Data file: " << data_file << std::endl;
    std::cout << "MC file: " << mc_file << std::endl;
    std::cout << "Output file: " << out_file << std::endl;

    // --- pT bin edges
    std::vector<double> pt_bins = {0.7, 1.0};
    const int nPtBins = pt_bins.size() - 1;
    std::cout << "Number of pT bins: " << nPtBins << std::endl;

    // --- 打开数据文件
    std::cout << "Opening data file..." << std::endl;
    TFile* fData = TFile::Open(data_file.c_str());
    if(!fData || fData->IsZombie()) {
        std::cerr << "Can't open data file: " << data_file << std::endl;
        return;
    }
    std::cout << "Data file opened successfully" << std::endl;

    std::cout << "Getting data directory..." << std::endl;
    TDirectoryFile* dirData = (TDirectoryFile*)fData->Get("default");
    if(!dirData) {
        std::cerr << "Cannot find 'default' directory in data file" << std::endl;
        return;
    }
    std::cout << "Getting data list..." << std::endl;
    TList* listData = (TList*)dirData->Get("ResultsList_default");
    if(!listData) {
        std::cerr << "Cannot find 'ResultsList_default' list in data file" << std::endl;
        return;
    }

    // 根据particle类型选择正确的直方图
    std::string hist_name = Form("h3_pt_dcaXY_%s", particle);
    std::cout << "Looking for histogram: " << hist_name << std::endl;
    TH3D* h3_data = (TH3D*)listData->FindObject(hist_name.c_str());
    if(!h3_data) {
        std::cerr << "Cannot find histogram: " << hist_name << std::endl;
        return;
    }
    std::cout << "Found data histogram" << std::endl;

    // --- 打开MC文件
    std::cout << "Opening MC file..." << std::endl;
    TFile* fMC = TFile::Open(mc_file.c_str());
    if(!fMC || fMC->IsZombie()) {
        std::cerr << "Can't open MC file: " << mc_file << std::endl;
        return;
    }
    std::cout << "MC file opened successfully" << std::endl;

    std::cout << "Getting MC directory..." << std::endl;
    TDirectoryFile* dirMC = (TDirectoryFile*)fMC->Get("default");
    if(!dirMC) {
        std::cerr << "Cannot find 'default' directory in MC file" << std::endl;
        return;
    }
    std::cout << "Getting MC list..." << std::endl;
    TList* listMC = (TList*)dirMC->Get("ResultsList_default");
    if(!listMC) {
        std::cerr << "Cannot find 'ResultsList_default' list in MC file" << std::endl;
        return;
    }

    // --- 新建输出文件
    std::cout << "Creating output file..." << std::endl;
    TFile* fout = new TFile(out_file.c_str(), "RECREATE");
    if(!fout || fout->IsZombie()) {
        std::cerr << "Cannot create output file: " << out_file << std::endl;
        return;
    }

    // 保存fraction结果到csv文件
    std::string csv_name = Form("./result/fraction_%s_%s.csv", particle, dataset);
    std::cout << "Creating CSV file: " << csv_name << std::endl;
    std::ofstream csv_file(csv_name);
    csv_file << "dataset,particle,cent_low,cent_high,pt_low,pt_high,frac_pri,frac_err,frac_sec,frac_secerr\n";

    // --- 中心度循环
    std::cout << "Starting centrality loop..." << std::endl;
    for(size_t icent = 0; icent < cent_bins.size(); ++icent) {
        double cent_low = cent_bins[icent].first;
        double cent_high = cent_bins[icent].second;
        std::cout << "\nProcessing centrality bin " << icent << " (" << cent_low << "-" << cent_high << "%)" << std::endl;

        // ---- source集中管理 ----
        std::vector<SourceInfo> sources = {
            // hist2d_name                    parname         init   step  min  max     color    particle    dataset    cent_bin
              {"h3_pt_dcaXY_origin_rc_" + std::string(particle), "primary", 0.9, 0.005, 0.0, 1.0, kRed+1, particle, dataset, static_cast<size_t>(icent)}
            , {"h3_pt_dcaXY_material_rc_" + std::string(particle), "secondary", 0.2, 0.005, 0.0, 1.0, kBlue+1, particle, dataset, static_cast<size_t>(icent)}
        };
        const int nTemplates = sources.size();
        std::cout << "Number of templates: " << nTemplates << std::endl;

        // 输出的fraction
        std::vector< std::vector<double> > fractions(nTemplates, std::vector<double>(nPtBins,0));
        std::vector< std::vector<double> > fraction_errors(nTemplates, std::vector<double>(nPtBins,0));

        // --- 逐个pT区间处理
        std::cout << "Starting pT bin loop..." << std::endl;
        for(int ipt=0; ipt<nPtBins; ++ipt){
            double ptlow = pt_bins[ipt];
            double pthigh = pt_bins[ipt+1];
            std::cout << "\nProcessing pT bin " << ipt << " (" << ptlow << "-" << pthigh << " GeV/c)" << std::endl;

            // 数据投影
            std::cout << "Projecting data..." << std::endl;
            TString name_proj_data = Form("proj_data_cent%zu_ptbin%d", icent, ipt);
            int binx1 = h3_data->GetXaxis()->FindBin(cent_low+1e-5);
            int binx2 = h3_data->GetXaxis()->FindBin(cent_high-1e-5);
            int biny1 = h3_data->GetYaxis()->FindBin(ptlow+1e-5);
            int biny2 = h3_data->GetYaxis()->FindBin(pthigh-1e-5);
            std::cout << "Projection bins: x[" << binx1 << "," << binx2 << "], y[" << biny1 << "," << biny2 << "]" << std::endl;

            TH1D* h_data_proj_raw = h3_data->ProjectionZ(name_proj_data, binx1, binx2, biny1, biny2);
            if(!h_data_proj_raw) {
                std::cerr << "Failed to project data histogram" << std::endl;
                continue;
            }
            std::cout << "Data projection successful" << std::endl;
            std::cout << "Raw data histogram bins: " << h_data_proj_raw->GetNbinsX() << std::endl;

            // MC的bin比data更细，直接使用原始data binning
            std::cout << "Using original data binning (MC will be rebinned to match)" << std::endl;
            TH1D* h_data_proj = (TH1D*)h_data_proj_raw->Clone(Form("final_data_cent%zu_ptbin%d", icent, ipt));

            // 输出数据统计信息
            double data_integral = h_data_proj->Integral();
            double data_entries = h_data_proj->GetEntries();
            std::cout << "Data statistics:" << std::endl;
            std::cout << "  Integral: " << data_integral << ", Entries: " << data_entries << std::endl;
            std::cout << "  Mean: " << h_data_proj->GetMean() << ", RMS: " << h_data_proj->GetRMS() << std::endl;

            // MC模板投影
            std::cout << "Projecting MC templates..." << std::endl;
            TObjArray* mc_templates = new TObjArray(nTemplates);
            std::vector<TH1D*> h_mc_proj(nTemplates, nullptr);
            for(int itmp=0; itmp<nTemplates; ++itmp){
                std::cout << "Processing template " << itmp << std::endl;
                TH1D* h_proj = nullptr;
                if(sources[itmp].parname == "primary") {
                    // 原始primary模板
                    std::string hist3d_name = get_mc_template_name("primary", particle);
                    std::cout << "Looking for MC template: " << hist3d_name << std::endl;
                    TH3D* h3tmp = (TH3D*)listMC->FindObject(hist3d_name.c_str());
                    if(!h3tmp) {
                        std::cerr << "Cannot find MC hist " << hist3d_name << std::endl;
                        continue;
                    }
                    h_proj = h3tmp->ProjectionZ(Form("proj_mc_primary_cent%zu_ptbin%d", icent, ipt),
                                              binx1, binx2, biny1, biny2);
                    // Rebin MC template to match data binning
                    std::cout << "Primary MC template bins before rebin: " << h_proj->GetNbinsX() << std::endl;
                    std::cout << "Data bins: " << h_data_proj->GetNbinsX() << std::endl;
                    if(h_proj->GetNbinsX() != h_data_proj->GetNbinsX()) {
                        int rebin_factor = h_proj->GetNbinsX() / h_data_proj->GetNbinsX();
                        std::cout << "Rebinning primary MC template by factor: " << rebin_factor << std::endl;
                        h_proj->Rebin(rebin_factor);
                    }
                } else if(sources[itmp].parname == "secondary") {
                    // 合并secondary
                    std::vector<std::string> sec_types = {"material", "lambda", "other"};
                    TH1D* h_sum = nullptr;
                    bool any_template_found = false;

                    for(const auto& type : sec_types) {
                        std::string hist3d_name = get_mc_template_name(type, particle);
                        std::cout << "Looking for secondary template: " << hist3d_name << std::endl;
                        TH3D* h3tmp = (TH3D*)listMC->FindObject(hist3d_name.c_str());
                        if(!h3tmp) {
                            std::cerr << "Error: Cannot find template " << hist3d_name << std::endl;
                            continue;
                        }
                        any_template_found = true;
                        TH1D* h_tmp_proj = h3tmp->ProjectionZ(Form("proj_mc_%s_cent%zu_ptbin%d",
                                                                  type.c_str(), icent, ipt),
                                                            binx1, binx2, biny1, biny2);
                        // Rebin secondary template to match data binning
                        std::cout << "Secondary template (" << type << ") bins before rebin: " << h_tmp_proj->GetNbinsX() << std::endl;
                        if(h_tmp_proj->GetNbinsX() != h_data_proj->GetNbinsX()) {
                            int rebin_factor = h_tmp_proj->GetNbinsX() / h_data_proj->GetNbinsX();
                            std::cout << "Rebinning secondary template (" << type << ") by factor: " << rebin_factor << std::endl;
                            h_tmp_proj->Rebin(rebin_factor);
                        }
                        if(!h_sum) {
                            h_sum = (TH1D*)h_tmp_proj->Clone(Form("proj_mc_secondary_cent%zu_ptbin%d",
                                                                icent, ipt));
                        } else {
                            h_sum->Add(h_tmp_proj);
                            delete h_tmp_proj;
                        }
                    }

                    if(!any_template_found) {
                        std::cerr << "Error: No secondary templates found!" << std::endl;
                        continue;
                    }

                    if(!h_sum) {
                        std::cerr << "Error: Failed to create secondary sum histogram" << std::endl;
                        continue;
                    }

                    h_proj = h_sum;
                }
                if(!h_proj) {
                    std::cerr << "Failed to create projection for template " << itmp << std::endl;
                    continue;
                }

                // 确保MC模板没有负值
                for(int i=1; i<=h_proj->GetNbinsX(); i++) {
                    if(h_proj->GetBinContent(i) < 0) {
                        h_proj->SetBinContent(i, 0);
                    }
                }

                // 检查直方图内容
                bool hasContent = false;
                for(int i=1; i<=h_proj->GetNbinsX(); i++) {
                    if(h_proj->GetBinContent(i) > 0) {
                        hasContent = true;
                        break;
                    }
                }
                if(!hasContent) {
                    std::cerr << "Warning: Template " << itmp << " has no content!" << std::endl;
                    continue;
                }
                h_mc_proj[itmp] = h_proj;
                mc_templates->Add(h_proj);

                // 输出统计信息
                double integral = h_proj->Integral();
                double entries = h_proj->GetEntries();
                std::cout << "Template " << itmp << " (" << sources[itmp].parname << ") statistics:" << std::endl;
                std::cout << "  Integral: " << integral << ", Entries: " << entries << std::endl;
                std::cout << "  Mean: " << h_proj->GetMean() << ", RMS: " << h_proj->GetRMS() << std::endl;

                // 检查模板质量
                if(integral < 50) {
                    std::cout << "Warning: Low template statistics for " << sources[itmp].parname
                             << " (integral < 50)" << std::endl;
                }
                if(integral < 5) {
                    std::cout << "Error: Extremely low template statistics for " << sources[itmp].parname
                             << " (integral < 5) - this may cause fit problems" << std::endl;
                }

                // 检查模板形状
                double template_mean = h_proj->GetMean();
                double template_rms = h_proj->GetRMS();
                if(template_rms < 0.01) {
                    std::cout << "Warning: Very narrow template distribution (RMS < 0.01)" << std::endl;
                }

                // 检查与数据的重叠度
                double data_mean = h_data_proj->GetMean();
                double data_rms = h_data_proj->GetRMS();
                double mean_diff = TMath::Abs(template_mean - data_mean);
                if(mean_diff > 2.0 * TMath::Max(template_rms, data_rms)) {
                    std::cout << "Warning: Template and data have very different distributions" << std::endl;
                    std::cout << "  Data mean: " << data_mean << " ± " << data_rms << std::endl;
                    std::cout << "  Template mean: " << template_mean << " ± " << template_rms << std::endl;
                }

                std::cout << "Template " << itmp << " processed successfully" << std::endl;
            }

            // 检查是否所有模板都成功处理
            if(mc_templates->GetEntries() != nTemplates) {
                std::cerr << "Error: Not all templates were processed successfully!" << std::endl;
                delete mc_templates;
                continue;
            }

            // --- 拟合前的最终检查
            std::cout << "Pre-fit template comparison:" << std::endl;
            bool fit_feasible = true;
            for(int itmp = 0; itmp < nTemplates; ++itmp) {
                if(h_mc_proj[itmp]) {
                    double ratio = h_mc_proj[itmp]->Integral() / h_data_proj->Integral();
                    std::cout << "  Template " << itmp << " / Data ratio: " << ratio << std::endl;

                    if(ratio < 0.001) {
                        std::cout << "  Warning: Template " << itmp << " has very low statistics relative to data" << std::endl;
                    }
                }
            }

            // 检查模板之间的区分度
            if(nTemplates >= 2 && h_mc_proj[0] && h_mc_proj[1]) {
                // 计算两个模板的KS测试
                double ks_prob = h_mc_proj[0]->KolmogorovTest(h_mc_proj[1]);
                std::cout << "Template separation (KS test p-value): " << ks_prob << std::endl;
                if(ks_prob > 0.5) {
                    std::cout << "Warning: Templates are very similar (KS p-value > 0.5) - fit may be unstable" << std::endl;
                    fit_feasible = false;
                }
            }

            if(!fit_feasible) {
                std::cout << "Warning: Fit conditions are not optimal - results may be unreliable" << std::endl;
            }

            // --- 拟合
            std::cout << "Starting fit..." << std::endl;
            TFractionFitter* fit = new TFractionFitter(h_data_proj, mc_templates);
            if(!fit) {
                std::cerr << "Failed to create TFractionFitter" << std::endl;
                continue;
            }
            ROOT::Fit::Fitter* fitter = fit->GetFitter();
            if(!fitter) {
                std::cerr << "Failed to get fitter" << std::endl;
                continue;
            }

            // 设置Minuit2最小化器的参数
            fitter->Config().MinimizerOptions().SetTolerance(1e-4);  // 提高精度
            fitter->Config().MinimizerOptions().SetMaxFunctionCalls(1000000);  // 适中的迭代次数
            fitter->Config().MinimizerOptions().SetStrategy(2);  // 使用更稳健的策略
            fitter->Config().MinimizerOptions().SetPrintLevel(1);  // 增加输出信息

            // 修改初始值和参数范围 - 根据pT调整初始值
            for(int itmp=0; itmp<nTemplates; ++itmp){
                double init_val, vmin, vmax, step;

                if(sources[itmp].parname == "primary") {
                    // Primary fraction随pT增加而减少
                    if(ptlow < 1.5) {
                        init_val = 0.8;  // 低pT区间primary占主导
                    } else if(ptlow < 3.0) {
                        init_val = 0.6;  // 中等pT
                    } else {
                        init_val = 0.3;  // 高pT区间primary比例较低
                    }
                    vmin = 0.0;
                    vmax = 1.0;
                    step = 0.01;
                } else {  // secondary
                    // Secondary fraction与primary相反
                    if(ptlow < 1.5) {
                        init_val = 0.2;
                    } else if(ptlow < 3.0) {
                        init_val = 0.4;
                    } else {
                        init_val = 0.7;  // 高pT区间secondary比例较高
                    }
                    vmin = 0.0;
                    vmax = 1.0;
                    step = 0.01;
                }

                std::cout << "Setting parameter " << itmp << " (" << sources[itmp].parname
                         << ") initial value: " << init_val << std::endl;

                fitter->Config().ParSettings(itmp).Set(
                    sources[itmp].parname.c_str(),
                    init_val,
                    step,
                    vmin,
                    vmax
                );
            }

            // 检查拟合前的模板质量和归一化
            std::cout << "Final template normalization check:" << std::endl;
            double total_template_integral = 0;
            for(int itmp = 0; itmp < nTemplates; ++itmp) {
                if(h_mc_proj[itmp]) {
                    double template_integral = h_mc_proj[itmp]->Integral();
                    total_template_integral += template_integral;
                    std::cout << "  Template " << itmp << " integral: " << template_integral << std::endl;

                    // 检查是否有空的bin
                    int empty_bins = 0;
                    for(int ibin = 1; ibin <= h_mc_proj[itmp]->GetNbinsX(); ibin++) {
                        if(h_mc_proj[itmp]->GetBinContent(ibin) <= 0) empty_bins++;
                    }
                    double empty_fraction = (double)empty_bins / h_mc_proj[itmp]->GetNbinsX();
                    if(empty_fraction > 0.5) {
                        std::cout << "  Warning: Template " << itmp << " has " << empty_fraction*100
                                 << "% empty bins" << std::endl;
                    }
                }
            }

            std::cout << "Data integral: " << h_data_proj->Integral() << std::endl;
            std::cout << "Total template integrals: " << total_template_integral << std::endl;

            if(total_template_integral < h_data_proj->Integral() * 0.1) {
                std::cout << "Warning: Template statistics much lower than data - fit may be unreliable" << std::endl;
            }



            // 执行拟合
            std::cout << "Executing fit..." << std::endl;
            int fit_status = fit->Fit();

            // 检查拟合收敛性
            std::cout << "Fit status: " << fit_status << std::endl;

            if(fit_status == 0){
                std::cout << "Fit converged successfully!" << std::endl;
                for(int itmp=0; itmp<nTemplates; ++itmp){
                    double val, err;
                    fit->GetResult(itmp, val, err);
                    fractions[itmp][ipt] = val;
                    fraction_errors[itmp][ipt] = err;
                    std::cout << Form("Template %d (%s): %.3f +/- %.3f",
                                    itmp, sources[itmp].parname.c_str(), val, err) << std::endl;
                }
                std::cout << Form("Cent %.0f-%.0f%%, pT %.1f-%.1f: Fit OK.",
                                cent_low, cent_high, ptlow, pthigh) << std::endl;
            }else{
                std::cout << Form("Cent %.0f-%.0f%%, pT %.1f-%.1f: Fit FAILED! Status: %d",
                                cent_low, cent_high, ptlow, pthigh, fit_status) << std::endl;
                for(int itmp=0; itmp<nTemplates; ++itmp){
                    fractions[itmp][ipt] = -1;
                    fraction_errors[itmp][ipt] = 0;
                }
            }

            // --- 手动计算拟合结果
            std::cout << "Calculating fit results..." << std::endl;
            TH1* h_fit = nullptr;
            if(fit_status == 0) {
                h_fit = (TH1D*)h_data_proj->Clone(Form("h_fit_cent%zu_ptbin%d", icent, ipt));
                h_fit->Reset();

                std::vector<TH1D*> h_mc_norm(nTemplates, nullptr);
                for(int itmp=0; itmp<nTemplates; ++itmp) {
                    h_mc_norm[itmp] = (TH1D*)h_mc_proj[itmp]->Clone(Form("h_mc_norm_%d_cent%zu_ptbin%d",
                                                                       itmp, icent, ipt));
                    if(h_mc_norm[itmp]->Integral() > 0) {
                        h_mc_norm[itmp]->Scale(1.0 / h_mc_norm[itmp]->Integral());
                    }
                }

                for(int itmp=0; itmp<nTemplates; ++itmp){
                    double val = fractions[itmp][ipt];
                    if(val > 0) {
                        h_fit->Add(h_mc_norm[itmp], val);
                    }
                }

                if(h_fit->Integral() > 0) {
                    h_fit->Scale(1.0 / h_fit->Integral());
                }

                for(auto h : h_mc_norm) delete h;
            }

            std::cout << "Drawing comparison histograms..." << std::endl;
            DrawFitCompareHistos(h_data_proj, h_mc_proj, sources, h_fit, ipt, ptlow, pthigh);

            // 保存每个投影和模板到输出root
            delete fit;
            delete mc_templates;
            std::cout << "Finished processing pT bin " << ipt << std::endl;
        }

        // 保存fraction结果到root文件
        std::cout << "Saving fraction results to ROOT file..." << std::endl;
        std::vector<double> pt_centers;
        for(int ipt=0; ipt<nPtBins; ++ipt){
            pt_centers.push_back(0.5*(pt_bins[ipt]+pt_bins[ipt+1]));
        }
        for(int itmp=0; itmp<nTemplates; ++itmp){
            TString grname = Form("fraction_%s_%s_%s_cent%zu",
                                sources[itmp].particle.c_str(),
                                sources[itmp].parname.c_str(),
                                sources[itmp].dataset.c_str(),
                                icent);
            TString grtitle = Form(";pT(Gev/c);Fraction");
            TGraphErrors* gr = new TGraphErrors(nPtBins, &pt_centers[0], &fractions[itmp][0],
                                              0, &fraction_errors[itmp][0]);
            gr->SetName(grname);
            gr->SetTitle(grtitle);
            gr->Write();
        }

        // 保存到CSV
        std::cout << "Saving results to CSV..." << std::endl;
        for(int ipt=0; ipt<nPtBins; ++ipt) {
            csv_file << dataset << ","
                     << particle << ","
                     << cent_low << ","
                     << cent_high << ","
                     << pt_bins[ipt] << ","
                     << pt_bins[ipt+1] << ","
                     << fractions[0][ipt] << ","
                     << fraction_errors[0][ipt] << ","
                     << fractions[1][ipt] << ","
                     << fraction_errors[1][ipt] << "\n";
        }
        std::cout << "Finished processing centrality bin " << icent << std::endl;
    }

    csv_file.close();
    fout->Close();
    std::cout << "All done! Output: " << out_file << " and " << csv_name << std::endl;
}
