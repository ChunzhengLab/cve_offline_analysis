/***************************************************************************
 *  Differential correlator figure generator
 *  -----------------------------------------------------------------
 *  • Reads data from "../csv_data_point/diff.csv"
 *  • Creates plots for differential correlators (DEta and SPt)
 *  • All styling configs are placed at the top for easy access
 *  • No shift axis functionality (markers show at actual positions)
 *
 *  Usage:
 *      root -l -q plot_diff.C
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>

#include "TROOT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveText.h"

// ──────────────────────────────────────────────────────────────────────────
//  Style configuration - direct color and marker settings for easy modification
// ──────────────────────────────────────────────────────────────────────────
// Canvas settings
int canvas_width = 1000;                  // Canvas width in pixels
int canvas_height = 750;                 // Canvas height in pixels

// Text and legend settings
double legend_text_size = 0.6;           // Legend text size (relative to gStyle)
double alice_text_size = 0.6;            // ALICE text size (relative to gStyle)

// TPaveText position for ALICE tag
double alice_x1 = 0.2;                   // Left position
double alice_y1 = 0.65;                  // Bottom position
double alice_x2 = 0.4;                   // Right position
double alice_y2 = 0.85;                  // Top position

// Legend positions for different modes
double legend_SPt_SSOS_x1 = 0.19;        // SPt SSOS legend left
double legend_SPt_SSOS_y1 = 0.19;       // SPt SSOS legend bottom
double legend_SPt_SSOS_x2 = 0.65;       // SPt SSOS legend right
double legend_SPt_SSOS_y2 = 0.40;       // SPt SSOS legend top

double legend_SPt_Del_x1 = 0.19;         // SPt Del legend left
double legend_SPt_Del_y1 = 0.44;        // SPt Del legend bottom
double legend_SPt_Del_x2 = 0.65;        // SPt Del legend right
double legend_SPt_Del_y2 = 0.65;        // SPt Del legend top

double legend_DEta_SSOS_x1 = 0.19;       // DEta SSOS legend left
double legend_DEta_SSOS_y1 = 0.19;      // DEta SSOS legend bottom
double legend_DEta_SSOS_x2 = 0.65;      // DEta SSOS legend right
double legend_DEta_SSOS_y2 = 0.40;      // DEta SSOS legend top

double legend_DEta_Del_x1 = 0.19;        // DEta Del legend left
double legend_DEta_Del_y1 = 0.44;       // DEta Del legend bottom
double legend_DEta_Del_x2 = 0.65;       // DEta Del legend right
double legend_DEta_Del_y2 = 0.65;       // DEta Del legend top

// Color settings for different kinematic bins
int color_bin0 = kRed+1;                // First bin color
int color_bin1 = kBlue+1;               // Second bin color
int color_bin2 = kGreen+3;              // Third bin color

// Fill transparency
double alpha_band = 0.2;                // Transparency for systematic error bands
int fill_style = 1000;                  // Fill style for error bands

// Marker settings for SS, OS and Del
int marker_SS_bin0 = kOpenCircle;       // SS marker for bin 0
int marker_OS_bin0 = kFullCircle;       // OS marker for bin 0
int marker_Del_bin0 = kFullCircle;      // Del marker for bin 0

int marker_SS_bin1 = kOpenSquare;       // SS marker for bin 1
int marker_OS_bin1 = kFullSquare;       // OS marker for bin 1
int marker_Del_bin1 = kFullSquare;      // Del marker for bin 1

int marker_SS_bin2 = kOpenDiamond;      // SS marker for bin 2
int marker_OS_bin2 = kFullDiamond;      // OS marker for bin 2
int marker_Del_bin2 = kFullDiamond;     // Del marker for bin 2

// Marker size
int markerSize = 2;                     // Size for all markers

// ──────────────────────────────────────────────────────────────────────────
//  Axis-range configuration
// ──────────────────────────────────────────────────────────────────────────
namespace FrameConfig {
struct Range { double xmin, ymin, xmax, ymax; };
constexpr Range SPt_Delta_SSOS { 0., -0.0195, 60., 0.0195 };
constexpr Range SPt_Delta_Del  { 0., -0.002, 60., 0.022 };
constexpr Range SPt_Gamma_SSOS { 0., -0.007, 60., 0.003 };
constexpr Range SPt_Gamma_Del  { 0., -0.00085, 60., 0.0085 };

constexpr Range DEta_Delta_SSOS { 0., -0.012, 60., 0.012 };
constexpr Range DEta_Delta_Del  { 0., -0.001, 60., 0.0145 };
constexpr Range DEta_Gamma_SSOS { 0., -0.004, 60., 0.001 };
constexpr Range DEta_Gamma_Del  { 0., -0.00032, 60., 0.0032 };

inline const Range& get(const TString& diff_type, const TString& obv, const TString& mode)
{
    if (diff_type=="SPt" && obv=="Delta" && mode=="SSOS") return SPt_Delta_SSOS;
    if (diff_type=="SPt" && obv=="Delta" && mode=="Del")  return SPt_Delta_Del;
    if (diff_type=="SPt" && obv=="Gamma" && mode=="SSOS") return SPt_Gamma_SSOS;
    if (diff_type=="SPt" && obv=="Gamma" && mode=="Del")  return SPt_Gamma_Del;

    if (diff_type=="DEta" && obv=="Delta" && mode=="SSOS") return DEta_Delta_SSOS;
    if (diff_type=="DEta" && obv=="Delta" && mode=="Del")  return DEta_Delta_Del;
    if (diff_type=="DEta" && obv=="Gamma" && mode=="SSOS") return DEta_Gamma_SSOS;
    if (diff_type=="DEta" && obv=="Gamma" && mode=="Del")  return DEta_Gamma_Del;

    throw std::runtime_error("FrameConfig::get – unsupported combination");
}
}

// ──────────────────────────────────────────────────────────────────────────
//  Legend-box configuration - helper function
// ──────────────────────────────────────────────────────────────────────────
inline void GetLegendBox(const TString& diff_type, const TString& mode, double& x1, double& y1, double& x2, double& y2)
{
    if (diff_type=="SPt" && mode=="SSOS") {
        x1 = legend_SPt_SSOS_x1; y1 = legend_SPt_SSOS_y1;
        x2 = legend_SPt_SSOS_x2; y2 = legend_SPt_SSOS_y2;
    }
    else if (diff_type=="SPt" && mode=="Del") {
        x1 = legend_SPt_Del_x1; y1 = legend_SPt_Del_y1;
        x2 = legend_SPt_Del_x2; y2 = legend_SPt_Del_y2;
    }
    else if (diff_type=="DEta" && mode=="SSOS") {
        x1 = legend_DEta_SSOS_x1; y1 = legend_DEta_SSOS_y1;
        x2 = legend_DEta_SSOS_x2; y2 = legend_DEta_SSOS_y2;
    }
    else if (diff_type=="DEta" && mode=="Del") {
        x1 = legend_DEta_Del_x1; y1 = legend_DEta_Del_y1;
        x2 = legend_DEta_Del_x2; y2 = legend_DEta_Del_y2;
    }
    else throw std::runtime_error("GetLegendBox – unsupported combination");
}

// ──────────────────────────────────────────────────────────────────────────
//  Globals (overridden per call)
// ──────────────────────────────────────────────────────────────────────────
TString DiffType  = "DEta";    // ["DEta","SPt"]
TString ObvName   = "Delta";   // ["Delta","Gamma"]
TString PlaneName = "TPC";     // kept for possible future use
TString SSOSDel   = "SSOS";    // ["SSOS","Del"]

// Constants for bin labels
const std::vector<std::string> SPtBinLabels = {
    "1 < #it{p}_{T}^{sum} < 3 GeV/#it{c}",
    "3 < #it{p}_{T}^{sum} < 5 GeV/#it{c}",
    "5 < #it{p}_{T}^{sum} < 8 GeV/#it{c}"
};

const std::vector<std::string> DEtaBinLabels = {
    "|#Delta#eta| < 0.6",
    "|#Delta#eta| > 0.6"
};

using namespace std;

// ───────────────────────────── Forward decls ─────────────────────────────
void SetStyle(Bool_t graypalette=kFALSE);
void fig_diff();
void plot_diff();

// ───────────────────────── Data helpers ──────────────────────────────
struct PlotData {
    string diff_type;
    double diff_bin;
    double centrality;
    string pair_type;
    double delta, delta_err, delta_syst_err;
    double gamma, gamma_err, gamma_syst_err;
};

// ShiftAxis工具，仿照fig_intg.C
void ShiftAxis(TH1D* h, double dx) {
    if(!h) return;
    TAxis* a = h->GetXaxis();
    a->Set(a->GetNbins(), a->GetXmin() + dx, a->GetXmax() + dx);
}

vector<PlotData> ReadCSVData(const string& fn)
{
    vector<PlotData> v; ifstream f(fn); string line; getline(f,line);
    while (getline(f,line)) {
        stringstream ss(line); string it; PlotData p;
        getline(ss,p.diff_type,',');
        getline(ss,it,','); p.diff_bin=stod(it);
        getline(ss,it,','); p.centrality=stod(it);
        getline(ss,p.pair_type,',');
        getline(ss,it,','); p.delta=stod(it);
        getline(ss,it,','); p.delta_err=stod(it);
        getline(ss,it,','); p.delta_syst_err=stod(it);
        getline(ss,it,','); p.gamma=stod(it);
        getline(ss,it,','); p.gamma_err=stod(it);
        getline(ss,it,','); p.gamma_syst_err=stod(it);
        v.push_back(p);
    }
    return v;
}

TH1D* CreateHist(const vector<PlotData>& data, const string& diff_type, double diff_bin,
                 const string& pair_type, const string& obs, const string& err)
{
    vector<double> cx, val, errv;
    for(const auto& p:data) {
        if(p.diff_type==diff_type && p.diff_bin==diff_bin && p.pair_type==pair_type) {
            cx.push_back(p.centrality);
            if(obs=="delta"){
                val.push_back(p.delta);
                errv.push_back(err=="stat"?p.delta_err:p.delta_syst_err);
            }
            else {
                val.push_back(p.gamma);
                errv.push_back(err=="stat"?p.gamma_err:p.gamma_syst_err);
            }
        }
    }

    if(cx.empty()) return nullptr;

    string n = diff_type + "_" + to_string(diff_bin) + "_" + pair_type + "_" + obs + "_" + err;
    auto* h = new TH1D(n.c_str(), n.c_str(), (int)cx.size(), 0, 60);

    for(size_t i=0; i<cx.size(); ++i) {
        int b = h->FindBin(cx[i]);
        h->SetBinContent(b, val[i]);
        h->SetBinError(b, errv[i]);
    }

    return h;
}

// No shift is applied to histograms on x-axis

// Style setting helpers for histograms
template<class TH> void SetStyle(TH& h, unsigned col, int mk, double ms, int ls=1, double lw=1) {
    h->SetLineColor(col); h->SetMarkerColor(col);
    h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->SetLineStyle(ls);   h->SetLineWidth(lw);
}

template<class TH> void SetStyle(TH& hS, TH& h) {
    hS->SetFillColorAlpha(h->GetLineColor(), alpha_band);
    hS->SetFillStyle(fill_style);
    hS->SetMarkerStyle(20);
    hS->SetMarkerSize(0);
}

// Get the correct marker style based on bin and pair type
int GetMarkerStyle(const string& pair_type, int bin) {
    if(pair_type == "SS") {
        if(bin == 0) return marker_SS_bin0;
        else if(bin == 1) return marker_SS_bin1;
        else return marker_SS_bin2;
    }
    else if(pair_type == "OS") {
        if(bin == 0) return marker_OS_bin0;
        else if(bin == 1) return marker_OS_bin1;
        else return marker_OS_bin2;
    }
    else { // Del
        if(bin == 0) return marker_Del_bin0;
        else if(bin == 1) return marker_Del_bin1;
        else return marker_Del_bin2;
    }
}

// Get the correct color based on bin
int GetColor(int bin) {
    if(bin == 0) return color_bin0;
    else if(bin == 1) return color_bin1;
    else return color_bin2;
}

// Get the correct bin label
string GetBinLabel(const string& diff_type, int bin) {
    if(diff_type == "SPt") {
        if(bin < SPtBinLabels.size()) return SPtBinLabels[bin];
    }
    else if(diff_type == "DEta") {
        if(bin < DEtaBinLabels.size()) return DEtaBinLabels[bin];
    }
    return "Unknown bin";
}

// ───────────────────────── Main plotting routine ─────────────────────────
void fig_diff()
{
    using FrameConfig::get;
    const auto& rng = FrameConfig::get(DiffType, ObvName, SSOSDel);

    SetStyle();
    TString fname = DiffType + "_" + ObvName + "_" + PlaneName + "_" + SSOSDel + ".pdf";
    auto* c = new TCanvas(fname, fname, canvas_width, canvas_height); c->SetTickx(1);
    TH1* frame = c->DrawFrame(rng.xmin, rng.ymin, rng.xmax, rng.ymax);
    if(!frame) {cerr<<"DrawFrame failed\n"; return;}

    frame->SetXTitle("Centrality (%)");
    if(ObvName=="Delta"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#delta} = #LTcos(#it{#phi}_{#alpha} - #it{#phi}_{#beta})#GT");
    if(ObvName=="Delta"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#delta} = #it{#delta}_{OS} - #it{#delta}_{SS}");
    if(ObvName=="Gamma"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#gamma} = #LTcos(#it{#phi}_{#alpha} + #it{#phi}_{#beta} - 2#Psi_{2})#GT");
    if(ObvName=="Gamma"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#gamma} = #it{#gamma}_{OS} - #it{#gamma}_{SS}");
    frame->GetYaxis()->SetMaxDigits(3);

    auto data = ReadCSVData("../csv_data_point/diff_exp.csv");
    string obs = (ObvName=="Gamma") ? "gamma" : "delta";

    // Determine the number of bins based on the differential type
    int nBins = (DiffType == "SPt") ? 3 : 2;

    // Create a map to store histograms for each bin
    map<string, vector<TH1D*>> histograms;

    // For SSOS mode, we need SS and OS histograms
    if(SSOSDel=="SSOS") {
        histograms["hSS"] = vector<TH1D*>(nBins, nullptr);
        histograms["hOS"] = vector<TH1D*>(nBins, nullptr);
        histograms["sSS"] = vector<TH1D*>(nBins, nullptr);
        histograms["sOS"] = vector<TH1D*>(nBins, nullptr);

        for(int i = 0; i < nBins; i++) {
            double bin_value = i + 0.5; // Convert to 0.5, 1.5, 2.5
            histograms["hSS"][i] = CreateHist(data, DiffType.Data(), bin_value, "SS", obs, "stat");
            histograms["hOS"][i] = CreateHist(data, DiffType.Data(), bin_value, "OS", obs, "stat");
            histograms["sSS"][i] = CreateHist(data, DiffType.Data(), bin_value, "SS", obs, "syst");
            histograms["sOS"][i] = CreateHist(data, DiffType.Data(), bin_value, "OS", obs, "syst");
        }
    } else {
        histograms["hDel"] = vector<TH1D*>(nBins, nullptr);
        histograms["sDel"] = vector<TH1D*>(nBins, nullptr);

        for(int i = 0; i < nBins; i++) {
            double bin_value = i + 0.5; // Convert to 0.5, 1.5, 2.5
            histograms["hDel"][i] = CreateHist(data, DiffType.Data(), bin_value, "Del", obs, "stat");
            histograms["sDel"][i] = CreateHist(data, DiffType.Data(), bin_value, "Del", obs, "syst");
        }
    }

    // Apply styles and colors
    for(int i = 0; i < nBins; i++) {
        int colorIndex = i;
        int markerIndex = i;

        if(SSOSDel=="SSOS") {
            if(histograms["hSS"][i]) {
                SetStyle(histograms["hSS"][i], GetColor(i), GetMarkerStyle("SS", i), markerSize);
            }
            if(histograms["hOS"][i]) {
                SetStyle(histograms["hOS"][i], GetColor(i), GetMarkerStyle("OS", i), markerSize);
            }
            if(histograms["sSS"][i] && histograms["hSS"][i]) {
                SetStyle(histograms["sSS"][i], histograms["hSS"][i]);
            }
            if(histograms["sOS"][i] && histograms["hOS"][i]) {
                SetStyle(histograms["sOS"][i], histograms["hOS"][i]);
            }
        } else {
            if(histograms["hDel"][i]) {
                SetStyle(histograms["hDel"][i], GetColor(i), GetMarkerStyle("Del", i), markerSize);
            }
            if(histograms["sDel"][i] && histograms["hDel"][i]) {
                SetStyle(histograms["sDel"][i], histograms["hDel"][i]);
            }
        }
    }

    // ====== 这里做ShiftAxis ======
    // 对于DEta，两个结果分别左右移动1.5
    // 对于SPt，三个结果中间不移动，左右移动-2，2
    if(DiffType=="DEta") {
        double dx_bin0 = -1.5; // 左移
        double dx_bin1 = 1.5;  // 右移

        if(SSOSDel=="SSOS") {
            ShiftAxis(histograms["hSS"][0], dx_bin0); ShiftAxis(histograms["sSS"][0], dx_bin0);
            ShiftAxis(histograms["hOS"][0], dx_bin0); ShiftAxis(histograms["sOS"][0], dx_bin0);
            ShiftAxis(histograms["hSS"][1], dx_bin1); ShiftAxis(histograms["sSS"][1], dx_bin1);
            ShiftAxis(histograms["hOS"][1], dx_bin1); ShiftAxis(histograms["sOS"][1], dx_bin1);
        } else {
            ShiftAxis(histograms["hDel"][0], dx_bin0); ShiftAxis(histograms["sDel"][0], dx_bin0);
            ShiftAxis(histograms["hDel"][1], dx_bin1); ShiftAxis(histograms["sDel"][1], dx_bin1);
        }
    } else if(DiffType=="SPt") {
        double dx_bin0 = -2.0; // 左移
        double dx_bin1 = 0.0;  // 中间不移动
        double dx_bin2 = 2.0;  // 右移

        if(SSOSDel=="SSOS") {
            ShiftAxis(histograms["hSS"][0], dx_bin0); ShiftAxis(histograms["sSS"][0], dx_bin0);
            ShiftAxis(histograms["hOS"][0], dx_bin0); ShiftAxis(histograms["sOS"][0], dx_bin0);
            // bin1不移动
            ShiftAxis(histograms["hSS"][2], dx_bin2); ShiftAxis(histograms["sSS"][2], dx_bin2);
            ShiftAxis(histograms["hOS"][2], dx_bin2); ShiftAxis(histograms["sOS"][2], dx_bin2);
        } else {
            ShiftAxis(histograms["hDel"][0], dx_bin0); ShiftAxis(histograms["sDel"][0], dx_bin0);
            // bin1不移动
            ShiftAxis(histograms["hDel"][2], dx_bin2); ShiftAxis(histograms["sDel"][2], dx_bin2);
        }
    }

    /* draw */
    c->cd();
    if(SSOSDel=="SSOS") {
        // Draw systematic error bands first
        for(size_t i = 0; i < histograms["sSS"].size(); i++) {
            if(histograms["sSS"][i]) histograms["sSS"][i]->Draw("E2,same");
        }
        for(size_t i = 0; i < histograms["sOS"].size(); i++) {
            if(histograms["sOS"][i]) histograms["sOS"][i]->Draw("E2,same");
        }

        // Then draw data points
        for(size_t i = 0; i < histograms["hSS"].size(); i++) {
            if(histograms["hSS"][i]) histograms["hSS"][i]->Draw("E X0,same");
        }
        for(size_t i = 0; i < histograms["hOS"].size(); i++) {
            if(histograms["hOS"][i]) histograms["hOS"][i]->Draw("E X0,same");
        }
    } else {
        // Draw systematic error bands first
        for(size_t i = 0; i < histograms["sDel"].size(); i++) {
            if(histograms["sDel"][i]) histograms["sDel"][i]->Draw("E2,same");
        }

        // Then draw data points
        for(size_t i = 0; i < histograms["hDel"].size(); i++) {
            if(histograms["hDel"][i]) histograms["hDel"][i]->Draw("E X0,same");
        }
    }

    /* ALICE Preliminary tag */
    auto* pt = new TPaveText(alice_x1, alice_y1, alice_x2, alice_y2, "brNDC");
    pt->SetFillColor(0); pt->SetBorderSize(0); pt->SetShadowColor(0);
    pt->SetTextSize(gStyle->GetTextSize() * alice_text_size); pt->SetTextAlign(12);
    pt->AddText("ALICE CVE 2025");
    pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pt->Draw();

    /* Legend positioned using global settings */
    double leg_x1, leg_y1, leg_x2, leg_y2;
    GetLegendBox(DiffType, SSOSDel, leg_x1, leg_y1, leg_x2, leg_y2);
    auto* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    leg->SetTextSize(legend_text_size * gStyle->GetTextSize());

    if(SSOSDel=="SSOS") {
        for(size_t i = 0; i < histograms["hOS"].size(); i++) {
            if(histograms["hOS"][i]) {
                string label = "OS " + GetBinLabel(DiffType.Data(), i);
                leg->AddEntry(histograms["hOS"][i], label.c_str(), "PE");
            }
        }
        for(size_t i = 0; i < histograms["hSS"].size(); i++) {
            if(histograms["hSS"][i]) {
                string label = "SS " + GetBinLabel(DiffType.Data(), i);
                leg->AddEntry(histograms["hSS"][i], label.c_str(), "PE");
            }
        }
    } else {
        for(size_t i = 0; i < histograms["hDel"].size(); i++) {
            if(histograms["hDel"][i]) {
                leg->AddEntry(histograms["hDel"][i],
                              (GetBinLabel(DiffType.Data(), i)).c_str(), "PE");
            }
        }
    }
    leg->Draw();

    c->SaveAs(fname); delete c;
    cout<<" → saved "<<fname<<'\n';
}

// ───────────────────────── ROOT style util ───────────────────────────────
void SetStyle(Bool_t graypalette)
{
    gStyle->Reset("Plain"); gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
    graypalette ? gStyle->SetPalette(8,0) : gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10); gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadColor(10); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15); gStyle->SetPadLeftMargin(0.15);
    gStyle->SetLabelSize(0.048,"xyz"); gStyle->SetTextFont(42);
    gStyle->SetTitleSize(0.05,"xyz"); gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite); gStyle->SetErrorX(0.075);
    gStyle->SetNdivisions(505,"y");
}

// ───────────────────────── Driver to build all canvases ─────────────────
void plot_diff()
{
    cout<<"Creating SPt Delta SSOS ...\n"; DiffType="SPt"; ObvName="Delta"; SSOSDel="SSOS"; fig_diff();
    cout<<"Creating SPt Delta Del ...\n";  DiffType="SPt"; ObvName="Delta"; SSOSDel="Del";  fig_diff();
    cout<<"Creating SPt Gamma SSOS ...\n"; DiffType="SPt"; ObvName="Gamma"; SSOSDel="SSOS"; fig_diff();
    cout<<"Creating SPt Gamma Del ...\n";  DiffType="SPt"; ObvName="Gamma"; SSOSDel="Del";  fig_diff();

    cout<<"Creating DEta Delta SSOS ...\n"; DiffType="DEta"; ObvName="Delta"; SSOSDel="SSOS"; fig_diff();
    cout<<"Creating DEta Delta Del ...\n";  DiffType="DEta"; ObvName="Delta"; SSOSDel="Del";  fig_diff();
    cout<<"Creating DEta Gamma SSOS ...\n"; DiffType="DEta"; ObvName="Gamma"; SSOSDel="SSOS"; fig_diff();
    cout<<"Creating DEta Gamma Del ...\n";  DiffType="DEta"; ObvName="Gamma"; SSOSDel="Del";  fig_diff();

    cout<<"All plots completed!\n";
}
