/***************************************************************************
 *  Differential correlator figure generator - 2023 vs 2025 comparison
 *  -----------------------------------------------------------------
 *  • Compares 2023 and 2025 results for differential correlators (DEta and SPt)
 *  • Only shows Del mode (no SS/OS separation)
 *  • 2023 results shifted left by -1.5, 2025 results shifted right by 1.5
 *  • Based on plot_diff.C structure with additional 2023 data handling
 *
 *  Usage:
 *      root -l -q plot_diff_compare2023.C
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
double legend_SPt_Del_x1 = 0.19;         // SPt Del legend left
double legend_SPt_Del_y1 = 0.44;        // SPt Del legend bottom
double legend_SPt_Del_x2 = 0.65;        // SPt Del legend right
double legend_SPt_Del_y2 = 0.65;        // SPt Del legend top

double legend_DEta_Del_x1 = 0.19;        // DEta Del legend left
double legend_DEta_Del_y1 = 0.44;       // DEta Del legend bottom
double legend_DEta_Del_x2 = 0.65;       // DEta Del legend right
double legend_DEta_Del_y2 = 0.65;       // DEta Del legend top

// Color settings for different years and bins
int color_2023_bin0 = kRed+1;           // 2023 first bin color
int color_2023_bin1 = kBlue+1;          // 2023 second bin color
int color_2023_bin2 = kGreen+3;         // 2023 third bin color

int color_2025_bin0 = kRed+3;           // 2025 first bin color
int color_2025_bin1 = kBlue+3;          // 2025 second bin color
int color_2025_bin2 = kGreen+1;         // 2025 third bin color

// Fill transparency
double alpha_band = 0.2;                // Transparency for systematic error bands
int fill_style = 1000;                  // Fill style for error bands

// Marker settings for different years
int marker_2023_bin0 = kFullCircle;     // 2023 bin 0 marker
int marker_2023_bin1 = kFullSquare;     // 2023 bin 1 marker
int marker_2023_bin2 = kFullDiamond;    // 2023 bin 2 marker

int marker_2025_bin0 = kOpenCircle;     // 2025 bin 0 marker
int marker_2025_bin1 = kOpenSquare;     // 2025 bin 1 marker
int marker_2025_bin2 = kOpenDiamond;    // 2025 bin 2 marker

// Marker size
int markerSize = 2;                     // Size for all markers

// ──────────────────────────────────────────────────────────────────────────
//  Axis-range configuration
// ──────────────────────────────────────────────────────────────────────────
namespace FrameConfig {
struct Range { double xmin, ymin, xmax, ymax; };
constexpr Range SPt_Delta_Del  { 0., -0.002, 60., 0.022 };
constexpr Range SPt_Gamma_Del  { 0., -0.00085, 60., 0.0085 };

constexpr Range DEta_Delta_Del  { 0., -0.001, 60., 0.0145 };
constexpr Range DEta_Gamma_Del  { 0., -0.00032, 60., 0.0032 };

inline const Range& get(const TString& diff_type, const TString& obv)
{
    if (diff_type=="SPt" && obv=="Delta") return SPt_Delta_Del;
    if (diff_type=="SPt" && obv=="Gamma") return SPt_Gamma_Del;
    if (diff_type=="DEta" && obv=="Delta") return DEta_Delta_Del;
    if (diff_type=="DEta" && obv=="Gamma") return DEta_Gamma_Del;

    throw std::runtime_error("FrameConfig::get – unsupported combination");
}
}

// ──────────────────────────────────────────────────────────────────────────
//  Legend-box configuration - helper function
// ──────────────────────────────────────────────────────────────────────────
inline void GetLegendBox(const TString& diff_type, double& x1, double& y1, double& x2, double& y2)
{
    if (diff_type=="SPt") {
        x1 = legend_SPt_Del_x1; y1 = legend_SPt_Del_y1;
        x2 = legend_SPt_Del_x2; y2 = legend_SPt_Del_y2;
    }
    else if (diff_type=="DEta") {
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
void fig_diff_compare2023();
void plot_diff_compare2023();

// ───────────────────────── Data helpers ──────────────────────────────
struct PlotData {
    string diff_type;
    double diff_bin;
    double centrality;
    string pair_type;
    double delta, delta_err, delta_syst_err;
    double gamma, gamma_err, gamma_syst_err;
};

// ShiftAxis工具，仿照fig_diff.C
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

// Get the correct marker style based on year and bin
int GetMarkerStyle(int year, int bin) {
    if(year == 2023) {
        if(bin == 0) return marker_2023_bin0;
        else if(bin == 1) return marker_2023_bin1;
        else return marker_2023_bin2;
    }
    else { // 2025
        if(bin == 0) return marker_2025_bin0;
        else if(bin == 1) return marker_2025_bin1;
        else return marker_2025_bin2;
    }
}

// Get the correct color based on year and bin
int GetColor(int year, int bin) {
    if(year == 2023) {
        if(bin == 0) return color_2023_bin0;
        else if(bin == 1) return color_2023_bin1;
        else return color_2023_bin2;
    }
    else { // 2025
        if(bin == 0) return color_2025_bin0;
        else if(bin == 1) return color_2025_bin1;
        else return color_2025_bin2;
    }
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
void fig_diff_compare2023()
{
    using FrameConfig::get;
    const auto& rng = FrameConfig::get(DiffType, ObvName);

    SetStyle();
    TString fname = DiffType + "_" + ObvName + "_" + PlaneName + "_Del_compare2023.pdf";
    auto* c = new TCanvas(fname, fname, canvas_width, canvas_height); c->SetTickx(1);
    TH1* frame = c->DrawFrame(rng.xmin, rng.ymin, rng.xmax, rng.ymax);
    if(!frame) {cerr<<"DrawFrame failed\n"; return;}

    frame->SetXTitle("Centrality (%)");
    if(ObvName=="Delta") frame->SetYTitle("#Delta#it{#delta} = #it{#delta}_{OS} - #it{#delta}_{SS}");
    if(ObvName=="Gamma") frame->SetYTitle("#Delta#it{#gamma} = #it{#gamma}_{OS} - #it{#gamma}_{SS}");
    frame->GetYaxis()->SetMaxDigits(3);

    // Read 2025 data from CSV
    auto data_2025 = ReadCSVData("../csv_data_point/diff_exp.csv");
    string obs = (ObvName=="Gamma") ? "gamma" : "delta";

    // Determine the number of bins based on the differential type
    int nBins = (DiffType == "SPt") ? 3 : 2;

    // Create histograms for 2025 data
    vector<TH1D*> h2025(nBins, nullptr);
    vector<TH1D*> s2025(nBins, nullptr);
    
    for(int i = 0; i < nBins; i++) {
        double bin_value = i + 0.5; // Convert to 0.5, 1.5, 2.5
        h2025[i] = CreateHist(data_2025, DiffType.Data(), bin_value, "Del", obs, "stat");
        s2025[i] = CreateHist(data_2025, DiffType.Data(), bin_value, "Del", obs, "syst");
    }

    // Load 2023 data from ROOT file (similar to fig_diff.C)
    TFile* inputFile = new TFile("finalResult_diff.root");
    bool has2023Data = false;
    if (!inputFile || inputFile->IsZombie()) {
        cout << "Warning: Cannot open 2023 data file 'finalResult_diff.root'" << endl;
        cout << "Will only show 2025 data" << endl;
    } else {
        has2023Data = true;
        cout << "Successfully opened 2023 data file" << endl;
    }

    vector<TString> vecKineBins;
    if (DiffType.EqualTo("SPt")) {
        vecKineBins = {"SPtkineBin0", "SPtkineBin1", "SPtkineBin2"};
    } else if (DiffType.EqualTo("DEta")) {
        vecKineBins = {"DEtakineBin0", "DEtakineBin1"};
    }

    vector<TH1D*> h2023(nBins, nullptr);
    vector<TH1D*> s2023(nBins, nullptr);

    if (has2023Data) {
        for (int i = 0; i < nBins; i++) {
            TString histName = "h_18qr_"+ ObvName + "_" + PlaneName + "_Del_" + vecKineBins[i] + "_default";
            h2023[i] = (TH1D*)inputFile->Get(histName);
            if (!h2023[i]) {
                cout << "Warning: Cannot find 2023 histogram: " << histName << endl;
                has2023Data = false;
                break;
            }
            
            TString systName = "h"+ObvName+"DelSystError_" + PlaneName + "_" + vecKineBins[i];
            s2023[i] = (TH1D*)inputFile->Get(systName);
            if (!s2023[i]) {
                cout << "Warning: Cannot find 2023 systematic error: " << systName << endl;
                has2023Data = false;
                break;
            }
            cout << "Found 2023 data for bin " << i << endl;
        }
    }

    // Debug: Check 2025 data
    cout << "2025 data check:" << endl;
    for(int i = 0; i < nBins; i++) {
        if(h2025[i]) {
            cout << "  Bin " << i << ": " << h2025[i]->GetEntries() << " entries" << endl;
        } else {
            cout << "  Bin " << i << ": No data found" << endl;
        }
    }

    // Apply styles and colors
    for(int i = 0; i < nBins; i++) {
        // 2023 data styling
        if(h2023[i] && has2023Data) {
            SetStyle(h2023[i], GetColor(2023, i), GetMarkerStyle(2023, i), markerSize);
        }
        if(s2023[i] && h2023[i] && has2023Data) {
            SetStyle(s2023[i], h2023[i]);
        }

        // 2025 data styling
        if(h2025[i]) {
            SetStyle(h2025[i], GetColor(2025, i), GetMarkerStyle(2025, i), markerSize);
        }
        if(s2025[i] && h2025[i]) {
            SetStyle(s2025[i], h2025[i]);
        }
    }

    // ====== 应用ShiftAxis ======
    // 2023年向左移动-1.5，2025年向右移动1.5
    double dx_2023 = -1.5;
    double dx_2025 = 1.5;
    
    for(int i = 0; i < nBins; i++) {
        if(h2023[i] && has2023Data) {
            ShiftAxis(h2023[i], dx_2023);
            ShiftAxis(s2023[i], dx_2023);
        }
        if(h2025[i]) {
            ShiftAxis(h2025[i], dx_2025);
            ShiftAxis(s2025[i], dx_2025);
        }
    }

    /* draw */
    c->cd();
    
    // Draw systematic error bands first (2023)
    if(has2023Data) {
        for(size_t i = 0; i < s2023.size(); i++) {
            if(s2023[i]) s2023[i]->Draw("E2,same");
        }
    }
    
    // Draw systematic error bands (2025)
    for(size_t i = 0; i < s2025.size(); i++) {
        if(s2025[i]) s2025[i]->Draw("E2,same");
    }

    // Draw data points (2023)
    if(has2023Data) {
        for(size_t i = 0; i < h2023.size(); i++) {
            if(h2023[i]) h2023[i]->Draw("E X0,same");
        }
    }
    
    // Draw data points (2025)
    for(size_t i = 0; i < h2025.size(); i++) {
        if(h2025[i]) h2025[i]->Draw("E X0,same");
    }

    /* ALICE Preliminary tag */
    auto* pt = new TPaveText(alice_x1, alice_y1, alice_x2, alice_y2, "brNDC");
    pt->SetFillColor(0); pt->SetBorderSize(0); pt->SetShadowColor(0);
    pt->SetTextSize(gStyle->GetTextSize() * alice_text_size); pt->SetTextAlign(12);
    if(has2023Data) {
        pt->AddText("ALICE CVE 2023 vs 2025");
    } else {
        pt->AddText("ALICE CVE 2025");
    }
    pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pt->Draw();

    /* Legend positioned using global settings */
    double leg_x1, leg_y1, leg_x2, leg_y2;
    GetLegendBox(DiffType, leg_x1, leg_y1, leg_x2, leg_y2);
    auto* leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    leg->SetTextSize(legend_text_size * gStyle->GetTextSize());

    // Add legend entries for 2023 data
    if(has2023Data) {
        for(size_t i = 0; i < h2023.size(); i++) {
            if(h2023[i]) {
                string label = "2023 " + GetBinLabel(DiffType.Data(), i);
                leg->AddEntry(h2023[i], label.c_str(), "PE");
            }
        }
    }
    
    // Add legend entries for 2025 data
    for(size_t i = 0; i < h2025.size(); i++) {
        if(h2025[i]) {
            string label = "2025 " + GetBinLabel(DiffType.Data(), i);
            leg->AddEntry(h2025[i], label.c_str(), "PE");
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
void plot_diff_compare2023()
{
    cout<<"Creating SPt Delta Del comparison ...\n"; DiffType="SPt"; ObvName="Delta"; fig_diff_compare2023();
    cout<<"Creating SPt Gamma Del comparison ...\n";  DiffType="SPt"; ObvName="Gamma"; fig_diff_compare2023();

    cout<<"Creating DEta Delta Del comparison ...\n"; DiffType="DEta"; ObvName="Delta"; fig_diff_compare2023();
    cout<<"Creating DEta Gamma Del comparison ...\n";  DiffType="DEta"; ObvName="Gamma"; fig_diff_compare2023();

    cout<<"All comparison plots completed!\n";
}
