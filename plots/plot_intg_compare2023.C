/***************************************************************************
 *  Integrated correlator figure generator - 2023 vs 2025 comparison
 *  -----------------------------------------------------------------
 *  • Compares 2023 and 2025 results for Lambda-Proton and Lambda-Hadron
 *  • 2023 results shifted left by -1.5, 2025 results shifted right by 1.5
 *  • Based on plot_intg.C structure with additional 2023 data handling
 *
 *  Usage:
 *      root -l -q plot_intg_compare2023.C
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

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

// Legend positions (all modes)
double legend_Delta_SSOS_x1 = 0.19;      // Delta SSOS legend left
double legend_Delta_SSOS_y1 = 0.19;      // Delta SSOS legend bottom
double legend_Delta_SSOS_x2 = 0.65;      // Delta SSOS legend right
double legend_Delta_SSOS_y2 = 0.40;      // Delta SSOS legend top

double legend_Delta_Del_x1 = 0.19;       // Delta Del legend left
double legend_Delta_Del_y1 = 0.44;       // Delta Del legend bottom
double legend_Delta_Del_x2 = 0.65;       // Delta Del legend right
double legend_Delta_Del_y2 = 0.65;       // Delta Del legend top

double legend_Gamma_SSOS_x1 = 0.19;      // Gamma SSOS legend left
double legend_Gamma_SSOS_y1 = 0.19;      // Gamma SSOS legend bottom
double legend_Gamma_SSOS_x2 = 0.65;      // Gamma SSOS legend right
double legend_Gamma_SSOS_y2 = 0.40;      // Gamma SSOS legend top

double legend_Gamma_Del_x1 = 0.19;       // Gamma Del legend left
double legend_Gamma_Del_y1 = 0.44;       // Gamma Del legend bottom
double legend_Gamma_Del_x2 = 0.65;       // Gamma Del legend right
double legend_Gamma_Del_y2 = 0.65;       // Gamma Del legend top

// Color settings for different years and particle pairs
int color_2023_LambdaProton = kRed+1;    // 2023 Lambda-Proton color
int color_2025_LambdaProton = kRed+3;    // 2025 Lambda-Proton color
int color_2023_LambdaHadron = kBlue+1;   // 2023 Lambda-Hadron color
int color_2025_LambdaHadron = kBlue+3;   // 2025 Lambda-Hadron color

// Fill transparency
double alpha_band = 0.2;                 // Transparency for error bands
int fill_style = 1000;                   // Fill style for error bands

// Marker settings - 2023 vs 2025
int marker_2023_LP_SS = kOpenSquare;     // 2023 Lambda-Proton Same Sign marker
int marker_2023_LP_OS = kFullSquare;     // 2023 Lambda-Proton Opposite Sign marker
int marker_2023_LP_Del = kFullSquare;    // 2023 Lambda-Proton Delta marker

int marker_2025_LP_SS = kOpenCircle;     // 2025 Lambda-Proton Same Sign marker
int marker_2025_LP_OS = kFullCircle;     // 2025 Lambda-Proton Opposite Sign marker
int marker_2025_LP_Del = kFullCircle;    // 2025 Lambda-Proton Delta marker

int marker_2023_LH_SS = kOpenTriangleUp; // 2023 Lambda-Hadron SS marker
int marker_2023_LH_OS = kFullTriangleDown; // 2023 Lambda-Hadron OS marker
int marker_2023_LH_Del = kFullTriangleUp; // 2023 Lambda-Hadron Delta marker

int marker_2025_LH_SS = kOpenTriangleDown; // 2025 Lambda-Hadron SS marker
int marker_2025_LH_OS = kFullTriangleUp;   // 2025 Lambda-Hadron OS marker
int marker_2025_LH_Del = kFullTriangleDown; // 2025 Lambda-Hadron Delta marker

// ──────────────────────────────────────────────────────────────────────────
//  Axis-range configuration
// ──────────────────────────────────────────────────────────────────────────
namespace FrameConfig {
struct Range { double xmin, ymin, xmax, ymax; };
constexpr Range Delta_SSOS { 0. , -0.01 , 60. ,  0.01};
constexpr Range Delta_Del  { 0. , -0.002 , 60. ,  0.015};
constexpr Range Gamma_SSOS { 0. , -0.0038,  60. ,  0.002};
constexpr Range Gamma_Del  { 0. , -0.0003, 60. ,  0.003};

inline const Range& get(const TString& obv, const TString& mode)
{
    if (obv=="Delta" && mode=="SSOS") return Delta_SSOS;
    if (obv=="Delta" && mode=="Del")  return Delta_Del;
    if (obv=="Gamma" && mode=="SSOS") return Gamma_SSOS;
    if (obv=="Gamma" && mode=="Del")  return Gamma_Del;
    throw std::runtime_error("FrameConfig::get – unsupported combination");
}
}

// ──────────────────────────────────────────────────────────────────────────
//  Legend-box configuration - helper function
// ──────────────────────────────────────────────────────────────────────────
inline void GetLegendBox(const TString& obv, const TString& mode, double& x1, double& y1, double& x2, double& y2)
{
    if (obv=="Delta" && mode=="SSOS") { x1 = legend_Delta_SSOS_x1; y1 = legend_Delta_SSOS_y1; x2 = legend_Delta_SSOS_x2; y2 = legend_Delta_SSOS_y2; }
    else if (obv=="Delta" && mode=="Del") { x1 = legend_Delta_Del_x1; y1 = legend_Delta_Del_y1; x2 = legend_Delta_Del_x2; y2 = legend_Delta_Del_y2; }
    else if (obv=="Gamma" && mode=="SSOS") { x1 = legend_Gamma_SSOS_x1; y1 = legend_Gamma_SSOS_y1; x2 = legend_Gamma_SSOS_x2; y2 = legend_Gamma_SSOS_y2; }
    else if (obv=="Gamma" && mode=="Del") { x1 = legend_Gamma_Del_x1; y1 = legend_Gamma_Del_y1; x2 = legend_Gamma_Del_x2; y2 = legend_Gamma_Del_y2; }
    else throw std::runtime_error("GetLegendBox – unsupported combination");
}

// ──────────────────────────────────────────────────────────────────────────
//  Globals (overridden per call)
// ──────────────────────────────────────────────────────────────────────────
TString ObvName  = "Gamma";   // ["Delta","Gamma"]
TString PlaneName= "TPC";     // kept for possible future use
TString SSOSDel  = "SSOS";    // ["SSOS","Del"]
int markerSize   = 2;

using namespace std;

// ───────────────────────────── Forward decls ─────────────────────────────
void SetStyle(Bool_t graypalette=kFALSE);
void fig_intg();
void plot_intg_compare2023();

// ───────────────────────── Data helpers ──────────────────────────────
struct PlotData {
    string particle; double centrality; string pair_type;
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
        getline(ss,p.particle,',');
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

TH1D* CreateHist(const vector<PlotData>& data,const string& part,const string& pair,
                 const string& obs,const string& err,double shift=0)
{
    vector<double> cx, val, errv;
    for(const auto& p:data) if(p.particle==part && p.pair_type==pair){
        cx.push_back(p.centrality);
        if(obs=="delta"){ val.push_back(p.delta);
            errv.push_back(err=="stat"?p.delta_err:p.delta_syst_err);}
        else {            val.push_back(p.gamma);
            errv.push_back(err=="stat"?p.gamma_err:p.gamma_syst_err);}
    }
    if(cx.empty()) return nullptr;
    string n=part+"_"+pair+"_"+obs+"_"+err;
    auto* h=new TH1D(n.c_str(),n.c_str(),(int)cx.size(),0,60);
    for(size_t i=0;i<cx.size();++i){
        int b=h->FindBin(cx[i]+shift);
        h->SetBinContent(b,val[i]); h->SetBinError(b,errv[i]); }
    return h;
}

// Style setting helpers for histograms
template<class TH> void SetStyle(TH& h, unsigned col, int mk, double ms, int ls=1, double lw=1){
    h->SetLineColor(col); h->SetMarkerColor(col);
    h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->SetLineStyle(ls);   h->SetLineWidth(lw);
}

// Systemic error styling
template<class TH> void SetStyle(TH& hS, TH& h){
    hS->SetFillColorAlpha(h->GetLineColor(), alpha_band);
    hS->SetFillStyle(fill_style);
    hS->SetMarkerStyle(20);
    hS->SetMarkerSize(0);
}

// ───────────────────────── Main plotting routine ─────────────────────────
void fig_intg()
{
    using FrameConfig::get;
    const auto& rng = FrameConfig::get(ObvName,SSOSDel);

    SetStyle();
    TString fname="fig_intg_compare2023_"+ObvName+"_"+PlaneName+"_"+SSOSDel+".pdf";
    auto* c=new TCanvas(fname,fname,canvas_width,canvas_height); c->SetTickx(1);
    TH1* frame=c->DrawFrame(rng.xmin,rng.ymin,rng.xmax,rng.ymax);
    if(!frame){cerr<<"DrawFrame failed\n";return;}

    frame->SetXTitle("Centrality (%)");
    if(ObvName=="Delta"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#delta} = #LTcos(#it{#phi}_{#alpha} - #it{#phi}_{#beta})#GT");
    if(ObvName=="Delta"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#delta} = #it{#delta}_{OS} - #it{#delta}_{SS}");
    if(ObvName=="Gamma"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#gamma} = #LTcos(#it{#phi}_{#alpha} + #it{#phi}_{#beta} - 2#Psi_{2})#GT");
    if(ObvName=="Gamma"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#gamma} = #it{#gamma}_{OS} - #it{#gamma}_{SS}");
    frame->GetYaxis()->SetMaxDigits(3);

    // Read 2025 data
    auto data_2025=ReadCSVData("../csv_data_point/intg_exp.csv");
    string obs=(ObvName=="Gamma")?"gamma":"delta";

    // 2025 Lambda–proton
    auto *hSS_2025=CreateHist(data_2025,"Proton","SS",obs,"stat"),
         *hOS_2025=CreateHist(data_2025,"Proton","OS",obs,"stat"),
         *hDel_2025=CreateHist(data_2025,"Proton","Del",obs,"stat"),
         *sSS_2025=CreateHist(data_2025,"Proton","SS",obs,"syst"),
         *sOS_2025=CreateHist(data_2025,"Proton","OS",obs,"syst"),
         *sDel_2025=CreateHist(data_2025,"Proton","Del",obs,"syst");

    // 2025 Lambda–Hadron
    auto *hSS_LH_2025=CreateHist(data_2025,"Hadron","SS",obs,"stat"),
         *hOS_LH_2025=CreateHist(data_2025,"Hadron","OS",obs,"stat"),
         *hDel_LH_2025=CreateHist(data_2025,"Hadron","Del",obs,"stat"),
         *sSS_LH_2025=CreateHist(data_2025,"Hadron","SS",obs,"syst"),
         *sOS_LH_2025=CreateHist(data_2025,"Hadron","OS",obs,"syst"),
         *sDel_LH_2025=CreateHist(data_2025,"Hadron","Del",obs,"syst");

    // 2023 data - Reading from ROOT files like fig_intg.C
    TH1D* hSS_2023 = nullptr, *hOS_2023 = nullptr, *hDel_2023 = nullptr;
    TH1D* sSS_2023 = nullptr, *sOS_2023 = nullptr, *sDel_2023 = nullptr;
    TH1D* hSS_LH_2023 = nullptr, *hOS_LH_2023 = nullptr, *hDel_LH_2023 = nullptr;
    TH1D* sSS_LH_2023 = nullptr, *sOS_LH_2023 = nullptr, *sDel_LH_2023 = nullptr;

    // Read 2023 data from ROOT files (like fig_intg.C)
    TFile* inputFile_2023 = new TFile("finalResult.root");
    if (inputFile_2023 && !inputFile_2023->IsZombie()) {
        hSS_2023 = (TH1D*)inputFile_2023->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_SS_DEta"  + "_default_integral");
        hOS_2023 = (TH1D*)inputFile_2023->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_OS_DEta"  + "_default_integral");
        hDel_2023 = (TH1D*)inputFile_2023->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_Del_DEta"  + "_default_integral");
        sSS_2023 = (TH1D*)inputFile_2023->Get("h"+ObvName+"SSSystError_" + PlaneName + "_DEta" );
        sOS_2023 = (TH1D*)inputFile_2023->Get("h"+ObvName+"OSSystError_" + PlaneName + "_DEta" );
        sDel_2023 = (TH1D*)inputFile_2023->Get("h"+ObvName+"DelSystError_" + PlaneName + "_DEta" );
    }

    TFile* inputFile_LambdaHadron_2023 = new TFile("finalResult_LambdaHadron.root");
    if (inputFile_LambdaHadron_2023 && !inputFile_LambdaHadron_2023->IsZombie()) {
        hSS_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_SS_"  + "default_integral");
        hOS_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_OS_"  + "default_integral");
        hDel_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_Del_"  + "default_integral");
        sSS_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("hLambdaHadron_"+ObvName+"SSSystError_" + PlaneName);
        sOS_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("hLambdaHadron_"+ObvName+"OSSystError_" + PlaneName);
        sDel_LH_2023 = (TH1D*)inputFile_LambdaHadron_2023->Get("hLambdaHadron_"+ObvName+"DelSystError_" + PlaneName);
    }

    /* Apply styles for 2025 data */
    if(hSS_2025){   SetStyle(hSS_2025 ,color_2025_LambdaProton ,marker_2025_LP_SS   ,markerSize); if(sSS_2025)   SetStyle(sSS_2025 ,hSS_2025 ); }
    if(hOS_2025){   SetStyle(hOS_2025 ,color_2025_LambdaProton ,marker_2025_LP_OS   ,markerSize); if(sOS_2025)   SetStyle(sOS_2025 ,hOS_2025 ); }
    if(hDel_2025){  SetStyle(hDel_2025,color_2025_LambdaProton ,marker_2025_LP_Del  ,markerSize); if(sDel_2025)  SetStyle(sDel_2025,hDel_2025); }

    if(hSS_LH_2025){SetStyle(hSS_LH_2025 ,color_2025_LambdaHadron,marker_2025_LH_SS ,markerSize); if(sSS_LH_2025)SetStyle(sSS_LH_2025,hSS_LH_2025);}
    if(hOS_LH_2025){SetStyle(hOS_LH_2025 ,color_2025_LambdaHadron,marker_2025_LH_OS ,markerSize); if(sOS_LH_2025)SetStyle(sOS_LH_2025,hOS_LH_2025);}
    if(hDel_LH_2025){SetStyle(hDel_LH_2025,color_2025_LambdaHadron,marker_2025_LH_Del ,markerSize); if(sDel_LH_2025)SetStyle(sDel_LH_2025,hDel_LH_2025);}

    /* Apply styles for 2023 data */
    if(hSS_2023){   SetStyle(hSS_2023 ,color_2023_LambdaProton ,marker_2023_LP_SS   ,markerSize); if(sSS_2023)   SetStyle(sSS_2023 ,hSS_2023 ); }
    if(hOS_2023){   SetStyle(hOS_2023 ,color_2023_LambdaProton ,marker_2023_LP_OS   ,markerSize); if(sOS_2023)   SetStyle(sOS_2023 ,hOS_2023 ); }
    if(hDel_2023){  SetStyle(hDel_2023,color_2023_LambdaProton ,marker_2023_LP_Del  ,markerSize); if(sDel_2023)  SetStyle(sDel_2023,hDel_2023); }

    if(hSS_LH_2023){SetStyle(hSS_LH_2023 ,color_2023_LambdaHadron,marker_2023_LH_SS ,markerSize); if(sSS_LH_2023)SetStyle(sSS_LH_2023,hSS_LH_2023);}
    if(hOS_LH_2023){SetStyle(hOS_LH_2023 ,color_2023_LambdaHadron,marker_2023_LH_OS ,markerSize); if(sOS_LH_2023)SetStyle(sOS_LH_2023,hOS_LH_2023);}
    if(hDel_LH_2023){SetStyle(hDel_LH_2023,color_2023_LambdaHadron,marker_2023_LH_Del ,markerSize); if(sDel_LH_2023)SetStyle(sDel_LH_2023,hDel_LH_2023);}

    // ====== 这里做ShiftAxis ======
    // 2023年结果向左移动-1.5，2025年向右移动1.5
    double dx_2023 = -1.5; // 2023年左移
    double dx_2025 = 1.5;  // 2025年右移
    
    // SSOS模式
    if(SSOSDel=="SSOS"){
        // 2023年数据
        ShiftAxis(hSS_2023, dx_2023); ShiftAxis(sSS_2023, dx_2023);
        ShiftAxis(hOS_2023, dx_2023); ShiftAxis(sOS_2023, dx_2023);
        ShiftAxis(hSS_LH_2023, dx_2023); ShiftAxis(sSS_LH_2023, dx_2023);
        ShiftAxis(hOS_LH_2023, dx_2023); ShiftAxis(sOS_LH_2023, dx_2023);
        
        // 2025年数据
        ShiftAxis(hSS_2025, dx_2025); ShiftAxis(sSS_2025, dx_2025);
        ShiftAxis(hOS_2025, dx_2025); ShiftAxis(sOS_2025, dx_2025);
        ShiftAxis(hSS_LH_2025, dx_2025); ShiftAxis(sSS_LH_2025, dx_2025);
        ShiftAxis(hOS_LH_2025, dx_2025); ShiftAxis(sOS_LH_2025, dx_2025);
    }else{
        // 2023年数据
        ShiftAxis(hDel_2023, dx_2023); ShiftAxis(sDel_2023, dx_2023);
        ShiftAxis(hDel_LH_2023, dx_2023); ShiftAxis(sDel_LH_2023, dx_2023);
        
        // 2025年数据
        ShiftAxis(hDel_2025, dx_2025); ShiftAxis(sDel_2025, dx_2025);
        ShiftAxis(hDel_LH_2025, dx_2025); ShiftAxis(sDel_LH_2025, dx_2025);
    }

    /* draw */
    c->cd();
    if(SSOSDel=="SSOS"){
        // Draw systematic error bands first (2023)
        if(sSS_2023) sSS_2023->Draw("E2,same"); if(sOS_2023) sOS_2023->Draw("E2,same");
        if(sSS_LH_2023) sSS_LH_2023->Draw("E2,same"); if(sOS_LH_2023) sOS_LH_2023->Draw("E2,same");
        
        // Draw systematic error bands (2025)
        if(sSS_2025) sSS_2025->Draw("E2,same"); if(sOS_2025) sOS_2025->Draw("E2,same");
        if(sSS_LH_2025) sSS_LH_2025->Draw("E2,same"); if(sOS_LH_2025) sOS_LH_2025->Draw("E2,same");

        // Draw data points (2023)
        if(hSS_2023) hSS_2023->Draw("E X0,same"); if(hOS_2023) hOS_2023->Draw("E X0,same");
        if(hSS_LH_2023) hSS_LH_2023->Draw("E X0,same"); if(hOS_LH_2023) hOS_LH_2023->Draw("E X0,same");
        
        // Draw data points (2025)
        if(hSS_2025) hSS_2025->Draw("E X0,same"); if(hOS_2025) hOS_2025->Draw("E X0,same");
        if(hSS_LH_2025) hSS_LH_2025->Draw("E X0,same"); if(hOS_LH_2025) hOS_LH_2025->Draw("E X0,same");
    }else{
        // Draw systematic error bands first (2023)
        if(sDel_2023) sDel_2023->Draw("E2,same");
        if(sDel_LH_2023) sDel_LH_2023->Draw("E2,same");
        
        // Draw systematic error bands (2025)
        if(sDel_2025) sDel_2025->Draw("E2,same");
        if(sDel_LH_2025) sDel_LH_2025->Draw("E2,same");

        // Draw data points (2023)
        if(hDel_2023) hDel_2023->Draw("E X0,same");
        if(hDel_LH_2023) hDel_LH_2023->Draw("E X0,same");
        
        // Draw data points (2025)
        if(hDel_2025) hDel_2025->Draw("E X0,same");
        if(hDel_LH_2025) hDel_LH_2025->Draw("E X0,same");
    }

    /* ALICE Preliminary tag */
    auto* pt=new TPaveText(alice_x1,alice_y1,alice_x2,alice_y2,"brNDC");
    pt->SetFillColor(0); pt->SetBorderSize(0); pt->SetShadowColor(0);
    pt->SetTextSize(gStyle->GetTextSize() * alice_text_size); pt->SetTextAlign(12);
    pt->AddText("ALICE CVE 2023 vs 2025");
    pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pt->Draw();

    /* Legend positioned using global settings */
    double leg_x1, leg_y1, leg_x2, leg_y2;
    GetLegendBox(ObvName,SSOSDel,leg_x1,leg_y1,leg_x2,leg_y2);
    auto* leg=new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
    leg->SetTextSize(legend_text_size*gStyle->GetTextSize());
    if(SSOSDel=="SSOS"){ leg->SetNColumns(2);
        // 2023 data
        if(hOS_2023) leg->AddEntry(hOS_2023,"#Lambda#minus#bar{p} 2023","PE");
        if(hSS_2023) leg->AddEntry(hSS_2023,"#Lambda#minusp 2023","PE");
        if(hOS_LH_2023) leg->AddEntry(hOS_LH_2023,"#Lambda#minush 2023","PE");
        if(hSS_LH_2023) leg->AddEntry(hSS_LH_2023,"#Lambda#minush 2023","PE");
        
        // 2025 data
        if(hOS_2025) leg->AddEntry(hOS_2025,"#Lambda#minus#bar{p} 2025","PE");
        if(hSS_2025) leg->AddEntry(hSS_2025,"#Lambda#minusp 2025","PE");
        if(hOS_LH_2025) leg->AddEntry(hOS_LH_2025,"#Lambda#minush 2025","PE");
        if(hSS_LH_2025) leg->AddEntry(hSS_LH_2025,"#Lambda#minush 2025","PE");
    }else{
        if(hDel_2023) leg->AddEntry(hDel_2023,"#Lambda#minusp 2023","PE");
        if(hDel_LH_2023) leg->AddEntry(hDel_LH_2023,"#Lambda#minush 2023","PE");
        if(hDel_2025) leg->AddEntry(hDel_2025,"#Lambda#minusp 2025","PE");
        if(hDel_LH_2025) leg->AddEntry(hDel_LH_2025,"#Lambda#minush 2025","PE");
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

// ───────────────────────── Driver to build all four canvases ─────────────
void plot_intg_compare2023()
{
    cout<<"Creating Delta SSOS comparison ...\n"; ObvName="Delta"; SSOSDel="SSOS"; fig_intg();
    cout<<"Creating Delta Del comparison ...\n";  ObvName="Delta"; SSOSDel="Del";  fig_intg();
    cout<<"Creating Gamma SSOS comparison ...\n"; ObvName="Gamma"; SSOSDel="SSOS"; fig_intg();
    cout<<"Creating Gamma Del comparison ...\n";  ObvName="Gamma"; SSOSDel="Del";  fig_intg();
    cout<<"All comparison plots completed!\n";
}
