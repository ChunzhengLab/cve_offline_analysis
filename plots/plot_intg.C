/***************************************************************************
 *  Integrated correlator figure generator  (legend-position refactor)
 *  -----------------------------------------------------------------
 *  • In addition to FrameConfig (axis limits), a new LegendConfig
 *    namespace now holds the TLegend box coordinates for every
 *    (observable, mode) combination.
 *  • fig_intg() fetches its legend position via LegendConfig::get().
 *  • All styling configs are now placed at the top for easy access
 *  • Hadron-hadron results are now included
 *  • Marker styles corrected for Lambda-Lambda (circles) and Lambda-Proton (squares)
 *
 *  Usage:
 *      root -l -q plot_intg_refactor.C
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "RtypesCore.h"
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

// Color settings for different particle pairs
int color_LambdaProton = kRed+1;         // Lambda-Proton color
int color_LambdaLambda = kGreen+3;       // Lambda-Lambda color
int color_LambdaHadron = kBlue+1;        // Lambda-Hadron color
int color_HadronHadron_fill = kOrange+1; // Hadron-Hadron fill color
int color_HadronHadron_OS = kCyan+2;     // Hadron-Hadron opposite sign color

// Fill transparency
double alpha_HH = 0.5;                   // Transparency for Hadron-Hadron fills
int fill_style = 1000;                   // Fill style for error bands

// Marker settings
int marker_LP_SS = kOpenSquare;          // Lambda-Proton Same Sign marker
int marker_LP_OS = kFullSquare;          // Lambda-Proton Opposite Sign marker
int marker_LP_Del = kFullSquare;         // Lambda-Proton Delta marker

int marker_LL_SS = kOpenCircle;          // Lambda-Lambda SS marker
int marker_LL_OS = kFullCircle;          // Lambda-Lambda OS marker
int marker_LL_Del = kFullCircle;         // Lambda-Lambda Delta marker

int marker_LH_SS = kOpenTriangleUp;      // Lambda-Hadron SS marker
int marker_LH_OS = kFullTriangleDown;    // Lambda-Hadron OS marker
int marker_LH_Del = kFullTriangleUp;     // Lambda-Hadron Delta marker

// Drawing options
bool draw_HH_in_SSOS = false;            // Whether to draw hadron-hadron in SSOS mode

// ──────────────────────────────────────────────────────────────────────────
//  Marker and color reference guide (modify as needed)
// ──────────────────────────────────────────────────────────────────────────
// Lambda-Proton:
//   - OS (Lambda-pbar, Lambdabar-p): kRed+1, kFullSquare
//   - SS (Lambda-p, Lambdabar-pbar): kRed+1, kOpenCircle
//   - Del (OS-SS): kRed+1, kFullSquare
//
// Lambda-Lambda:
//   - OS (Lambda-Lambdabar): kGreen+3, kFullCircle
//   - SS (Lambda-Lambda, Lambdabar-Lambdabar): kGreen+3, kOpenCircle
//   - Del (OS-SS): kGreen+3, kFullCircle
//
// Lambda-Hadron:
//   - OS (Lambda-h-, Lambdabar-h+): kBlue+1, kFullTriangleDown
//   - SS (Lambda-h+, Lambdabar-h-): kBlue+1, kOpenTriangleUp
//   - Del (OS-SS): kBlue+1, kFullTriangleUp

// ──────────────────────────────────────────────────────────────────────────
//  Axis-range configuration
// ──────────────────────────────────────────────────────────────────────────
namespace FrameConfig {
struct Range { double xmin, ymin, xmax, ymax; };
constexpr Range Delta_SSOS { 0. , -0.015 , 60. ,  0.0150};
constexpr Range Delta_Del  { 0. , -0.002 , 60. ,  0.029};
constexpr Range Gamma_SSOS { 0. , -0.0038,  60. ,  0.0025};
constexpr Range Gamma_Del  { 0. , -0.0005, 60. ,  0.0048};

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
void plot_intg_refactor();

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

// No shift is applied to histograms on x-axis

// Style setting helpers for histograms
template<class TH> void SetStyle(TH& h, unsigned col, int mk, double ms, int ls=1, double lw=1){
    h->SetLineColor(col); h->SetMarkerColor(col);
    h->SetMarkerStyle(mk); h->SetMarkerSize(ms);
    h->SetLineStyle(ls);   h->SetLineWidth(lw);
}

// Systemic error styling
template<class TH> void SetStyle(TH& hS, TH& h){
    hS->SetFillColorAlpha(h->GetLineColor(), 0.2);
    hS->SetFillStyle(1000);
    hS->SetMarkerStyle(20);
    hS->SetMarkerSize(0);
}

// ───────────────────────── Main plotting routine ─────────────────────────
void fig_intg()
{
    using FrameConfig::get;
    const auto& rng = FrameConfig::get(ObvName,SSOSDel);

    SetStyle();
    TString fname="fig_intg_"+ObvName+"_"+PlaneName+"_"+SSOSDel+".pdf";
    auto* c=new TCanvas(fname,fname,canvas_width,canvas_height); c->SetTickx(1);
    TH1* frame=c->DrawFrame(rng.xmin,rng.ymin,rng.xmax,rng.ymax);
    if(!frame){cerr<<"DrawFrame failed\n";return;}

    frame->SetXTitle("Centrality (%)");
    if(ObvName=="Delta"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#delta} = #LTcos(#it{#phi}_{#alpha} - #it{#phi}_{#beta})#GT");
    if(ObvName=="Delta"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#delta} = #it{#delta}_{OS} - #it{#delta}_{SS}");
    if(ObvName=="Gamma"&&SSOSDel=="SSOS") frame->SetYTitle("#it{#gamma} = #LTcos(#it{#phi}_{#alpha} + #it{#phi}_{#beta} - 2#Psi_{2})#GT");
    if(ObvName=="Gamma"&&SSOSDel=="Del" ) frame->SetYTitle("#Delta#it{#gamma} = #it{#gamma}_{OS} - #it{#gamma}_{SS}");
    frame->GetYaxis()->SetMaxDigits(3);

    auto data=ReadCSVData("../csv_data_point/intg_exp.csv");
    string obs=(ObvName=="Gamma")?"gamma":"delta";

    // Lambda–proton
    auto *hSS=CreateHist(data,"Proton","SS",obs,"stat"),
         *hOS=CreateHist(data,"Proton","OS",obs,"stat"),
         *hDel=CreateHist(data,"Proton","Del",obs,"stat"),
         *sSS=CreateHist(data,"Proton","SS",obs,"syst"),
         *sOS=CreateHist(data,"Proton","OS",obs,"syst"),
         *sDel=CreateHist(data,"Proton","Del",obs,"syst");

    // Lambda–Lambda
    auto *hSS_LL=CreateHist(data,"Lambda","SS",obs,"stat"),
         *hOS_LL=CreateHist(data,"Lambda","OS",obs,"stat"),
         *hDel_LL=CreateHist(data,"Lambda","Del",obs,"stat"),
         *sSS_LL=CreateHist(data,"Lambda","SS",obs,"syst"),
         *sOS_LL=CreateHist(data,"Lambda","OS",obs,"syst"),
         *sDel_LL=CreateHist(data,"Lambda","Del",obs,"syst");

    // Lambda–Hadron
    auto *hSS_LH=CreateHist(data,"Hadron","SS",obs,"stat"),
         *hOS_LH=CreateHist(data,"Hadron","OS",obs,"stat"),
         *hDel_LH=CreateHist(data,"Hadron","Del",obs,"stat"),
         *sSS_LH=CreateHist(data,"Hadron","SS",obs,"syst"),
         *sOS_LH=CreateHist(data,"Hadron","OS",obs,"syst"),
         *sDel_LH=CreateHist(data,"Hadron","Del",obs,"syst");

    // Hadron–Hadron - Reading from external files
    TGraphAsymmErrors *hSS_HH = nullptr, *hOS_HH = nullptr, *hDel_HH = nullptr;
    TFile* inputFile_HH = nullptr;

    if (ObvName.EqualTo("Delta")) {
        inputFile_HH = new TFile("HEPData-ins1798528-v1-Table_1.root");
        if (inputFile_HH && !inputFile_HH->IsZombie()) {
            inputFile_HH->cd("Table 1");
            hSS_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y2");
            hOS_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y1");
            hDel_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y3");
        }
    } else if (ObvName.EqualTo("Gamma")) {
        inputFile_HH = new TFile("HEPData-ins1798528-v1-Table_5.root");
        if (inputFile_HH && !inputFile_HH->IsZombie()) {
            inputFile_HH->cd("Table 5");
            hSS_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y2");
            hOS_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y1");
            hDel_HH = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y3");
        }
    }

    // Remove the last point if HH graphs were successfully loaded
    if (hSS_HH && hOS_HH && hDel_HH) {
        hSS_HH->RemovePoint(hSS_HH->GetN() - 1);
        hOS_HH->RemovePoint(hOS_HH->GetN() - 1);
        hDel_HH->RemovePoint(hDel_HH->GetN() - 1);
    }

    /* Apply styles */
    if(hSS){   SetStyle(hSS ,color_LambdaProton ,marker_LP_SS   ,markerSize); if(sSS)   SetStyle(sSS ,hSS ); }
    if(hOS){   SetStyle(hOS ,color_LambdaProton ,marker_LP_OS   ,markerSize); if(sOS)   SetStyle(sOS ,hOS ); }
    if(hDel){  SetStyle(hDel,color_LambdaProton ,marker_LP_Del  ,markerSize); if(sDel)  SetStyle(sDel,hDel); }

    if(hSS_LL){SetStyle(hSS_LL ,color_LambdaLambda,marker_LL_SS ,markerSize); if(sSS_LL)SetStyle(sSS_LL,hSS_LL);}
    if(hOS_LL){SetStyle(hOS_LL ,color_LambdaLambda,marker_LL_OS ,markerSize); if(sOS_LL)SetStyle(sOS_LL,hOS_LL);}
    if(hDel_LL){SetStyle(hDel_LL,color_LambdaLambda,marker_LL_Del,markerSize); if(sDel_LL)SetStyle(sDel_LL,hDel_LL);}

    if(hSS_LH){SetStyle(hSS_LH ,color_LambdaHadron,marker_LH_SS ,markerSize); if(sSS_LH)SetStyle(sSS_LH,hSS_LH);}
    if(hOS_LH){SetStyle(hOS_LH ,color_LambdaHadron,marker_LH_OS ,markerSize); if(sOS_LH)SetStyle(sOS_LH,hOS_LH);}
    if(hDel_LH){SetStyle(hDel_LH,color_LambdaHadron,marker_LH_Del ,markerSize); if(sDel_LH)SetStyle(sDel_LH,hDel_LH);}

    // ====== 这里做ShiftAxis ======
    // LambdaProton不移动，LambdaHadron左移，LambdaLambda右移，HadronHadron不移动
    double dx_LH = -2.0; // LambdaHadron左移
    double dx_LL = 2.0;  // LambdaLambda右移
    // SSOS模式
    if(SSOSDel=="SSOS"){
        ShiftAxis(hSS_LH, dx_LH); ShiftAxis(sSS_LH, dx_LH);
        ShiftAxis(hOS_LH, dx_LH); ShiftAxis(sOS_LH, dx_LH);
        ShiftAxis(hSS_LL, dx_LL); ShiftAxis(sSS_LL, dx_LL);
        ShiftAxis(hOS_LL, dx_LL); ShiftAxis(sOS_LL, dx_LL);
    }else{
        ShiftAxis(hDel_LH, dx_LH); ShiftAxis(sDel_LH, dx_LH);
        ShiftAxis(hDel_LL, dx_LL); ShiftAxis(sDel_LL, dx_LL);
    }

    // HadronHadron不移动

    // Hadron-Hadron styling
    if(hSS_HH) {
        hSS_HH->SetFillColorAlpha(color_HadronHadron_fill, alpha_HH);
        hSS_HH->SetFillStyle(fill_style);
        hSS_HH->SetLineColorAlpha(color_HadronHadron_fill, alpha_HH);
    }
    if(hOS_HH) {
        hOS_HH->SetFillColorAlpha(color_HadronHadron_OS, alpha_HH);
        hOS_HH->SetFillStyle(fill_style);
        hOS_HH->SetLineColorAlpha(color_HadronHadron_OS, alpha_HH);
    }
    if(hDel_HH) {
        hDel_HH->SetFillColorAlpha(color_HadronHadron_fill, alpha_HH);
        hDel_HH->SetFillStyle(fill_style);
        hDel_HH->SetLineColorAlpha(color_HadronHadron_fill, alpha_HH);
    }

    /* draw */
    c->cd();
    if(SSOSDel=="SSOS"){
        // Only draw hadron-hadron if configured to do so
        if(draw_HH_in_SSOS) {
            if(hSS_HH) hSS_HH->Draw("E3,same");
            if(hOS_HH) hOS_HH->Draw("E3,same");
        }

        if(sSS) sSS->Draw("E2,same"); if(sOS) sOS->Draw("E2,same");
        if(sSS_LL) sSS_LL->Draw("E2,same"); if(sOS_LL) sOS_LL->Draw("E2,same");
        if(sSS_LH) sSS_LH->Draw("E2,same"); if(sOS_LH) sOS_LH->Draw("E2,same");

        if(hSS) hSS->Draw("E X0,same"); if(hOS) hOS->Draw("E X0,same");
        if(hSS_LL) hSS_LL->Draw("E X0,same"); if(hOS_LL) hOS_LL->Draw("E X0,same");
        if(hSS_LH) hSS_LH->Draw("E X0,same"); if(hOS_LH) hOS_LH->Draw("E X0,same");
    }else{
        if(hDel_HH) hDel_HH->Draw("E3,same");

        if(sDel) sDel->Draw("E2,same");
        if(sDel_LL) sDel_LL->Draw("E2,same");
        if(sDel_LH) sDel_LH->Draw("E2,same");

        if(hDel) hDel->Draw("E X0,same");
        if(hDel_LL) hDel_LL->Draw("E X0,same");
        if(hDel_LH) hDel_LH->Draw("E X0,same");
    }

    /* ALICE Preliminary tag */
    auto* pt=new TPaveText(alice_x1,alice_y1,alice_x2,alice_y2,"brNDC");
    pt->SetFillColor(0); pt->SetBorderSize(0); pt->SetShadowColor(0);
    pt->SetTextSize(gStyle->GetTextSize() * alice_text_size); pt->SetTextAlign(12);
    pt->AddText("ALICE CVE 2025");
    pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    pt->Draw();

    /* Legend positioned using global settings */
    double leg_x1, leg_y1, leg_x2, leg_y2;
    GetLegendBox(ObvName,SSOSDel,leg_x1,leg_y1,leg_x2,leg_y2);
    auto* leg=new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
    leg->SetTextSize(legend_text_size*gStyle->GetTextSize());
    if(SSOSDel=="SSOS"){ leg->SetNColumns(2);
        if(hOS_LL) leg->AddEntry(hOS_LL,"#Lambda#minus#bar{#Lambda}","PE");
        if(hSS_LL) leg->AddEntry(hSS_LL,"#Lambda#minus#Lambda#oplus#bar{#Lambda}#minus#bar{#Lambda}","PE");
        if(hOS) leg->AddEntry(hOS,"#Lambda#minus#bar{p}#oplus#bar{#Lambda}#minusp","PE");
        if(hSS) leg->AddEntry(hSS,"#Lambda#minusp#oplus#bar{#Lambda}#minus#bar{p}","PE");
        if(hOS_LH) leg->AddEntry(hOS_LH,"#Lambda#minush^{#minus}#oplus#bar{#Lambda}#minush^{#plus}","PE");
        if(hSS_LH) leg->AddEntry(hSS_LH,"#Lambda#minush^{#plus}#oplus#bar{#Lambda}#minush^{#minus}","PE");
        // Only add hadron-hadron to legend if configured to draw it
        if(draw_HH_in_SSOS) {
            if(hOS_HH) leg->AddEntry(hOS_HH,"h^{+}#minush^{#minus}","F");
            if(hSS_HH) leg->AddEntry(hSS_HH,"h^{+}#minush^{+} #oplus h^{#minus}#minush^{#minus}","F");
        }
    }else{
        if(hDel_HH) leg->AddEntry(hDel_HH,"h#minush [ALICE, JHEP 160 (2020)]","F");
        if(hDel) leg->AddEntry(hDel,"#Lambda#minusp","PE");
        if(hDel_LL) leg->AddEntry(hDel_LL,"#Lambda#minus#Lambda","PE");
        if(hDel_LH) leg->AddEntry(hDel_LH,"#Lambda#minush","PE");
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
void plot_intg()
{
    cout<<"Creating Delta SSOS ...\n"; ObvName="Delta"; SSOSDel="SSOS"; fig_intg();
    cout<<"Creating Delta Del ...\n";  ObvName="Delta"; SSOSDel="Del";  fig_intg();
    cout<<"Creating Gamma SSOS ...\n"; ObvName="Gamma"; SSOSDel="SSOS"; fig_intg();
    cout<<"Creating Gamma Del ...\n";  ObvName="Gamma"; SSOSDel="Del";  fig_intg();
    cout<<"All plots completed!\n";
}
