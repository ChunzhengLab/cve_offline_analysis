#include <iostream>
#include <string>

#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TStyle.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPaveText.h"


//[Delta][Gamma]
TString ObvName = "Gamma";
//[TPC][V0C]
TString PlaneName = "TPC";
//[SSOS][Del]
TString SSOSDel = "SSOS";

int markerSize = 2;


using namespace std;

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

// 设置样式
template <class TH>
void SetStyle(TH& hist, unsigned int color = kBlack, unsigned int markerStyle = kFullCircle, double markerSize = 1, int lineStyle = 1, double lineWidth = 1) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(markerSize);
  hist->SetLineStyle(lineStyle);
  hist->SetLineWidth(lineWidth);
}

template <class TH>
void SetStyle(TH& h_syst, TH& h_stat) {
  h_syst->SetFillColorAlpha(h_stat->GetLineColor(), 0.2);
  h_syst->SetFillStyle(1000);
  h_syst->SetMarkerStyle(20);
  h_syst->SetMarkerSize(0);
}

void ShiftAxis(TH1D* h1, TH1D* h2) {
  // calculate the "x axis split"
  double dx = 1.;
  TAxis* a1 = h1->GetXaxis();
  a1->Set(a1->GetNbins(), a1->GetXmin() - dx, a1->GetXmax() - dx);
  TAxis* a2 = h2->GetXaxis();
  a2->Set(a2->GetNbins(), a2->GetXmin() + dx, a2->GetXmax() + dx);
}

void fig_intg() {
  // 颜色
  int ci[3];
  ci[0] = TColor::GetColor("#1F77B4");  // 深蓝色（中心值）
  ci[1] = TColor::GetColor("#FF7F0E");  // 橙色（中心值）
  ci[2] = TColor::GetColor("#D62728");  // 红色（中心值）

  int ci_fill[3];
  ci_fill[2] = TColor::GetColor("#AEC7E8");  // 浅蓝色（系统误差）
  ci_fill[1] = TColor::GetColor("#FFBB78");  // 浅橙色（系统误差）
  ci_fill[0] = TColor::GetColor("#FAC8C3");  // 粉红色（系统误差）

  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TString figname = "fig_intg_" + ObvName + "_" + PlaneName + "_" + SSOSDel + ".pdf";
  TCanvas *cfig = new TCanvas(figname, figname, 800, 600); 
  //抗锯齿
  cfig->SetTickx(1);
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  //cfig->SetLogy();
  // Set Titles etc..
  TH1 * h;
  if(ObvName == "Delta" && SSOSDel == "SSOS") {
    h = cfig->DrawFrame(0,-0.009,60,0.009);
  }
  if(ObvName == "Delta" && SSOSDel == "Del") {
    h = cfig->DrawFrame(0,-0.002,60,0.0145);
  }
  if(ObvName == "Gamma" && SSOSDel == "SSOS") {
    h = cfig->DrawFrame(0,-0.0035,60,0.0005);
  }
  if(ObvName == "Gamma" && SSOSDel == "Del") {
    h = cfig->DrawFrame(0,-0.0003,60,0.003);
  }
  if(!h) {
    cout << "Error: No frame drawn" << endl;
    return;
  }

  // === Commonly used x/ titles: ===
  // pt invariant yields
  const char *  texPtX="#it{p}_{T} (GeV/#it{c})";
  const char *  texPtY="1/#it{N}_{ev} 1/(2#pi#it{p}_{T}) d#it{N}/(d#it{p}_{T}d#it{y}) ((GeV/#it{c})^{-2})";
  // mt invariant yields
  const char *  texMtX="#it{m}_{T} (GeV/#it{c}^{2})";
  const char *  texMtY="1/#it{N}_{ev} 1/(2#pi#it{m}_{T}) d#it{N}/(d#it{m}_{T}d#it{y}) ((GeV/#it{c}^{2})^{-2}) "; 
  // Invariant mass with decay products K and pi
  const char *  texMassX="#it{M}_{K#pi} (GeV/#it{c}^{2})";
  const char *  texMassY="d#it{N}/(d#it{M}_{K#pi})";
  // <pt>, npart
  const char * texMeanPt    = "#LT#it{p}_{T}#GT (GeV/#it{c})";
  const char * texMeanNpart = "#LT#it{N}_{part}#GT";

  // Centrality
  const char * texCentrality = "Centrality (%)";
  // delta for CVE
  const char * texDeltaSSOS = "#it{#delta} = #LTcos(#it{#phi}_{#alpha} - #it{#phi}_{#beta})#GT";
  const char * texGammaSSOS = "#it{#gamma} = #LTcos(#it{#phi}_{#alpha} + #it{#phi}_{#beta} - 2#Psi_{2})#GT";
  const char * texDeltaDel = "#Delta#it{#delta} = #it{#delta}_{OS} - #it{#delta}_{SS}";
  const char * texGammaDel = "#Delta#it{#gamma} = #it{#gamma}_{OS} - #it{#gamma}_{SS}";


  // Set titles
  h->SetXTitle(texCentrality);
  // Please be consistent on the y label
  if (ObvName == "Delta"&& SSOSDel == "SSOS") h->SetYTitle(texDeltaSSOS);
  if (ObvName == "Delta"&& SSOSDel == "Del") h->SetYTitle(texDeltaDel);
  if (ObvName == "Gamma"&& SSOSDel == "SSOS") h->SetYTitle(texGammaSSOS);
  if (ObvName == "Gamma"&& SSOSDel == "Del") h->SetYTitle(texGammaDel);
  h->GetYaxis()->SetMaxDigits(3);



  TFile* inputFile = new TFile("finalResult.root");
  if (!inputFile) {
    cout << "Cannot open file" << endl;
    return;
  }

  TFile* inputFile_LambdaHadron = new TFile("finalResult_LambdaHadron.root");
  if (!inputFile_LambdaHadron) {
    cout << "Cannot open file" << endl;
    return;
  }

  TH1D* hObvSS;
  TH1D* hObvOS;
  TH1D* hObvDel;

  TH1D* hSysSS;
  TH1D* hSysOS;
  TH1D* hSysDel;

  TH1D* hObvSS_LambdaHadron;
  TH1D* hObvOS_LambdaHadron;
  TH1D* hObvDel_LambdaHadron;

  TH1D* hSysSS_LambdaHadron;
  TH1D* hSysOS_LambdaHadron;
  TH1D* hSysDel_LambdaHadron;


  ////////////////Read in LambdaProton/////////////////////
  hObvSS = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_SS_DEta"  + "_default_integral");
  if (!hObvSS) {
    cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_SS_DEta"  + "_default" << endl;
    return;
  }
  hObvOS = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_OS_DEta"  + "_default_integral");
  if (!hObvOS) {
    cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_OS_DEta"  + "_default" << endl;
    return;
  }
  hObvDel = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_Del_DEta"  + "_default_integral");
  if (!hObvDel) {
    cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_Del_DEta"  + "_default" << endl;
    return;
  }

  hSysSS = (TH1D*)inputFile->Get("h"+ObvName+"SSSystError_" + PlaneName + "_DEta" );
  if (!hSysSS) {
    cout << "Cannot find h"+ObvName+"SSSystError_" + PlaneName + "_DEta"  << endl;
    return;
  }
  hSysOS = (TH1D*)inputFile->Get("h"+ObvName+"OSSystError_" + PlaneName + "_DEta" );
  if (!hSysOS) {
    cout << "Cannot find h"+ObvName+"OSSystError_" + PlaneName + "_DEta"  << endl;
    return;
  }
  hSysDel = (TH1D*)inputFile->Get("h"+ObvName+"DelSystError_" + PlaneName + "_DEta" );
  if (!hSysDel) {
    cout << "Cannot find h"+ObvName+"DelSystError_" + PlaneName + "_DEta"  << endl;
    return;
  }  


  ////////////////Read in LambdaHadron/////////////////////
  hObvSS_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_SS_"  + "default_integral");
  if (!hObvSS_LambdaHadron) {
    cout << "Cannot find h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_SS_"  + "default" << endl;
    return;
  }
  hObvOS_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_OS_"  + "default_integral");
  if (!hObvOS_LambdaHadron) {
    cout << "Cannot find h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_OS_"  + "default" << endl;
    return;
  }
  hObvDel_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_Del_"  + "default_integral");
  if (!hObvDel_LambdaHadron) {
    cout << "Cannot find h_LambdaHadron_18qr_"+ ObvName + "_" + PlaneName + "_Del_"  + "default" << endl;
    return;
  }

  hSysSS_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("hLambdaHadron_"+ObvName+"SSSystError_" + PlaneName);
  if (!hSysSS_LambdaHadron) {
    cout << "Cannot find hLambdaHadron_"+ObvName+"SSSystError_" + PlaneName << endl;
    return;
  }
  hSysOS_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("hLambdaHadron_"+ObvName+"OSSystError_" + PlaneName);
  if (!hSysOS_LambdaHadron) {
    cout << "Cannot find hLambdaHadron_"+ObvName+"OSSystError_" + PlaneName << endl;
    return;
  }
  hSysDel_LambdaHadron = (TH1D*)inputFile_LambdaHadron->Get("hLambdaHadron_"+ObvName+"DelSystError_" + PlaneName);
  if (!hSysDel_LambdaHadron) {
    cout << "Cannot find hLambdaHadron_"+ObvName+"DelSystError_" + PlaneName << endl;
    return;
  }

  ////////////////Read in HadronHadron/////////////////////
  TFile* inputFile_HadronHadron;
  if (ObvName.EqualTo("Delta")) {
   inputFile_HadronHadron = new TFile("HEPData-ins1798528-v1-Table_1.root");
  inputFile_HadronHadron->cd("Table 1");
  } else if (ObvName.EqualTo("Gamma")) {
    inputFile_HadronHadron = new TFile("HEPData-ins1798528-v1-Table_5.root");
    inputFile_HadronHadron->cd("Table 5");
  } else return;

  TGraphAsymmErrors* hObvSS_HadronHadron = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y2");
  if (!hObvSS_HadronHadron) {
    cout << "Cannot find Graph1D_y1" << endl;
    return;
  }
  TGraphAsymmErrors* hObvOS_HadronHadron = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y1");
  if (!hObvOS_HadronHadron) {
    cout << "Cannot find Graph1D_y2" << endl;
    return;
  }
  TGraphAsymmErrors* hObvDel_HadronHadron = (TGraphAsymmErrors*)gDirectory->Get("Graph1D_y3");
  if (!hObvDel_HadronHadron) {
    cout << "Cannot find Graph1D_y3" << endl;
    return;
  }
  hObvSS_HadronHadron->RemovePoint(hObvSS_HadronHadron->GetN() - 1);
  hObvOS_HadronHadron->RemovePoint(hObvOS_HadronHadron->GetN() - 1);
  hObvDel_HadronHadron->RemovePoint(hObvDel_HadronHadron->GetN() - 1);

  int colorIndex = -999;
  int markerIndex = -999;

  SetStyle(hObvSS, colors[1], kOpenCircle, markerSize, 1);
  SetStyle(hSysSS, hObvSS);

  SetStyle(hObvOS, colors[1], kFullSquare, markerSize, 1);
  SetStyle(hSysOS, hObvOS);

  SetStyle(hObvDel, colors[1], kFullSquare, markerSize, 1);
  SetStyle(hSysDel, hObvDel);

  SetStyle(hObvSS_LambdaHadron, colors[2], kOpenTriangleUp, markerSize, 1);
  SetStyle(hSysSS_LambdaHadron, hObvSS_LambdaHadron);

  SetStyle(hObvOS_LambdaHadron, colors[2], kFullTriangleDown, markerSize, 1);
  SetStyle(hSysOS_LambdaHadron, hObvOS_LambdaHadron);

  SetStyle(hObvDel_LambdaHadron, colors[2], kFullTriangleUp, markerSize, 1);
  SetStyle(hSysDel_LambdaHadron, hObvDel_LambdaHadron);

  ShiftAxis(hObvSS, hObvSS_LambdaHadron);
  ShiftAxis(hObvOS, hObvOS_LambdaHadron);
  ShiftAxis(hSysSS, hSysSS_LambdaHadron);
  ShiftAxis(hSysOS, hSysOS_LambdaHadron);
  ShiftAxis(hObvDel, hObvDel_LambdaHadron);
  ShiftAxis(hSysDel, hSysDel_LambdaHadron);

  hObvSS_HadronHadron->SetFillColorAlpha(kOrange +1, 0.5); // 设置浅蓝色，并使用alpha透明度参数设置透明度
  hObvSS_HadronHadron->SetFillStyle(1000);            // 设置填充样式
  hObvSS_HadronHadron->SetLineColorAlpha(kOrange +1, 0.5);            // 设置填充样式

  hObvOS_HadronHadron->SetFillColorAlpha(kCyan+2, 0.5); // 设置浅蓝色，并使用alpha透明度参数设置透明度
  hObvOS_HadronHadron->SetFillStyle(1000);            // 设置填充样式
  hObvOS_HadronHadron->SetLineColorAlpha(kCyan+2, 0.5);            // 设置填充样式

  hObvDel_HadronHadron->SetFillColorAlpha(kOrange +1, 0.5); // 设置浅蓝色，并使用alpha透明度参数设置透明度
  hObvDel_HadronHadron->SetFillStyle(1000);            // 设置填充样式
  hObvDel_HadronHadron->SetLineColorAlpha(kOrange +1, 0.5);            // 设置填充样式

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  //DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with

  
  cfig->cd();
  if (SSOSDel.EqualTo("SSOS")) {
    hSysSS->Draw("E2,same");
    hSysOS->Draw("E2,same");
    hSysOS_LambdaHadron->Draw("E2,same");
    hSysSS_LambdaHadron->Draw("E2,same");
    hObvSS->Draw("E X0,same");
    hObvOS->Draw("E X0,same");
    hObvOS_LambdaHadron->Draw("E X0,same");
    hObvSS_LambdaHadron->Draw("E X0,same");
  } else if (SSOSDel.EqualTo("Del")) {
    hObvDel->Draw("E X0,same");
    hSysDel->Draw("E2,same");
    hObvDel_LambdaHadron->Draw("E X0,same");
    hSysDel_LambdaHadron->Draw("E2,same");
    hObvDel_HadronHadron->Draw("E3,same");
  } else {
    cout << "Wrong SSOSDel" << endl;
    return;
  }

  // TLatex * texALIEC = new TLatex(0.4,0.4, "ALICE Preliminary");
  // texALIEC->SetNDC();
  // texALIEC->SetTextFont(42);
  // texALIEC->Draw();

  // TLatex * text = new TLatex (0.5,0.5,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // text->SetNDC();
  // text->SetTextFont(42);

  TPaveText *pt = new TPaveText(0.2, 0.65, 0.4, 0.85, "brNDC");
  pt->SetTextSize(gStyle->GetTextSize());
  pt->SetTextAlign(12);
  pt->SetFillColor(0);
  pt->SetBorderSize(0);  // 设置边框大小为0，移除边框
  pt->SetShadowColor(0);  // 设置阴影颜色为0（透明），移除阴影
  pt->AddText("ALICE Preliminary");
  pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  pt->Draw();

  //text->Draw();
  //TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  //text2->SetTextSizePixels(24);
  //text2->Draw();

  TLegend * leg = new TLegend(0.19,  0.19,  0.65, 0.4);
  leg->SetTextSize(0.9 * gStyle->GetTextSize());
  //leg->SetHeader("ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "C");
  if(SSOSDel == "SSOS") {
    leg->SetNColumns(2);
    leg->AddEntry(hObvOS, "#Lambda#minus#bar{p}  #oplus #bar{#Lambda}#minusp", "PE");
    leg->AddEntry(hObvSS, "#Lambda#minusp  #oplus #bar{#Lambda}#minus#bar{p}   ", "PE");
    leg->AddEntry(hObvOS_LambdaHadron, "#Lambda#minush^{#minus} #oplus #bar{#Lambda}#minush^{#plus}", "PE");
    leg->AddEntry(hObvSS_LambdaHadron, "#Lambda#minush^{#plus} #oplus #bar{#Lambda}#minush^{#minus}   ", "PE");
    // leg->AddEntry(hObvSS_HadronHadron, "h^{+}-h^{+} #oplus h^{-}-h^{-}  ", "F");
    // leg->AddEntry(hObvOS_HadronHadron, "h^{+}-h^{-}", "F");
  } else if(SSOSDel == "Del") {
    leg->AddEntry(hObvDel_HadronHadron, "h#minush [ALICE, JHEP 160 (2020)]", "F");
    leg->AddEntry(hObvDel, "#Lambda#minusp", "PE");
    leg->AddEntry(hObvDel_LambdaHadron, "#Lambda#minush", "PE");
  }
  // leg->AddEntry(, "#Lambda-h^{+} and #bar{#Lambda}-h^{-}", "LPE");
  // leg->AddEntry(, "#Lambda-h^{-} and #bar{#Lambda}-h^{+}", "LPE");
  leg->Draw();






  //leg->AddEntry(hstat,     "0-5\%, stat errors",   "LPE");
  //leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");

  //leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  // leg->SetFillColor(0);
  // leg->SetTexgStyle->GetTextSize()*0.8);
  
  // leg->Draw();

  // Save to HEP data

  // AliHEPDataParser * hepParser = new AliHEPDataParser(hstat, hsyst);
  // hepParser->SetTitle("pt distribution of pi+-, arXiv:XXXX.YYYY");
  // hepParser->SetName("1/Nev 1/p_T 1/2pi d^2N/(dp_Tdy) (GeV/c)^{-1}"); 
  // hepParser->SetXaxisName("PT IN GEV/c");
  // hepParser->SetReaction("RE: P PB --> PI + X");
  // hepParser->SetEnergy("SQRT(SNN) : 5020.0 GeV");
  // hepParser->SetRapidityRange("YRAP : -0.5 - +0.5");
  // hepParser->SaveHEPDataFile("figTemplateHEPData.txt");    // it must be specified explicity if graphs are to be used






}

//________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE maybe Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.048,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetTitleOffset(1,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetErrorX(0.075);
  gStyle->SetPadLeftMargin(0.09);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetNdivisions(505,"y");
}