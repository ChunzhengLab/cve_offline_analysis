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

//[DEta][SPt]
TString KineName = "SPt";
//[Delta][Gamma]
TString ObvName = "Delta";
//[TPC][V0C]
TString PlaneName = "TPC";
//[SSOS][Del]
TString SSOSDel = "Del";

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


const int markerFullStyle[] = {kFullCircle, kFullSquare, kFullCrossX};
const int markerOpenStyle[] = {kOpenCircle, kOpenSquare, kOpenCrossX};

vector<TString> vecKineSPtBins = {"SPtkineBin0", "SPtkineBin1", "SPtkineBin2"};
vector<TString> vecKineDEtaBins = {"DEtakineBin0", "DEtakineBin1"};

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
  double dx = 1;
  TAxis* a1 = h1->GetXaxis();
  a1->Set(a1->GetNbins(), a1->GetXmin() - dx, a1->GetXmax() - dx);
  TAxis* a2 = h2->GetXaxis();
  a2->Set(a2->GetNbins(), a2->GetXmin() + dx, a2->GetXmax() + dx);
}
void ShiftAxis(TH1D* h1, TH1D* h2, TH1D* h3) {
  // calculate the "x axis split"
  double dx = 1.5;
  TAxis* a1 = h1->GetXaxis();
  a1->Set(a1->GetNbins(), a1->GetXmin() - dx, a1->GetXmax() - dx);
  TAxis* a3 = h3->GetXaxis();
  a3->Set(a3->GetNbins(), a3->GetXmin() + dx, a3->GetXmax() + dx);
}

void fig_diff() {
  // 颜色
  int ci[3];
  ci[0] = TColor::GetColor("#1F77B4");  // 深蓝色（中心值）
  ci[1] = TColor::GetColor("#FF7F0E");  // 橙色（中心值）
  ci[2] = TColor::GetColor("#D62728");  // 红色（中心值）

  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TString figname = KineName + "_" + ObvName + "_" + PlaneName + "_" + SSOSDel;

  TCanvas *cfig = new TCanvas(figname, figname, 800, 600); 
  //抗锯齿
  cfig->SetTickx(1);
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  //cfig->SetLogy();
  // Set Titles etc..
  TH1 * h;
  if(ObvName == "Delta" && SSOSDel == "SSOS") {
    if(KineName == "SPt") h = cfig->DrawFrame(0,-0.0195,60,0.0195);
    if(KineName == "DEta") h = cfig->DrawFrame(0,-0.012,60,0.012);
  }
  if(ObvName == "Delta" && SSOSDel == "Del") {
    if(KineName == "SPt") h = cfig->DrawFrame(0,-0.002,60,0.022);
    if(KineName == "DEta") h = cfig->DrawFrame(0,-0.001,60,0.0145);
  }
  if(ObvName == "Gamma" && SSOSDel == "SSOS") {
    if(KineName == "SPt") h = cfig->DrawFrame(0,-0.007,60,0.003);
    if(KineName == "DEta") h = cfig->DrawFrame(0,-0.004,60,0.001);
  }
  if(ObvName == "Gamma" && SSOSDel == "Del") {
    if(KineName == "SPt") h = cfig->DrawFrame(0,-0.00085,60,0.0085);
    if(KineName == "DEta") h = cfig->DrawFrame(0,-0.00032,60,0.0032);
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
  h->GetYaxis()->SetTitleSize(gStyle->GetTextSize());
  h->GetXaxis()->SetTitleSize(gStyle->GetTextSize());



  TFile* inputFile = new TFile("finalResult_diff.root");
  if (!inputFile) {
    cout << "Cannot open file" << endl;
    return;
  }

  int nkineBins = 0;
  vector<TString> vecKineBins;
  
  if (KineName.EqualTo("SPt")) {
    nkineBins = 3;
    vecKineBins.assign(vecKineSPtBins.begin(), vecKineSPtBins.end());
  } else if (KineName.EqualTo("DEta")) {
    nkineBins = 2;
    vecKineBins.assign(vecKineDEtaBins.begin(), vecKineDEtaBins.end());
  } else {
    cout << "Wrong KineName" << endl;
    return;
  }
  TH1D* hObvSS[nkineBins];
  TH1D* hObvOS[nkineBins];
  TH1D* hObvDel[nkineBins];

  TH1D* hSysSS[nkineBins];
  TH1D* hSysOS[nkineBins];
  TH1D* hSysDel[nkineBins];

  for (int i = 0; i < nkineBins; i++) {
    hObvSS[i] = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_SS_" + vecKineBins[i] + "_default");
    if (!hObvSS[i]) {
      cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_SS_" + vecKineBins[i] + "_default" << endl;
      return;
    }
    hObvOS[i] = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_OS_" + vecKineBins[i] + "_default");
    if (!hObvOS[i]) {
      cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_OS_" + vecKineBins[i] + "_default" << endl;
      return;
    }
    hObvDel[i] = (TH1D*)inputFile->Get("h_18qr_"+ ObvName + "_" + PlaneName + "_Del_" + vecKineBins[i] + "_default");
    if (!hObvDel[i]) {
      cout << "Cannot find h_18qr_"+ ObvName + "_" + PlaneName + "_Del_" + vecKineBins[i] + "_default" << endl;
      return;
    }

    hSysSS[i] = (TH1D*)inputFile->Get("h"+ObvName+"SSSystError_" + PlaneName + "_" + vecKineBins[i]);
    if (!hSysSS[i]) {
      cout << "Cannot find h"+ObvName+"SSSystError_" + PlaneName + "_" + vecKineBins[i] << endl;
      return;
    }
    hSysOS[i] = (TH1D*)inputFile->Get("h"+ObvName+"OSSystError_" + PlaneName + "_" + vecKineBins[i]);
    if (!hSysOS[i]) {
      cout << "Cannot find h"+ObvName+"OSSystError_" + PlaneName + "_" + vecKineBins[i] << endl;
      return;
    }
    hSysDel[i] = (TH1D*)inputFile->Get("h"+ObvName+"DelSystError_" + PlaneName + "_" + vecKineBins[i]);
    if (!hSysDel[i]) {
      cout << "Cannot find h"+ObvName+"DelSystError_" + PlaneName + "_" + vecKineBins[i] << endl;
      return;
    }
  }

  ///////Shift Axis////////
  if (KineName.EqualTo("SPt")) {
    ShiftAxis(hObvSS[0],hObvSS[1],hObvSS[2]);
    ShiftAxis(hObvOS[0],hObvOS[1],hObvOS[2]);
    ShiftAxis(hSysSS[0],hSysSS[1],hSysSS[2]);
    ShiftAxis(hSysOS[0],hSysOS[1],hSysOS[2]);
    ShiftAxis(hObvDel[0],hObvDel[1],hObvDel[2]);
    ShiftAxis(hSysDel[0],hSysDel[1],hSysDel[2]);
  } else if (KineName.EqualTo("DEta")) {
    ShiftAxis(hObvSS[0],hObvSS[1]);
    ShiftAxis(hObvOS[0],hObvOS[1]);
    ShiftAxis(hSysSS[0],hSysSS[1]);
    ShiftAxis(hSysOS[0],hSysOS[1]);
    ShiftAxis(hObvDel[0],hObvDel[1]);
    ShiftAxis(hSysDel[0],hSysDel[1]);
  }
  cout<<"debug"<<endl;

  int colorIndex = -999;
  int markerIndex = -999;

  for (int i = 0; i < nkineBins; i++) {
    if      (vecKineBins[i] == "SPtkineBin0")  colorIndex = 0, markerIndex = 0;
    else if (vecKineBins[i] == "SPtkineBin1")  colorIndex = 1, markerIndex = 1;
    else if (vecKineBins[i] == "SPtkineBin2")  colorIndex = 2, markerIndex = 2;
    else if (vecKineBins[i] == "DEtakineBin0") colorIndex = 0, markerIndex = 0;
    else if (vecKineBins[i] == "DEtakineBin1") colorIndex = 1, markerIndex = 1;
    else {
      cout << "Error: kineBin is not in the list!" << endl;
      return;
    }

    SetStyle(hObvSS[i], ci[colorIndex], markerOpenStyle[markerIndex], markerSize, 1);
    SetStyle(hSysSS[i], hObvSS[i]);

    SetStyle(hObvOS[i], ci[colorIndex], markerFullStyle[markerIndex], markerSize, 1);
    SetStyle(hSysOS[i], hObvOS[i]);

    SetStyle(hObvDel[i], ci[colorIndex], markerFullStyle[markerIndex], markerSize, 1);
    SetStyle(hSysDel[i], hObvDel[i]);
  }

  if (SSOSDel.EqualTo("SSOS")) {
    for (int i = 0; i < nkineBins; i++) {
      hSysSS[i]->Draw("E2,same");
      hSysOS[i]->Draw("E2,same");
      hObvSS[i]->Draw("E X0,same");
      hObvOS[i]->Draw("E X0,same");
    }
  } else if (SSOSDel.EqualTo("Del")) {
    for (int i = 0; i < nkineBins; i++) {
      hSysDel[i]->Draw("E2,same");
      hObvDel[i]->Draw("E X0,same");
    }
  } else {
    cout << "Wrong SSOSDel" << endl;
    return;
  }



  TPaveText *pt = new TPaveText(0.3, 0.6, 0.5, 0.85, "brNDC");
  pt->SetTextSize(gStyle->GetTextSize());
  pt->SetTextAlign(12);
  pt->SetFillColor(0);
  pt->SetBorderSize(0);  // 设置边框大小为0，移除边框
  pt->SetShadowColor(0);  // 设置阴影颜色为0（透明），移除阴影
  pt->AddText("ALICE Preliminary");
  pt->AddText("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (ObvName == "Delta") {
    pt ->AddText("#Lambda#minusp correlation");
  } else if (ObvName == "Gamma") {
    pt ->AddText("#Lambda#minusp correlation");
  }
  pt->Draw();
  

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  //DrawLogo(1, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  //TLatex * text = new TLatex (0.5,0.5,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  //text->SetNDC();
  //text->Draw();
  //TLatex * text2 = new TLatex (0.55,55,"V0A Multiplicity Classes (Pb-Side)");
  //text2->SetTextSizePixels(24);
  //text2->Draw();
 
  //Legend, if needed
  TLegend * leg = new TLegend(0.19,  0.19,  0.57, 0.42);
  //leg->SetHeader("ALICE Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "C");
  TString legTitle;
  leg->SetTextSize(0.85 * gStyle->GetTextSize());
  if (SSOSDel.EqualTo("SSOS")) {
    leg->SetNColumns(2);
    for (int i = 0; i < nkineBins; i++) {
      if      (vecKineBins[i] == "SPtkineBin0")  legTitle = "1 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 3(GeV/#it{c})";
      else if (vecKineBins[i] == "SPtkineBin1")  legTitle = "3 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 5(GeV/#it{c})";
      else if (vecKineBins[i] == "SPtkineBin2")  legTitle = "5 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 8(GeV/#it{c})";
      else if (vecKineBins[i] == "DEtakineBin0") legTitle = "|#it{#eta}(#Lambda)#minus#it{#eta}(p)| < 0.6";
      else if (vecKineBins[i] == "DEtakineBin1") legTitle = "|#it{#eta}(#Lambda)#minus#it{#eta}(p)| > 0.6";
      leg->AddEntry(hObvOS[i], "OS   ",  "PE");
      leg->AddEntry(hObvSS[i], TString("SS   ") + legTitle,  "LPE");
    }
  } else if (SSOSDel.EqualTo("Del")) {
    for (int i = 0; i < nkineBins; i++) {
      if      (vecKineBins[i] == "SPtkineBin0")  legTitle = "1 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 3(GeV/#it{c})";
      else if (vecKineBins[i] == "SPtkineBin1")  legTitle = "3 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 5(GeV/#it{c})";
      else if (vecKineBins[i] == "SPtkineBin2")  legTitle = "5 < #it{p}_{T}(#Lambda)#plus#it{p}_{T}(p) < 8(GeV/#it{c})";
      else if (vecKineBins[i] == "DEtakineBin0") legTitle = "|#it{#eta} (#Lambda)#minus#it{#eta} (p)| < 0.6";
      else if (vecKineBins[i] == "DEtakineBin1") legTitle = "|#it{#eta} (#Lambda)#minus#it{#eta} (p)| > 0.6";
      leg->AddEntry(hObvDel[i], legTitle,  "PE");
    }
  } else {
    cout << "Wrong SSOSDel" << endl;
    return;
  }

  //leg->AddEntry(hstat,     "0-5\%, stat errors",   "LPE");
  //leg->AddEntry(hsyst,     "syst error (Uncorrelated)",  "F");

  //leg->AddEntry(hsystCorr, "syst error (Correlated)",    "F" );
  leg->SetFillColor(0);
  leg->Draw();

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

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
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