#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

void comparePlane() {
  TFile* f[2];
  f[0] = TFile::Open("../train_output/Final_newPtCut/AnalysisResults_CVE2025_18r_TPC_proton.root");
  f[1] = TFile::Open("../train_output/TPC_RC/AnalysisResults_CVE2025_18q_TPC_Proton.root");

  auto getHist = [](TFile* file) -> TH2D* {
      if (!file) return nullptr;
      TList* list = (TList*)file->Get("default/ListResults_default");
      if (!list) return nullptr;
      return (TH2D*)list->FindObject("fHist2Psi2");
  };

  auto h2NM = getHist(f[0]);
  auto h2RC = getHist(f[1]);

  int n = 0;
  auto profile = [&n](TH2D * h2, int bin) -> TH1D* {
    if(h2) {
      return h2->ProjectionY(Form("hist_%d_%d",bin,n++), bin, bin);
    } else {
      return nullptr;
    }
  };

  // 创建画布
  auto c = new TCanvas("c", "TPC plane distribution", 800, 400);
  c->Divide(4, 2);

  // 中心度范围数组
  const char* centRanges[] = {
    "0-10%", "10-20%", "20-30%", "30-40%",
    "40-50%", "50-60%", "60-70%"
  };

  auto legend = new TLegend(0.1, 0.1, 0.9, 0.9);
  legend->SetBorderSize(0);


  // 创建并绘制7个子图
  for(int i = 1; i <= 7; i++) {
    c->cd(i);
    
    // 创建投影
    auto h1NM = profile(h2NM, i);
    auto h1RC = profile(h2RC, i);
    
    // 归一化
    if(h1NM) h1NM->Scale(h1NM->GetNbinsX()/h1NM->Integral());
    if(h1RC) h1RC->Scale(h1RC->GetNbinsX()/h1RC->Integral());
    
    // 设置颜色
    if(h1NM) h1NM->SetLineColor(kBlue);
    if(h1RC) h1RC->SetLineColor(kRed);
    
    // 设置线宽
    if(h1NM) h1NM->SetLineWidth(1);
    if(h1RC) h1RC->SetLineWidth(1);
    
    // 绘制框架
    gPad->DrawFrame(0, 0.85, TMath::Pi(), 1.15, 
                   Form("LHC18q (pass3) %s;#Psi_{2};Counts", centRanges[i-1]));
    
    // 绘制直方图
    if(h1NM) h1NM->Draw("same Hist");
    if(h1RC) h1RC->Draw("same Hist");
    
    if(i == 1) {
      legend->AddEntry(h1NM, "Before calibration", "l");
      legend->AddEntry(h1RC, "After calibration", "l");
    }
  }
  c->cd(8);
  legend->Draw();


  
  // 保存图片
  c->SaveAs("plane_comparison.pdf");
}
