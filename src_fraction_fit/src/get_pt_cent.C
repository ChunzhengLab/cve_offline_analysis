#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TList.h"
#include "TGraph.h"

void get_pt_cent() {
    // 打开文件
    TFile* f18q = TFile::Open("../data/data_18q.root", "READ");
    TFile* f18r = TFile::Open("../data/data_18r.root", "READ");
    TFile* f_nua = TFile::Open("../../nua_pt/eff_pt_calib_cent.root", "READ");
    
    // 获取TList
    TList* list18q = (TList*)f18q->Get("default/ResultsList_default");
    TList* list18r = (TList*)f18r->Get("default/ResultsList_default");
    TList* nua_list = (TList*)f_nua->Get("fListNUENUA");
    
    // 获取直方图
    TH2D* h2_pt_proton_18q = (TH2D*)list18q->FindObject("h2_pt_proton");
    TH2D* h2_pt_antiproton_18q = (TH2D*)list18q->FindObject("h2_pt_antiproton");
    TH2D* h2_pt_lambda_18q = (TH2D*)list18q->FindObject("h2_pt_lambda");
    TH2D* h2_pt_antilambda_18q = (TH2D*)list18q->FindObject("h2_pt_antilambda");
    
    TH2D* h2_pt_proton_18r = (TH2D*)list18r->FindObject("h2_pt_proton");
    TH2D* h2_pt_antiproton_18r = (TH2D*)list18r->FindObject("h2_pt_antiproton");
    TH2D* h2_pt_lambda_18r = (TH2D*)list18r->FindObject("h2_pt_lambda");
    TH2D* h2_pt_antilambda_18r = (TH2D*)list18r->FindObject("h2_pt_antilambda");

    // 创建输出文件
    TFile* fout = new TFile("pt_dist.root", "RECREATE");
    
    // 粒子类型数组
    const char* particles[] = {"proton", "antiproton", "lambda", "antilambda"};
    // 中心度范围
    const char* centRanges[] = {"cent0", "cent1", "cent2", "cent3", "cent4", "cent5", "cent6"};
    
    // 处理18q数据
    TH2D* h2_18q[] = {h2_pt_proton_18q, h2_pt_antiproton_18q, h2_pt_lambda_18q, h2_pt_antilambda_18q};
    
    for(int i = 0; i < 4; i++) {  // 4种粒子
        for(int j = 1; j <= 7; j++) {  // 7个中心度bin
            if(h2_18q[i]) {
                TH1D* h1 = h2_18q[i]->ProjectionY(Form("h_pt_%s_18q_%s", particles[i], centRanges[j-1]), j, j);
                h1->Write();
                
                // 获取对应中心度的效率修正TGraph
                TGraph* eff_graph = (TGraph*)nua_list->FindObject(Form("nue_pt_proton_18q_%s", centRanges[j-1]));
                if(eff_graph) {
                    // 创建效率修正后的直方图
                    TH1D* h1_eff = (TH1D*)h1->Clone(Form("h_pt_%s_18q_%s_eff_col", particles[i], centRanges[j-1]));
                    
                    // 应用效率修正
                    for(int bin = 1; bin <= h1_eff->GetNbinsX(); bin++) {
                        double pt = h1_eff->GetBinCenter(bin);
                        double efficiency = eff_graph->Eval(pt);
                        h1_eff->SetBinContent(bin, h1_eff->GetBinContent(bin) * efficiency);
                    }
                    h1_eff->Write();
                } else {
                    cout << "Warning: Efficiency graph not found for " << Form("nue_pt_proton_18q_%s", centRanges[j-1]) << endl;
                }
            }
        }
    }
    
    // 处理18r数据
    TH2D* h2_18r[] = {h2_pt_proton_18r, h2_pt_antiproton_18r, h2_pt_lambda_18r, h2_pt_antilambda_18r};
    
    for(int i = 0; i < 4; i++) {  // 4种粒子
        for(int j = 1; j <= 7; j++) {  // 7个中心度bin
            if(h2_18r[i]) {
                TH1D* h1 = h2_18r[i]->ProjectionY(Form("h_pt_%s_18r_%s", particles[i], centRanges[j-1]), j, j);
                h1->Write();
                
                // 获取对应中心度的效率修正TGraph
                TGraph* eff_graph = (TGraph*)nua_list->FindObject(Form("nue_pt_proton_18r_%s", centRanges[j-1]));
                if(eff_graph) {
                    // 创建效率修正后的直方图
                    TH1D* h1_eff = (TH1D*)h1->Clone(Form("h_pt_%s_18r_%s_eff_col", particles[i], centRanges[j-1]));
                    
                    // 应用效率修正
                    for(int bin = 1; bin <= h1_eff->GetNbinsX(); bin++) {
                        double pt = h1_eff->GetBinCenter(bin);
                        double efficiency = eff_graph->Eval(pt);
                        h1_eff->SetBinContent(bin, h1_eff->GetBinContent(bin) * efficiency);
                    }
                    h1_eff->Write();
                } else {
                    cout << "Warning: Efficiency graph not found for " << Form("nue_pt_proton_18r_%s", centRanges[j-1]) << endl;
                }
            }
        }
    }
    
    // 关闭文件
    fout->Close();
    f18q->Close();
    f18r->Close();
    f_nua->Close();
    
    cout << "Created pt_dist.root with all projections and efficiency corrections" << endl;
}
