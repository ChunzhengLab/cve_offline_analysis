// 绘制MC中不同source的DCA分布（pT积分）
void DrawMCDCAvsPtBins(
    const char* dataset = "18q",      // "18q" 或 "18r"
    const char* particle = "proton")   // "proton" 或 "antiproton"
{
    // 打开MC文件
    std::string mc_file = Form("mc_%s.root", dataset);
    TFile* fMC = TFile::Open(mc_file.c_str());
    if(!fMC || fMC->IsZombie()) { std::cerr << "Can't open MC file: " << mc_file << std::endl; return; }
    TDirectoryFile* dirMC = (TDirectoryFile*)fMC->Get("MyTask");
    TList* listMC = (TList*)dirMC->Get("ResultsList");

    // 创建画布
    TCanvas* c1 = new TCanvas("c_mc_dca_ptbins", "MC DCA_{xy} Distribution", 800, 600);
    gPad->SetLogy();
    
    // 创建图例
    TLegend* leg = new TLegend(0.75, 0.75, 0.95, 0.95);
    
    // 获取所有source的直方图
    std::vector<const char*> source_names = {
        Form("h_pt_dca_%s_origin", particle),
        Form("h_pt_dca_%s_lambdadecay", particle),
        Form("h_pt_dca_%s_otherdecay", particle),
        Form("h_pt_dca_%s_material", particle)
    };
    
    // 颜色数组
    int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1};
    const char* source_labels[] = {"Primary", "#Lambda decay", "Other decay", "Material", "Total Secondary"};
    
    // 为每个source创建投影
    bool first = true;
    TH1D* h_total_secondary = nullptr;
    
    // 先创建并累加所有的secondary sources
    for(size_t isrc = 1; isrc < source_names.size(); ++isrc) {  // 从1开始，跳过primary
        TH2D* h2 = (TH2D*)listMC->FindObject(source_names[isrc]);
        if(!h2) continue;
        
        // 获取pT积分的投影
        TString name_proj = Form("proj_%s", source_names[isrc]);
        TH1D* h_proj = h2->ProjectionY(name_proj);
        
        // 累加到总的secondary
        if(!h_total_secondary) {
            h_total_secondary = (TH1D*)h_proj->Clone("h_total_secondary");
        } else {
            h_total_secondary->Add(h_proj);
        }
        delete h_proj;  // 清理临时直方图
    }
    
    // 归一化总的secondary
    if(h_total_secondary && h_total_secondary->Integral() > 0) {
        h_total_secondary->Scale(1.0 / h_total_secondary->Integral());
    }
    
    // 绘制各个source
    for(size_t isrc = 0; isrc < source_names.size(); ++isrc) {
        TH2D* h2 = (TH2D*)listMC->FindObject(source_names[isrc]);
        if(!h2) continue;
        
        // 获取pT积分的投影
        TString name_proj = Form("proj_%s", source_names[isrc]);
        TH1D* h_proj = h2->ProjectionY(name_proj);
        
        // 归一化
        if(h_proj->Integral() > 0) h_proj->Scale(1.0 / h_proj->Integral());
        
        // 设置直方图属性
        h_proj->SetLineColor(colors[isrc]);
        h_proj->SetLineWidth(2);
        h_proj->SetMarkerStyle(20);
        h_proj->SetMarkerColor(colors[isrc]);
        
        // 设置坐标轴范围
        h_proj->GetXaxis()->SetRangeUser(0, 0.093);
        h_proj->GetYaxis()->SetRangeUser(1e-6, 9);
        
        // 绘制
        if(first) {
            h_proj->SetTitle(Form("%s %s MC Sources; DCA_{xy} (cm); Normalized Entries", 
                                dataset, particle));
            h_proj->Draw("E");
            first = false;
        } else {
            h_proj->Draw("E SAME");
        }
        
        // 添加到图例
        leg->AddEntry(h_proj, source_labels[isrc], "pl");
    }
    
    // 绘制总的secondary分布
    if(h_total_secondary) {
        h_total_secondary->SetLineColor(colors[4]);  // 使用橙色
        h_total_secondary->SetLineWidth(3);  // 加粗
        h_total_secondary->SetLineStyle(2);  // 虚线
        h_total_secondary->Draw("HIST SAME");
        leg->AddEntry(h_total_secondary, source_labels[4], "l");
    }
    
    // 绘制图例
    leg->Draw();
    
    // 保存图片
    c1->SaveAs(Form("mc_dca_%s_%s.pdf", dataset, particle));
    delete c1;
    if(h_total_secondary) delete h_total_secondary;
}